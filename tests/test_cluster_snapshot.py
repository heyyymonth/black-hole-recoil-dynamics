"""Synthetic checks of snapshot units, population selection and compression."""
import gzip
import math
import numpy as np
import pytest

from recoil_dynamics.cluster_snapshot import (
    load_snapshot, read_conversion, compress_shells, cdf_distance, check_upstream_summary,
    empirical_mass_cdf_distance)


def make_fixture(tmp_path):
    units = tmp_path/'units.sh'
    units.write_text('lengthunitparsec=2\nlengthunitcgs=6e18\n'
                     'nbtimeunitcgs=3e12\ntimeunitsmyr=1000\n'
                     'echo THIS_MUST_NEVER_EXECUTE\n')
    header = '#2:m[MSUN] #3:r #4:vr #5:vt #8:binflag #9:m0[MSUN] #10:m1[MSUN] '
    header += '#13:a[AU] #14:e #15:startype #18:bin_startype0 #19:bin_startype1\n'
    data = np.zeros((4, 19))
    data[:, 1] = [1, 10, 25, 11]
    data[:, 2:5] = [[1, .1, .2], [2, -.2, .3], [3, .3, .4], [4, -.4, .5]]
    data[:, 7] = [-100, -100, 1, 1]
    data[:, 14] = [1, 14, -100, -100]
    data[2:, 8:10] = [[15, 10], [10, 1]]
    data[2:, 12:14] = [[.2, .3], [1, .5]]
    data[2:, 17:19] = [[14, 14], [14, 1]]
    snapshot = tmp_path/'snapshot.gz'
    with gzip.open(snapshot, 'wt') as stream:
        stream.write('# t=2 [code units]\n'+header)
        np.savetxt(stream, data)
    return snapshot, units


def test_units_and_companion_definitions(tmp_path):
    snap, units = make_fixture(tmp_path)
    s = load_snapshot(snap, units)
    assert s.age_gyr == 2
    np.testing.assert_allclose(s.radius, [2, 4, 6, 8])
    # Velocity = length/N-body time = 20km/s. Age time is not used here.
    np.testing.assert_allclose(s.vr, [2, -4, 6, -8])
    assert s.total_mass == 47  # Already MSUN; no extra code-mass factor.
    assert s.half_mass_radius == 6
    assert {k: int(v.sum()) for k, v in s.cohorts().items()} == {
        'all_systems': 4, 'all_binaries': 2, 'single_black_holes': 1,
        'bh_bearing_binaries': 2, 'binary_black_holes': 1}


def test_binary_mass_accounting_rejects_inconsistent_input(tmp_path):
    snap, units = make_fixture(tmp_path)
    with gzip.open(snap, 'rt') as stream:
        lines = stream.readlines()
    values = lines[4].split()
    values[1] = '99'
    lines[4] = ' '.join(values)+'\n'
    with gzip.open(snap, 'wt') as stream:
        stream.writelines(lines)
    with pytest.raises(ValueError, match='COM mass'):
        load_snapshot(snap, units)


@pytest.mark.parametrize('corruption', [None, 'count', 'mass', 'age', 'missing', 'duplicate', 'nonfinite'])
def test_upstream_accounting_check(tmp_path, corruption):
    snap, units = make_fixture(tmp_path)
    reference = tmp_path/'PARAM.txt'
    text = 'Mtot 47\nMbh 45\nage 2\nNsys 4\nNobj 6\nNbin 2\nN_14 4\nNbin_14_14 1\n'
    if corruption == 'count':
        text = text.replace('N_14 4', 'N_14 3')
    elif corruption == 'mass':
        text = text.replace('Mtot 47', 'Mtot 4700')
    elif corruption == 'age':
        text = text.replace('age 2', 'age 2000')
    elif corruption == 'missing':
        text = text.replace('Nbin_14_14 1\n', '')
    elif corruption == 'duplicate':
        text += 'Nsys 4\n'
    elif corruption == 'nonfinite':
        text = text.replace('Mbh 45', 'Mbh nan')
    reference.write_text(text)
    snapshot = load_snapshot(snap, units)
    if corruption is None:
        checks = check_upstream_summary(snapshot, reference)
        assert len(checks) == 8
        assert all(check['passed'] for check in checks)
    else:
        with pytest.raises(ValueError):
            check_upstream_summary(snapshot, reference)


def test_conversion_parser_rejects_duplicates_and_executable_values(tmp_path):
    _, units = make_fixture(tmp_path)
    units.write_text(units.read_text()+'nbtimeunitcgs=2e12\n')
    with pytest.raises(ValueError, match='ambiguous'):
        read_conversion(units)
    units.write_text('lengthunitparsec=$(echo 2)\n')
    with pytest.raises(ValueError):
        read_conversion(units)


def test_compression_preserves_total_mass_and_central_potential():
    radius = np.geomspace(.002, 200, 1001)
    mass = np.linspace(.01, 1, len(radius))
    for count in [1, 16, 128, 2048]:
        r, m = compress_shells(radius[::-1], mass[::-1], count)
        assert m.sum() == pytest.approx(mass.sum(), rel=3e-14)
        assert np.sum(m/r) == pytest.approx(np.sum(mass/radius), rel=3e-14)
        assert np.all(np.diff(r) > 0)
        assert len(r) <= count


def test_cdf_distance_and_empty_sample():
    assert cdf_distance([.25, .75], lambda x: x) == .25
    assert cdf_distance([], lambda x: x) is None
    assert cdf_distance([.5, .5], lambda x: x) == .5


def test_discrete_cdf_distances_include_simultaneous_jumps():
    assert empirical_mass_cdf_distance([1, 2], [1, 2], [1, 1]) == 0
    assert empirical_mass_cdf_distance([1, 2], [1, 2], [9, 1]) == pytest.approx(.4)
    assert empirical_mass_cdf_distance([1, 1, 3], [1, 2, 3], [2, 1, 1]) == pytest.approx(1/6)
    assert empirical_mass_cdf_distance([3], [1, 2], [1, 1]) == 1
    assert empirical_mass_cdf_distance([], [1, 2], [1, 1]) is None


@pytest.mark.parametrize('r,m,n', [([], [], 2), ([0], [1], 2),
    ([1], [-1], 2), ([1], [math.inf], 2), ([1], [1], 0), ([1], [1], 1.5)])
def test_invalid_compression(r, m, n):
    with pytest.raises(ValueError):
        compress_shells(r, m, n)
