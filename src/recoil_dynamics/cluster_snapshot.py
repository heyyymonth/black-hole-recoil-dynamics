"""Read an explicitly pinned CMC text snapshot and audit its target population.

No upstream code or unit-conversion shell script is executed. Masses labeled
MSUN and separations labeled AU are already physical. Velocities use N-body
time units; snapshot age uses the distinct CMC relaxation-time units.
"""
from dataclasses import dataclass
import gzip
import math
from pathlib import Path
import re

import numpy as np
from scipy.special import betainc


@dataclass
class ClusterSnapshot:
    age_gyr: float
    mass: np.ndarray
    radius: np.ndarray
    vr: np.ndarray
    vt: np.ndarray
    binary: np.ndarray
    m0: np.ndarray
    m1: np.ndarray
    semimajor: np.ndarray
    eccentricity: np.ndarray
    kind: np.ndarray
    kind0: np.ndarray
    kind1: np.ndarray

    @property
    def speed(self):
        return np.hypot(self.vr, self.vt)

    @property
    def total_mass(self):
        return float(self.mass.sum())

    @property
    def half_mass_radius(self):
        order = np.argsort(self.radius)
        index = np.searchsorted(np.cumsum(self.mass[order]), self.total_mass/2)
        return float(self.radius[order[index]])

    def cohorts(self):
        return {
            'all_systems': np.ones(len(self.mass), dtype=bool),
            'all_binaries': self.binary,
            'single_black_holes': ~self.binary & (self.kind == 14),
            'bh_bearing_binaries': self.binary & ((self.kind0 == 14) | (self.kind1 == 14)),
            'binary_black_holes': self.binary & (self.kind0 == 14) & (self.kind1 == 14),
        }


def read_conversion(path):
    text = Path(path).read_text()
    needed = ['lengthunitparsec', 'lengthunitcgs', 'nbtimeunitcgs', 'timeunitsmyr']
    result = {}
    for key in needed:
        matches = re.findall(r'^'+key+r'=([+\-0-9.eE]+)\s*$', text, re.MULTILINE)
        if len(matches) != 1:
            raise ValueError(f'Missing or ambiguous conversion: {key}')
        value = float(matches[0])
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'Invalid conversion: {key}')
        result[key] = value
    return result


def load_snapshot(snapshot_path, conversion_path):
    units = read_conversion(conversion_path)
    with gzip.open(snapshot_path, 'rt') as stream:
        header = stream.readline()
        columns = stream.readline()
    time = re.search(r'^# t=([+\-0-9.eE]+)', header)
    required = ['#2:m[MSUN]', '#3:r ', '#4:vr ', '#5:vt ', '#8:binflag',
                '#9:m0[MSUN]', '#10:m1[MSUN]', '#13:a[AU]', '#14:e ',
                '#15:startype', '#18:bin_startype0', '#19:bin_startype1']
    if time is None or not all(token in columns for token in required):
        raise ValueError('Unsupported CMC snapshot header or column convention')
    age = float(time.group(1))*units['timeunitsmyr']/1000
    if not math.isfinite(age) or age < 0:
        raise ValueError('Invalid snapshot age')
    # Explicit indices refer to the checked numbered upstream header.
    raw = np.loadtxt(snapshot_path, usecols=[1, 2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18], ndmin=2)
    # The supplied text snapshot uses -100 for singles; the toolkit defines
    # binaries by flag==1, unlike newer outputs which may use flag==0 for singles.
    if not np.isfinite(raw).all() or not np.isin(raw[:, 4], [-100, 0, 1]).all():
        raise ValueError('Nonfinite data or invalid binary flags')
    if np.any(raw[:, :2] <= 0) or np.any(raw[:, 3] < 0):
        raise ValueError('Masses/radii must be positive and tangential speeds nonnegative')
    binary = raw[:, 4] == 1
    if not np.allclose(raw[binary, 0], raw[binary, 5]+raw[binary, 6], rtol=2e-6, atol=1e-6):
        raise ValueError('Binary COM mass differs from sum of components')
    if np.any(raw[binary, 5:8] <= 0) or np.any((raw[binary, 8] < 0) | (raw[binary, 8] >= 1)):
        raise ValueError('Invalid binary component masses, separation or eccentricity')
    velocity_unit = units['lengthunitcgs']/units['nbtimeunitcgs']/1e5
    return ClusterSnapshot(age, raw[:, 0], raw[:, 1]*units['lengthunitparsec'],
                           raw[:, 2]*velocity_unit, raw[:, 3]*velocity_unit,
                           binary, raw[:, 5], raw[:, 6], raw[:, 7], raw[:, 8],
                           raw[:, 9], raw[:, 10], raw[:, 11])


def check_upstream_summary(snapshot, path):
    """Cross-check parsing against the authors' separately tabulated summary.

    This checks the same simulation state, not independent astrophysical data.
    Rounded binary COM masses require a small relative mass tolerance.
    """
    s = snapshot
    single_bh = ~s.binary & (s.kind == 14)
    component0_bh = s.binary & (s.kind0 == 14)
    component1_bh = s.binary & (s.kind1 == 14)
    actual = {
        'Mtot': s.total_mass,
        'Mbh': float(s.mass[single_bh].sum()+s.m0[component0_bh].sum()
                     +s.m1[component1_bh].sum()),
        'age': s.age_gyr,
        'Nsys': len(s.mass),
        'Nobj': len(s.mass)+int(s.binary.sum()),
        'Nbin': int(s.binary.sum()),
        'N_14': int(single_bh.sum()+component0_bh.sum()+component1_bh.sum()),
        'Nbin_14_14': int((component0_bh & component1_bh).sum()),
    }
    text = Path(path).read_text()
    checks = []
    for key, value in actual.items():
        matches = re.findall(r'^'+key+r'\s+([^\s]+)\s*$', text, re.MULTILINE)
        if len(matches) != 1:
            raise ValueError(f'Missing or ambiguous upstream summary: {key}')
        reference = float(matches[0])
        rtol = 0 if key.startswith('N') else (1e-7 if key.startswith('M') else 1e-12)
        if not math.isfinite(reference) or not math.isclose(value, reference, rel_tol=rtol, abs_tol=0):
            raise ValueError(f'Upstream summary mismatch: {key}: {value} vs {reference}')
        checks.append({'quantity': key, 'derived': value, 'reference': reference,
                       'relative_tolerance': rtol, 'passed': True})
    return checks


def cdf_distance(values, reference_cdf):
    """ECDF distance from a continuous reference, not a significance test."""
    values = np.sort(np.asarray(values))
    if len(values) == 0:
        return None
    reference = np.asarray(reference_cdf(values))
    n = len(values)
    return float(max(np.max(np.arange(1, n+1)/n-reference),
                     np.max(reference-np.arange(n)/n)))


def empirical_mass_cdf_distance(values, reference_radii, reference_masses):
    """Exact descriptive distance between two step CDFs, including tied jumps."""
    values = np.sort(np.asarray(values))
    if len(values) == 0:
        return None
    order = np.argsort(reference_radii)
    radii = np.asarray(reference_radii)[order]
    mass = np.asarray(reference_masses)[order]
    cumulative = np.r_[0., np.cumsum(mass)/mass.sum()]
    jumps = np.union1d(values, radii)
    target = np.searchsorted(values, jumps, side='right')/len(values)
    reference = cumulative[np.searchsorted(radii, jumps, side='right')]
    return float(np.max(np.abs(target-reference)))


def compress_shells(radius, mass, count):
    """Approximately equal-mass bins; harmonic radii preserve central Phi."""
    radius, mass = np.asarray(radius), np.asarray(mass)
    if radius.ndim != 1 or mass.shape != radius.shape or len(radius) == 0:
        raise ValueError('Expected nonempty one-dimensional radius/mass arrays')
    if (not np.isfinite(radius).all() or not np.isfinite(mass).all()
            or np.any(radius <= 0) or np.any(mass <= 0)
            or int(count) != count or count < 1):
        raise ValueError('Invalid shell compression inputs')
    order = np.argsort(radius)
    r, m = radius[order], mass[order]
    midpoint_mass = np.cumsum(m)-m/2
    groups = np.minimum((midpoint_mass/m.sum()*count).astype(int), count-1)
    shell_mass = np.bincount(groups, weights=m, minlength=count)
    inverse_sum = np.bincount(groups, weights=m/r, minlength=count)
    used = shell_mass > 0
    return shell_mass[used]/inverse_sum[used], shell_mass[used]


def population_audit(snapshot, gravitational_constant):
    mass = snapshot.total_mass
    rh = snapshot.half_mass_radius
    b = rh*math.sqrt(2**(2/3)-1)
    velocity = math.sqrt(gravitational_constant*mass/b)
    def plummer_cdf(r):
        x = r/b
        return x**3/(1+x*x)**1.5

    def speed_cdf(y):
        return betainc(1.5, 4.5, np.clip(y, 0, 1)**2)

    rows = []
    for name, mask in snapshot.cohorts().items():
        count = int(mask.sum())
        r, m = snapshot.radius[mask], snapshot.mass[mask]
        escape = np.sqrt(2)*velocity*(1+(r/b)**2)**(-.25)
        y = snapshot.speed[mask]/escape
        row = {'population': name, 'count': count,
               'median_radius_pc': float(np.median(r)) if count else None,
               'median_radius_over_host_rh': float(np.median(r)/rh) if count else None,
               'median_system_mass_msun': float(np.median(m)) if count else None,
               'minimum_system_mass_msun': float(np.min(m)) if count else None,
               'maximum_system_mass_msun': float(np.max(m)) if count else None,
               'radial_cdf_distance_from_host_mass': empirical_mass_cdf_distance(
                   r, snapshot.radius, snapshot.mass),
               'radial_cdf_distance_from_plummer': cdf_distance(r, plummer_cdf),
               'speed_fraction_cdf_distance_from_plummer': cdf_distance(y, speed_cdf),
               'mean_squared_speed_fraction': float(np.mean(y*y)) if count else None,
               'count_above_matched_plummer_escape': int(np.sum(y > 1)),
               'median_binary_semimajor_axis_au': None,
               'minimum_binary_semimajor_axis_au': None,
               'maximum_binary_semimajor_axis_au': None}
        if count and np.all(snapshot.binary[mask]):
            a = snapshot.semimajor[mask]
            row.update(median_binary_semimajor_axis_au=float(np.median(a)),
                       minimum_binary_semimajor_axis_au=float(np.min(a)),
                       maximum_binary_semimajor_axis_au=float(np.max(a)))
        rows.append(row)
    scales = {'age_gyr': snapshot.age_gyr, 'total_mass_msun': mass,
              'half_mass_radius_pc': rh, 'matched_plummer_scale_pc': b,
              'velocity_unit_kms': velocity, 'matched_plummer_escape_kms': math.sqrt(2)*velocity}
    return scales, rows
