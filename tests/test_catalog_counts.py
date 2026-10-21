"""Synthetic accounting and age-selection checks; no real catalog outcomes."""
import math

import pytest

from recoil_dynamics.catalog_counts import (
    BINARY_COUNT_KEYS, TYPE_COUNT_KEYS, load_param, parse_param, summarize_models)


def fixture_values():
    values = {key: 0 for key in TYPE_COUNT_KEYS + BINARY_COUNT_KEYS}
    values.update({'N_0': 10, 'N_1': 7, 'N_14': 5,
                   'Nbin_1_0': 1, 'Nbin_14_0': 1, 'Nbin_14_14': 1,
                   'Nsys': 19, 'Nobj': 22, 'Nbin': 3, 'Mtot': 100.0, 'age': 12.0})
    return values


def fixture_text(values=None):
    values = fixture_values() if values is None else values
    return '\n'.join(f'{key} {value}' for key, value in values.items()) + '\n'


def parse_fixture(values=None, suffix='0', model='model-a'):
    return parse_param(fixture_text(values), model=model, path=f'{model}/PARAM{suffix}.txt')


def test_population_accounting_and_source_identity():
    row = parse_fixture()
    assert row == {'model': 'model-a', 'path': 'model-a/PARAM0.txt', 'age_gyr': 12.0,
                   'system_count': 19, 'object_count': 22, 'binary_count': 3,
                   'bh_count': 5, 'bbh_count': 1, 'bh_nonbbh_binary_count': 1,
                   'single_bh_count': 2, 'host_mass_msun': 100.0}


def test_zero_black_holes_are_real_zeros():
    values = fixture_values()
    values.update(N_14=0, Nbin_14_0=0, Nbin_14_14=0, Nobj=17, Nbin=1, Nsys=16)
    row = parse_fixture(values)
    assert (row['bh_count'], row['bbh_count'], row['single_bh_count']) == (0, 0, 0)
    assert summarize_models([row])[0]['selected_bbh_count'] == 0
    assert summarize_models([row])[0]['selected_meets_count_screen'] is False


def test_exact_scientific_notation_and_undefined_unused_dispersion(tmp_path):
    values = fixture_values()
    values.update(Nsys='1.9e1', N_14='5.0000000000000000', sig_14='nan')
    path = tmp_path/'PARAM.txt'
    path.write_text('# synthetic\n' + fixture_text(values))
    assert load_param(path, model='a', source_path='a/PARAM.txt')['bh_count'] == 5


@pytest.mark.parametrize('key', ['N_14', 'Nbin_14_14', 'Nsys', 'Nbin_15_15'])
@pytest.mark.parametrize('bad', ['nan', 'inf', '-1', '1.5', '1.00000000000000001', 'bad'])
def test_all_required_counts_reject_invalid_values(key, bad):
    values = fixture_values()
    values[key] = bad
    with pytest.raises(ValueError, match='[Cc]ount'):
        parse_fixture(values)


@pytest.mark.parametrize('key', ['Mtot', 'age'])
@pytest.mark.parametrize('bad', ['nan', 'inf', '-1', '0', 'bad'])
def test_mass_and_age_must_be_finite_positive(key, bad):
    values = fixture_values()
    values[key] = bad
    with pytest.raises(ValueError, match='[Qq]uantity'):
        parse_fixture(values)


@pytest.mark.parametrize('key', ['N_0', 'N_15', 'Nbin_15_15', 'Nobj', 'age', 'Mtot'])
def test_missing_fields_are_not_imputed(key):
    values = fixture_values()
    del values[key]
    with pytest.raises(ValueError, match='Missing'):
        parse_fixture(values)


@pytest.mark.parametrize('duplicate', ['N_14 5', 'age 12', 'sig_14 nan\nsig_14 nan'])
def test_duplicate_keys_are_rejected(duplicate):
    with pytest.raises(ValueError, match='Duplicate'):
        parse_param(fixture_text() + duplicate, model='a', path='a/PARAM.txt')


@pytest.mark.parametrize('key,value,match', [
    ('N_0', 11, 'Stellar-type'),
    ('Nbin_0_0', 1, 'Binary-type'),
    ('Nsys', 20, 'Nsys'),
])
def test_global_count_identities(key, value, match):
    values = fixture_values()
    values[key] = value
    with pytest.raises(ValueError, match=match):
        parse_fixture(values)


def test_bh_membership_cannot_exceed_bh_objects():
    values = fixture_values()
    values.update(N_14=2, N_0=13)  # Preserve global counts, break BH membership.
    with pytest.raises(ValueError, match='BH binary membership'):
        parse_fixture(values)


def test_bh_massless_companion_is_included_in_other_binaries():
    values = fixture_values()
    values.update(Nbin_14_0=0, Nbin_15_14=1, N_15=1, N_0=9)
    assert parse_fixture(values)['bh_nonbbh_binary_count'] == 1


@pytest.mark.parametrize('extra', ['N_16 0', 'Nbin_0_14 0', 'Nbin_16_0 0'])
def test_unsupported_or_reversed_count_keys_fail(extra):
    with pytest.raises(ValueError, match='Unsupported'):
        parse_param(fixture_text() + extra, model='a', path='a/PARAM.txt')


def test_unknown_count_is_still_validated():
    with pytest.raises(ValueError, match='[Cc]ount'):
        parse_param(fixture_text() + 'Nmsp nan\n', model='a', path='a/PARAM.txt')


def selection_row(model, path, age, bbh):
    return {'model': model, 'path': path, 'age_gyr': age, 'bbh_count': bbh}


def test_common_age_selection_ties_use_earlier_age_then_path():
    rows = [selection_row('a', 'z', 12.125, 100),
            selection_row('a', 'b', 11.875, 4),
            selection_row('a', 'a', 11.875, 5),
            selection_row('a', 'early', 10.5, 0)]
    result = summarize_models(rows)[0]
    assert result['selected_path'] == 'a'
    assert result['selected_age_gyr'] == 11.875
    assert result['selected_age_offset_gyr'] == -0.125
    assert result['selected_bbh_count'] == 5
    assert result['snapshot_count'] == 4
    assert (result['minimum_age_gyr'], result['maximum_age_gyr']) == (10.5, 12.125)
    assert (result['minimum_bbh_count'], result['median_bbh_count'],
            result['maximum_bbh_count']) == (0, 4.5, 100)


def test_missing_coverage_is_none_not_zero_and_screen_is_inclusive():
    rows = [selection_row('missing', 'missing', 11.0, 100),
            selection_row('zero', 'zero', 12.0, 0),
            selection_row('screen', 'boundary', 12.25, 30)]
    results = {row['model']: row for row in summarize_models(rows)}
    missing = results['missing']
    for field in ('selected_path', 'selected_age_gyr', 'selected_bbh_count',
                  'selected_age_offset_gyr', 'selected_meets_count_screen'):
        assert missing[field] is None
    assert results['zero']['selected_bbh_count'] == 0
    assert results['screen']['selected_meets_count_screen'] is True
    assert results['screen']['selected_path'] == 'boundary'
    assert results['screen']['selected_age_offset_gyr'] == 0.25


def test_saved_state_count_screen_does_not_count_distinct_binaries():
    rows = [selection_row('a', 'a1', 11.0, 29),
            selection_row('a', 'a2', 11.5, 30),
            selection_row('a', 'a3', 12.0, 31),
            selection_row('b', 'b1', 12.0, 0)]
    results = {row['model']: row for row in summarize_models(rows)}
    assert results['a']['saved_states_meeting_count_screen'] == 2
    assert results['b']['saved_states_meeting_count_screen'] == 0
    assert summarize_models(rows, minimum_bbh_count=31)[0][
        'saved_states_meeting_count_screen'] == 1


def test_duplicate_source_records_cannot_inflate_history():
    row = selection_row('a', 'same', 12, 4)
    with pytest.raises(ValueError, match='Duplicate source'):
        summarize_models([row, row])


def test_sorting_is_deterministic_and_no_interpolation():
    rows = [selection_row('b', 'b', 12.4, 30), selection_row('a', 'a', 11.6, 2)]
    assert summarize_models(rows) == summarize_models(rows[::-1])
    assert [row['model'] for row in summarize_models(rows)] == ['a', 'b']
    assert all(row['selected_path'] is None for row in summarize_models(rows))


@pytest.mark.parametrize('kwargs', [
    {'reference_age_gyr': math.nan}, {'reference_age_gyr': 0},
    {'tolerance_gyr': -0.1}, {'tolerance_gyr': math.inf},
    {'minimum_bbh_count': 0}, {'minimum_bbh_count': 30.5}, {'minimum_bbh_count': True},
])
def test_invalid_selection_settings(kwargs):
    with pytest.raises(ValueError):
        summarize_models([], **kwargs)
