"""Reject changed insertion designs even when rewards remain self-consistent."""
import numpy as np
import pytest

from analysis.precision.verify import verify_record_arrays
from analysis.response.test_verify import fixture as weighted_fixture


def fixture(selection):
    arrays, record = weighted_fixture()
    seed = [2026110501, 40, 0, 0, 0]
    selection_seed = [2026110501, 41, 0, 0, 0]
    protocol = {'step_strata': [0, 1, 2, 3, 4]}
    record['n_steps'] = 4
    record['job'].update({'method': {'selection': selection},
                          'seed': seed, 'selection_seed': selection_seed})
    for branch in record['branches']:
        arrays[branch+'__evaluated_steps'][:] = 4
    if selection == 'uniform':
        rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(5)[4])
        arrays['replacement_steps'] = rng.integers(0, 4, size=4)
        arrays['insertion_stratum'] = np.full(4, -1, dtype=np.int64)
    else:
        arrays['replacement_steps'] = np.arange(4, dtype=np.int64)
        arrays['insertion_stratum'] = np.arange(4, dtype=np.int64)
    record['array_inventory'] = {
        key: {'shape': list(value.shape), 'dtype': str(value.dtype)}
        for key, value in arrays.items()}
    return arrays, record, protocol


@pytest.mark.parametrize('selection', ['uniform', 'stratified'])
def test_correct_retained_insertion_design_passes(selection):
    arrays, record, protocol = fixture(selection)
    verify_record_arrays(arrays, record, protocol)


@pytest.mark.parametrize('selection', ['uniform', 'stratified'])
def test_changed_selected_step_is_rejected_despite_unchanged_zero_scores(selection):
    arrays, record, protocol = fixture(selection)
    arrays['replacement_steps'][0] = (arrays['replacement_steps'][0]+1) % 4
    with pytest.raises(AssertionError):
        verify_record_arrays(arrays, record, protocol)


def test_permuted_stratum_labels_cannot_reassign_observations():
    arrays, record, protocol = fixture('stratified')
    arrays['insertion_stratum'] = arrays['insertion_stratum'][::-1].copy()
    with pytest.raises(AssertionError):
        verify_record_arrays(arrays, record, protocol)


def test_uniform_selection_cannot_be_relabelled_as_a_stratum():
    arrays, record, protocol = fixture('uniform')
    arrays['insertion_stratum'][0] = 0
    with pytest.raises(AssertionError, match='uniform rows have stratum labels'):
        verify_record_arrays(arrays, record, protocol)


def test_active_insertion_cannot_claim_it_was_never_evaluated():
    arrays, record, protocol = fixture('stratified')
    arrays['diffusion_plus__evaluated_steps'][3] = 3
    with pytest.raises(AssertionError, match='active insertion was never evaluated'):
        verify_record_arrays(arrays, record, protocol)
