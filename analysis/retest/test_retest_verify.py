"""Rejection controls for the independent response retest's evidence gates."""
from copy import deepcopy
import numpy as np
import pytest

from analysis.retest.verify import verify_chronology, verify_job, verify_record_arrays
from analysis.response.test_verify import fixture as weighted_fixture


def fixture():
    a, r = weighted_fixture()
    r['job']['seed'] = [2026110601, 20, 2, 0, 0]
    rng = np.random.default_rng(np.random.SeedSequence(r['job']['seed']).spawn(5)[4])
    a['replacement_steps'] = rng.integers(0, r['n_steps'], size=r['job']['paths'])
    return a, r


def test_valid_weighted_accounting_and_seeded_insertions_pass():
    a, r = fixture()
    verify_record_arrays(a, r)


def test_changed_seed_cannot_preserve_original_insertions():
    a, r = fixture()
    r['job']['seed'] = [2026110601, 20, 2, 0, 999]
    with pytest.raises(AssertionError):
        verify_record_arrays(a, r)


def test_changed_insertions_rejected_even_when_zero_scores_are_consistent():
    a, r = fixture()
    a['replacement_steps'][0] = 1-a['replacement_steps'][0]
    with pytest.raises(AssertionError):
        verify_record_arrays(a, r)


def job_fixture():
    expected = {'name':'prediction-example','stage':'prediction', 'method':'weighted',
        'paths':65536, 'seed':[2026110601,20,2,0,0], 'worker':0}
    identity = {'source_sha256':{'analysis/retest/run.py':'source'}}
    dependencies = {'execution_design_sha256':'design'}
    record = {'job':deepcopy(expected),'identity':deepcopy(identity),
              'stage_dependencies':deepcopy(dependencies)}
    return record, expected, identity, dependencies


@pytest.mark.parametrize('change', ['seed', 'count', 'worker', 'identity', 'dependency'])
def test_run_cannot_silently_change_registered_job(change):
    r, job, ids, deps = job_fixture()
    verify_job(r, job, ids, deps)
    if change=='seed':r['job']['seed'][-1] += 1
    elif change=='count':r['job']['paths'] //= 2
    elif change=='worker':r['job']['worker'] = 1
    elif change=='identity':r['identity']['source_sha256']['analysis/retest/run.py']='changed'
    else:r['stage_dependencies']['execution_design_sha256']='changed'
    with pytest.raises(AssertionError):
        verify_job(r, job, ids, deps)


def timeline():
    day = '2026-11-06T'
    def row(name, stage, worker, start, end):
        return {'job':{'name':name,'stage':stage,'worker':worker},
            'started_utc':day+start+'+00:00','completed_utc':day+end+'+00:00'}
    records = [row('p0','prediction',0,'01:00:00','02:00:00'),
               row('p1','prediction',1,'01:30:00','02:30:00'),
               row('v0','validation',0,'03:00:00','04:00:00'),
               row('v1','validation',1,'03:30:00','04:30:00')]
    design={'created_utc':day+'00:00:00+00:00'}
    register={'created_utc':day+'02:45:00+00:00'}
    return records, design, register


def test_independent_worker_overlap_is_allowed():
    verify_chronology(*timeline())


@pytest.mark.parametrize('change', ['early_register', 'late_register', 'late_design', 'overlap'])
def test_inconsistent_stage_or_worker_chronology_is_rejected(change):
    records, design, register = timeline()
    if change=='early_register':register['created_utc']='2026-11-06T02:15:00+00:00'
    elif change=='late_register':register['created_utc']='2026-11-06T03:15:00+00:00'
    elif change=='late_design':design['created_utc']='2026-11-06T01:15:00+00:00'
    else:records[1]['job']['worker']=0
    with pytest.raises(AssertionError):
        verify_chronology(records, design, register)
