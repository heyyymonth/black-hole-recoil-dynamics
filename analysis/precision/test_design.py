"""Insertion selection and independent-stream accounting controls."""
import numpy as np
from analysis.precision.run import read, PROTOCOL, jobs, selections


def test_fixed_inventory_has_disjoint_simulation_and_selection_streams():
    p=read(PROTOCOL); items=list(jobs(p))
    assert len(items)==36 and len({x['name'] for x in items})==36
    seeds=[tuple(x[k]) for x in items for k in ('seed','selection_seed')]
    assert len(set(seeds))==72
    assert {tuple(x['seed']) for x in items}.isdisjoint(
        {tuple(x['selection_seed']) for x in items})


def test_explicit_strata_cover_all_scheduled_steps_with_correct_allocation():
    p=read(PROTOCOL)
    for job in jobs(p):
        steps,strata=selections(p,job)
        if job['method']['selection']=='uniform':
            assert steps is None and (strata==-1).all()
        else:
            np.testing.assert_array_equal(np.bincount(strata),np.full(4,256))
            for label,(lo,hi) in enumerate(zip(p['step_strata'][:-1],p['step_strata'][1:])):
                assert ((steps[strata==label]>=lo)&(steps[strata==label]<hi)).all()
            again,_=selections(p,job)
            np.testing.assert_array_equal(steps,again)
