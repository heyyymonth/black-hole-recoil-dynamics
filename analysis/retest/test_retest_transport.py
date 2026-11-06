"""Complete-transition null control for the selected validation coupling."""
import numpy as np
from analysis.precision.coupling import make_transport
from analysis.precision.test_coupling import StateLaw
from analysis.response.run import flatten_result
from analysis.retest.run import jobs,read,PROTOCOL
from analysis.retest.verify import expected_jobs


def test_zero_mixture_symmetric_validation_pairs_are_exactly_equal():
    solver=make_transport('symmetric')
    initial=np.random.default_rng(771).normal(scale=.3,size=(256,6))
    result=solver.run_paired_mixture(StateLaw(),lambda x:-x,initial,1.03,.1,7712,
        mixture_probability=0.,spatial_stop=1.2,velocity_stop=1.4,region_radius=.4,
        density=lambda x:np.exp(-np.sum(x*x,axis=1)))
    arrays,_=flatten_result(result,'paired_mixture',1.03,.1)
    np.testing.assert_array_equal(arrays['paired_difference'],np.zeros(256))
    for field in ('states','residence','exposure','stop_cause','end_time'):
        np.testing.assert_array_equal(result['baseline'][field],result['mixture'][field])


def test_execution_and_independent_audit_agree_on_fixed_inventory():
    p=read(PROTOCOL)
    actual=list(jobs(p,'prediction'))+list(jobs(p,'validation'))
    assert actual==expected_jobs(p)
    assert sum(j['paths'] for j in actual)==368640
    assert len({tuple(j['seed']) for j in actual})==24
