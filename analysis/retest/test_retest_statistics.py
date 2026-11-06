"""Analytic checks on independent errors and bounded continuation uncertainty."""
import math
from scipy.stats import beta, norm

import numpy as np
import pytest

from analysis.retest.statistics import (
    binomial_upper, classify, continuation_upper, contrast, mean_interval, normal_interval,
    residual, tail_diagnostics,
)


def test_mean_unbiased_sample_variance_and_family_quantile():
    r = mean_interval([1., 2., 3.], confidence=.95, family_size=4)
    assert r['estimate'] == 2.
    assert r['sample_variance'] == 1.
    assert r['variance_of_estimate'] == 1/3
    critical = norm.ppf(1-.05/(2*4))
    np.testing.assert_allclose(r['interval'], [2-critical/math.sqrt(3), 2+critical/math.sqrt(3)])


def test_residual_uses_a_squared_and_independent_variance_sum():
    # Deliberately different sample counts; a=.25 must scale slope variance by1/16.
    r = residual([-.2, .0, .2], [-2., 0., 2., 4.], .25, family_size=1)
    assert r['estimate'] == -.25
    np.testing.assert_allclose(r['variance_components']['validation_shift'], .04/3)
    np.testing.assert_allclose(r['variance_components']['response_prediction'], (1/16)*(20/3)/4)
    np.testing.assert_allclose(r['variance_of_estimate'], .04/3+(20/3)/64)


def test_contrast_has_no_artificial_paired_covariance():
    r = contrast([1., 2., 3.], [1., 2., 3.], family_size=1)
    assert r['estimate'] == 0
    assert r['variance_of_estimate'] == 2/3


def test_shared_prediction_can_be_used_without_independent_residual_claim():
    scores = np.array([-1., 0., 1., 2.])
    first = residual([-.1, .1], scores, .1)
    second = residual([-.2, .2], scores, .2)
    np.testing.assert_allclose(second['variance_components']['response_prediction'],
                               4*first['variance_components']['response_prediction'])


@pytest.mark.parametrize('interval,tolerance,expected', [
    ([-.1,.1],.1,'within_tolerance'),
    ([.1,.1],.1,'within_tolerance'),
    ([-.2,-.1],.1,'inconclusive'),
    ([.1,.2],.1,'inconclusive'),
    ([.10001,.2],.1,'outside_tolerance'),
    ([-.2,-.10001],.1,'outside_tolerance'),
    ([0,0],0,'within_tolerance'),
])
def test_classification_respects_closed_tolerance_band(interval,tolerance,expected):
    assert classify(interval,tolerance)==expected


def test_zero_observed_exits_do_not_imply_zero_population_bound():
    r = continuation_upper(np.zeros(65536),2000.,delta=.05/24)
    assert r['sample_mean_bound']==0
    assert r['nonzero_rows']==0
    assert r['upper_bound'] > .0002
    np.testing.assert_allclose(r['upper_bound'],2000*(-math.expm1(math.log(.05/24)/65536)))


def test_bound_scales_linearly_and_caps_at_deterministic_maximum():
    x=np.array([0.,.1,.2,.4]*100)
    a=continuation_upper(x,1.,delta=.01)
    b=continuation_upper(7*x,7.,delta=.01)
    np.testing.assert_allclose(b['upper_bound'],7*a['upper_bound'])
    short=continuation_upper([0.,0.],2.,delta=.01)
    assert 0 < short['upper_bound'] < 2.
    all_positive=continuation_upper([.1,.1],2.,delta=.01)
    assert all_positive['upper_bound']==2.


def test_smaller_failure_probability_cannot_reduce_upper_bound():
    x=np.array([0.,0.,.01]*100)
    assert continuation_upper(x,2.,delta=.001)['upper_bound'] > continuation_upper(x,2.,delta=.05)['upper_bound']


def test_tail_diagnostic_keeps_rare_score_and_constant_variance_undefined():
    r=tail_diagnostics(np.r_[np.zeros(99),10.])
    assert r['largest_one_percent_absolute_score_share']==1.
    np.testing.assert_allclose(r['largest_single_variance_share'],.99)
    assert tail_diagnostics([2.,2.])['largest_single_variance_share'] is None


@pytest.mark.parametrize('values', [[1.], [[1.,2.]], [1.,np.nan], [1.,np.inf]])
def test_invalid_samples_rejected(values):
    with pytest.raises(ValueError):mean_interval(values)


@pytest.mark.parametrize('confidence,family', [(1.,1), (0.,1), (True,1),(.95,0),(.95,True),(.95,2.5),(.95,10**100)])
def test_invalid_confidence_family_rejected(confidence,family):
    with pytest.raises(ValueError):mean_interval([0.,1.],confidence,family)


@pytest.mark.parametrize('a',[0.,-1.,1.01,np.nan,True,[.5]])
def test_invalid_mixture_fraction_rejected(a):
    with pytest.raises(ValueError):residual([0.,1.],[0.,1.],a)


@pytest.mark.parametrize('values,maximum,delta', [
    ([-.001,0.],2.,.05),([0.,2.001],2.,.05),([0.,1.],0.,.05),
    ([0.,1.],np.inf,.05),([0.,1.],True,.05),([0.,1.],2.,0.),
    ([0.,1.],2.,1.),([0.,1.],2.,True),
])
def test_invalid_continuation_domain_rejected(values,maximum,delta):
    with pytest.raises(ValueError):continuation_upper(values,maximum,delta=delta)


@pytest.mark.parametrize('interval,tol', [([1.,0.],1.),([0.,np.inf],1.),([0.],1.),([0.,1.],-.1)])
def test_invalid_classification_inputs_rejected(interval,tol):
    with pytest.raises(ValueError):classify(interval,tol)


def test_negative_or_nonfinite_estimate_variance_rejected():
    with pytest.raises(ValueError):normal_interval(0.,-.1)
    with pytest.raises(ValueError):normal_interval(np.inf,1.)


@pytest.mark.parametrize('n,k,delta', [
    (1,0,.05),(1,1,.05),(20,1,.01),(20,10,.025/24),(20,19,.001),
    (65536,0,.025/24),(65536,1,.025/24),(65536,3,.025/24),
    (65536,32768,.025/24),(65536,65535,.025/24),
    (8192,25,.025/24),(8192,8192,.025/24),
])
def test_clopper_pearson_against_independent_beta_quantile(n,k,delta):
    r=binomial_upper(k,n,delta=delta)
    expected=1. if k==n else beta.ppf(1-delta,k+1,n-k)
    assert r['upper_probability'] >= expected-2e-15
    np.testing.assert_allclose(r['upper_probability'],expected,rtol=0,atol=2e-11)


def test_count_bound_uses_rows_not_branch_counts_and_ignores_magnitude():
    low=continuation_upper([0.,0.,.001,0.],2000.,delta=.001)
    high=continuation_upper([0.,0.,1000.,0.],2000.,delta=.001)
    assert low['nonzero_rows']==1
    assert low['upper_bound']==high['upper_bound']


@pytest.mark.parametrize('n,k', [(0,0),(3,-1),(3,4),(3,True),(2.5,1)])
def test_invalid_binomial_counts_rejected(n,k):
    with pytest.raises(ValueError):binomial_upper(k,n,delta=.01)
