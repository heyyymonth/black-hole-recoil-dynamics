import json
import math

import numpy as np
import pytest

from analysis.response.statistics import (prediction_error,
    replacement_response_samples, summarize_response)


def test_replacement_weights_and_paired_slope_variance():
    baseline = np.array([.1, .2, .3, .4])
    replaced = baseline+np.array([0., .01, -.01, .02])
    samples = replacement_response_samples(baseline, replaced, 100)
    np.testing.assert_allclose(samples, [0, 1, -1, 2])
    result = summarize_response(samples)
    assert result['estimate'] == pytest.approx(.5)
    assert result['sample_variance'] == pytest.approx(5/3)
    assert result['standard_error'] == pytest.approx(math.sqrt(5/12))
    assert result['zero_contribution_count'] == 1


def test_independent_error_variance_has_all_three_contributions():
    baseline = np.array([.1, .2, .3, .4])
    mixture = np.array([.3, .4, .5, .6, .7])
    slope = np.array([-.3, .1, .2, .4, -.2, .5])
    a = .25
    result = prediction_error(mixture, baseline, slope, a)
    expected = mixture.mean()-baseline.mean()-a*slope.mean()
    variance = mixture.var(ddof=1)/5+baseline.var(ddof=1)/4+a*a*slope.var(ddof=1)/6
    assert result['estimate'] == pytest.approx(expected)
    assert result['variance_of_estimate'] == pytest.approx(variance)
    assert result['standard_error'] == pytest.approx(math.sqrt(variance))
    assert result['validation']['paired_sample_covariance'] is None


def test_paired_validation_uses_covariance_without_pairing_response_batch():
    baseline = np.array([.1, .2, .3, .4])
    mixture = baseline+.02
    slope = [-.1, .1, .3, .5]
    paired = prediction_error(mixture, baseline, slope, .1, paired_validation=True)
    independent = prediction_error(mixture, baseline, slope, .1)
    assert paired['observed_shift_standard_error'] < 1e-16
    assert paired['variance_of_estimate'] == pytest.approx(.1**2*np.var(slope, ddof=1)/4)
    assert paired['standard_error'] < independent['standard_error']
    assert paired['validation']['paired_sample_covariance'] == pytest.approx(baseline.var(ddof=1))


def test_bonferroni_reused_baseline_does_not_require_independent_comparisons():
    baseline = [.1, .2, .3, .4]
    slope = [.1, .2, .3, .4]
    mixture = [.11, .22, .33, .44]
    first = prediction_error(mixture, baseline, slope, .1, family_size=4)
    second = prediction_error(mixture, baseline, slope, .2, family_size=4)
    marginal = prediction_error(mixture, baseline, slope, .1)
    assert first['component_confidence'] == pytest.approx(1-.05/4)
    assert first['standard_error'] == marginal['standard_error']
    assert first['interval'][0] < marginal['interval'][0]
    assert first['interval'][1] > marginal['interval'][1]
    assert first['family_size'] == second['family_size'] == 4


def test_absolute_interval_mapping_and_zero_variance_caveat():
    result = prediction_error([.2]*3, [.1]*3, [.1]*4, 1.)
    assert result['estimate'] == pytest.approx(0, abs=1e-15)
    assert result['absolute_interval'][0] <= 1e-15
    large = prediction_error([.8]*3, [.1]*3, [.1]*4, 1.)
    assert large['absolute_interval'] == pytest.approx([.6, .6])
    summary = summarize_response([0.]*10)
    assert summary['zero_contribution_fraction'] == 1
    assert 'unseen' in summary['rare_contribution_caveat']


def test_json_serializable_and_no_invented_accuracy_verdict():
    result = prediction_error([.2, .4], [.1, .3], [.1, .2], .5,
                              family_size=8, paired_validation=True)
    json.dumps(result, allow_nan=False)
    json.dumps(summarize_response([0., 1.]), allow_nan=False)
    assert not any('verdict' in key or key == 'passed' for key in result)


@pytest.mark.parametrize('values', [[], [1.], [0., np.nan], [0., np.inf], [[0., 1.]]])
def test_invalid_response_sample(values):
    with pytest.raises(ValueError):
        summarize_response(values)


@pytest.mark.parametrize('family', [0, -1, 1.5, True])
def test_invalid_family_size(family):
    with pytest.raises(ValueError):
        summarize_response([0., 1.], family_size=family)


@pytest.mark.parametrize('a', [0, -1, 1.01, np.nan, np.inf])
def test_invalid_mixture_parameter(a):
    with pytest.raises(ValueError):
        prediction_error([.1, .2], [.1, .2], [0., 1.], a)


@pytest.mark.parametrize('confidence', [0, 1, np.nan])
def test_invalid_confidence(confidence):
    with pytest.raises(ValueError):
        summarize_response([0., 1.], confidence=confidence)


def test_invalid_reward_normalization_and_alignment():
    with pytest.raises(ValueError, match='normalized'):
        prediction_error([.1, 1.2], [.1, .2], [0., 1.], .5)
    with pytest.raises(ValueError, match='equal shape'):
        prediction_error([.1, .2, .3], [.1, .2], [0., 1.], .5, paired_validation=True)
    with pytest.raises(ValueError, match='aligned pairs'):
        replacement_response_samples([.1, .2], [.1, .2, .3], 2)
    with pytest.raises(ValueError, match='positive integer'):
        replacement_response_samples([.1, .2], [.1, .2], 1.5)
