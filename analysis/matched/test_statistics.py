import json
import math

import numpy as np
import pytest

from analysis.matched.statistics import (compare_samples, empirical_cdf_distance,
                                        mean_summary, summarize_samples,
                                        wilson_probability)


def test_mean_and_independent_contrast_use_unbiased_variance():
    first = np.array([0., 1., 2., 3.])
    second = first+10
    result = compare_samples(first, second)
    assert result['first']['sample_variance'] == pytest.approx(5/3)
    assert result['first']['standard_error'] == pytest.approx(math.sqrt(5/12))
    difference = result['mean_difference']
    assert difference['estimate'] == -10
    assert difference['standard_error'] == pytest.approx(math.sqrt(5/6))
    # Equal-index differences are constant; a paired SE would incorrectly vanish.
    assert difference['standard_error'] > 0
    half = 1.959963984540054*math.sqrt(5/6)
    assert difference['interval'] == pytest.approx([-10-half, -10+half])
    assert difference['absolute_interval'] == pytest.approx([10-half, 10+half])


def test_absolute_interval_crossing_zero_and_constant_samples():
    result = compare_samples([0., 1.], [0., 1.])['mean_difference']
    assert result['absolute_interval'][0] == 0
    assert result['absolute_interval'][1] == pytest.approx(result['interval'][1])
    constant = compare_samples([.2]*5, [.7]*5)['mean_difference']
    assert constant['standard_error'] == 0
    assert constant['absolute_interval'] == pytest.approx([.5, .5])


def test_cdf_sup_with_ties_and_unequal_sample_sizes():
    # At zero the CDFs are 2/3 and 1/4; at one they are 1 and 1/2.
    result = empirical_cdf_distance([0., 0., 1.], [0., 1., 2., 2.])
    assert result['estimate'] == .5
    assert empirical_cdf_distance([0., 0., 1.], [0., 0., 1.])['estimate'] == 0
    assert empirical_cdf_distance([0., 0.], [1., 1.])['estimate'] == 1


def test_union_dkw_allocation_and_clipping():
    result = empirical_cdf_distance(np.zeros(1000), np.ones(2000), .9)
    eps_a = math.sqrt(math.log(40)/2000)
    eps_b = math.sqrt(math.log(40)/4000)
    assert result['epsilon_first'] == pytest.approx(eps_a)
    assert result['epsilon_second'] == pytest.approx(eps_b)
    assert result['interval'] == pytest.approx([1-eps_a-eps_b, 1])
    equal = empirical_cdf_distance(np.zeros(1000), np.zeros(2000), .9)
    assert equal['interval'] == pytest.approx([0, eps_a+eps_b])


def test_wilson_zero_and_all_events_have_nontrivial_limits():
    zero, all_events = wilson_probability(0, 100), wilson_probability(100, 100)
    z2 = 1.959963984540054**2
    assert zero['interval'][0] == pytest.approx(0, abs=1e-16)
    assert zero['interval'][1] == pytest.approx(z2/(100+z2))
    assert all_events['interval'][0] == pytest.approx(100/(100+z2))
    assert all_events['interval'][1] == pytest.approx(1)
    assert wilson_probability(50, 100)['interval'] == pytest.approx([.4038315303659956, .5961684696340044])


def fixture_arrays():
    return dict(residence=[.3, .2, .1, .4], exposure=[.5, .3, .2, .6],
                stop_time=[10., 4., 2., 10.],
                stop_cause=['horizon', 'spatial_domain', 'excluded_channel', 'horizon'],
                events=[1, 0, 2, 1], horizon=10.)


def test_stopped_summary_remaining_time_bound_and_json_serialization():
    result = summarize_samples(**fixture_arrays())
    missing = result['numerical_domain_missing_residence']
    # Only numerical-domain stop contributes; excluded-channel kill is modeled.
    assert missing['remaining_time_bound']['mean'] == pytest.approx(.6/4)
    assert missing['remaining_time_bound']['mean'] <= missing['stop_probability']['estimate']
    assert missing['probability_based_mean_upper_bound'] >= missing['stop_probability']['estimate']
    assert result['stop_probabilities']['velocity_domain']['count'] == 0
    assert result['stop_probabilities']['velocity_domain']['interval'][1] > 0
    assert result['event_counts']['total'] == 4
    json.dumps(result, allow_nan=False)
    json.dumps(compare_samples([0., 1.], [1., 2.]), allow_nan=False)


def test_no_observed_exit_does_not_imply_zero_missing_probability():
    args = fixture_arrays()
    args['stop_time'] = [10.]*4
    args['stop_cause'] = ['horizon']*4
    result = summarize_samples(**args)['numerical_domain_missing_residence']
    assert result['remaining_time_bound']['mean'] == 0
    assert result['probability_based_mean_upper_bound'] > 0


def test_remaining_time_bound_respects_endpoint_roundoff():
    args = fixture_arrays()
    args['stop_time'][1] = -1e-12
    args['residence'][1] = 0.
    missing = summarize_samples(**args)['numerical_domain_missing_residence']
    assert missing['remaining_time_bound']['mean'] == missing['stop_probability']['estimate']


@pytest.mark.parametrize('sample', [[], [1.], [0., np.nan], [0., np.inf], [[1., 2.]]])
def test_invalid_numeric_samples(sample):
    with pytest.raises(ValueError):
        mean_summary(sample)
    with pytest.raises(ValueError):
        compare_samples([0., 1.], sample)


@pytest.mark.parametrize('confidence', [0, 1, -1, np.nan, np.inf])
def test_invalid_confidence(confidence):
    with pytest.raises(ValueError):
        compare_samples([0., 1.], [0., 1.], confidence)


@pytest.mark.parametrize('changes', [
    {'horizon': 0}, {'stop_time': [9., 4., 2., 10.]},
    {'stop_time': [10., 4., -2., 10.]}, {'residence': [.3, .5, .1, .4]},
    {'exposure': [.5, -.3, .2, .6]}, {'events': [1., 1.5, 2., 1.]},
    {'events': [1., -1., 2., 1.]}, {'stop_cause': ['horizon', '', 'excluded_channel', 'horizon']},
    {'stop_cause': ['horizon', 'escape', 'excluded_channel', 'horizon']},
    {'events': [1., 0.]}, {'numerical_stop_causes': ('unknown',)},
])
def test_invalid_path_accounting(changes):
    args = fixture_arrays()
    args.update(changes)
    with pytest.raises(ValueError):
        summarize_samples(**args)


@pytest.mark.parametrize('count,n', [(0, 0), (-1, 3), (4, 3), (.5, 3), (True, 3)])
def test_invalid_binomial_counts(count, n):
    with pytest.raises(ValueError):
        wilson_probability(count, n)
