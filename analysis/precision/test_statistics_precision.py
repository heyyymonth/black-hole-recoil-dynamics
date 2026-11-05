"""Independent design and uncertainty controls for precision-pilot statistics."""
import numpy as np
import pytest

from analysis.precision.statistics import compare, efficiency_bootstrap, estimate


def test_uniform_estimate_matches_hand_calculation():
    result = estimate([1., 2., 3., 4.])
    assert result['mean'] == 2.5
    assert result['variance_of_mean'] == pytest.approx(5/12)
    assert result['equivalent_score_variance'] == pytest.approx(5/3)



def test_absolute_tail_diagnostics_keep_both_score_signs():
    result = estimate([-10., -2., 0., 2., 10.])
    assert result['absolute_score_quantiles']['0.5'] == 2.
    assert result['absolute_score_quantiles']['0.9'] == 10.
    assert result['zero_count'] == 1
    assert result['largest_one_percent_absolute_score_share'] == pytest.approx(10/24)


def test_fixed_strata_mean_does_not_use_empirical_allocation_as_probability():
    # Strata have sizes 2, 3, 4, 5, but each has population weight 1/4.
    groups = [np.array([-1., 1.]), np.array([8., 10., 12.]),
              np.array([17., 19., 21., 23.]), np.array([26., 28., 30., 32., 34.])]
    values = np.concatenate(groups)
    labels = np.repeat(np.arange(4), [2, 3, 4, 5])
    result = estimate(values, labels)
    assert result['mean'] == 15.
    # The four sample variances are 2, 4, 20/3, 10.
    expected = (2/2+4/3+(20/3)/4+10/5)/16
    assert result['variance_of_mean'] == pytest.approx(expected)
    assert result['mean'] != pytest.approx(values.mean())


def test_stratification_eliminates_only_between_stratum_sampling_noise():
    values = np.repeat([0., 10., -20., 30.], 8)
    labels = np.repeat(np.arange(4), 8)
    result = estimate(values, labels)
    assert result['variance_of_mean'] == 0.
    assert result['mean'] == 5.
    assert estimate(values)['variance_of_mean'] > 0


def test_independent_difference_adds_variances_and_preserves_sign():
    result = compare({'mean': 5., 'variance_of_mean': 4.},
                     {'mean': 2., 'variance_of_mean': 5.})
    assert result['difference'] == 3.
    assert result['variance_of_difference'] == 9.
    assert result['standard_error'] == 3.
    assert result['interval'][0] < 0 < result['interval'][1]
    # Increased family size widens the component interval.
    adjusted = compare({'mean': 5., 'variance_of_mean': 4.},
                       {'mean': 2., 'variance_of_mean': 5.}, family_size=6)
    assert adjusted['interval'][0] < result['interval'][0]
    assert adjusted['interval'][1] > result['interval'][1]


def test_bootstrap_is_seeded_scale_invariant_and_cpu_cost_is_explicit():
    a = np.arange(32., dtype=float)
    b = np.arange(32., dtype=float)**1.1
    kwargs = dict(reps=128, seed=31, cost_a=20., cost_b=10.)
    first = efficiency_bootstrap(a, None, b, None, **kwargs)
    assert first == efficiency_bootstrap(a, None, b, None, **kwargs)
    scaled = efficiency_bootstrap(-7*a, None, -7*b, None, **kwargs)
    assert scaled['variance_ratio']['estimate'] == pytest.approx(first['variance_ratio']['estimate'])
    np.testing.assert_allclose(scaled['variance_ratio']['interval'], first['variance_ratio']['interval'])
    assert first['cost_adjusted']['estimate'] == pytest.approx(2*first['variance_ratio']['estimate'])
    np.testing.assert_allclose(first['cost_adjusted']['interval'],
                              2*np.array(first['variance_ratio']['interval']))
    assert first['invalid_bootstrap_denominators'] == 0


def test_bootstrap_resamples_within_strata_without_mixing_their_means():
    labels = np.repeat(np.arange(4), 8)
    a = np.tile(np.arange(8.), 4)
    b = a+.3*np.sin(a)
    offsets = np.repeat([0., 100., -200., 1000.], 8)
    base = efficiency_bootstrap(a, labels, b, labels, reps=128, seed=53)
    shifted = efficiency_bootstrap(a+offsets, labels, b+offsets, labels, reps=128, seed=53)
    assert shifted['variance_ratio']['estimate'] == pytest.approx(base['variance_ratio']['estimate'])
    np.testing.assert_allclose(shifted['variance_ratio']['interval'], base['variance_ratio']['interval'])


def test_observed_zero_denominator_is_explicitly_unavailable():
    result = efficiency_bootstrap(np.arange(8.), None, np.ones(8), None, reps=32)
    assert result['variance_ratio'] == {'estimate': None, 'interval': None}
    assert result['unavailable_reason'] == 'Observed denominator variance is zero'


def test_invalid_bootstrap_replicates_are_not_silently_dropped():
    # With two observations, an iid bootstrap has positive mass at variance zero.
    result = efficiency_bootstrap([0., 2.], None, [0., 1.], None,
                                  reps=128, seed=2, cost_a=1., cost_b=2.)
    assert result['variance_ratio']['estimate'] == 4.
    assert result['invalid_bootstrap_denominators'] > 0
    assert result['variance_ratio']['interval'] is None
    assert result['cost_adjusted']['interval'] is None


@pytest.mark.parametrize('scores,strata', [
    ([], None), ([1.], None), ([1., np.nan], None),
    (np.arange(8.), np.repeat(np.arange(4), 2).astype(float)),
    (np.arange(8.), np.zeros(8, dtype=int)),
    (np.arange(7.), np.array([0, 1, 1, 2, 2, 3, 3])),
])
def test_invalid_sampling_designs_rejected(scores, strata):
    with pytest.raises(ValueError):
        estimate(scores, strata)


def test_efficiency_requires_matching_counts_and_both_positive_costs():
    with pytest.raises(ValueError, match='equal total'):
        efficiency_bootstrap([0., 1.], None, [0., 1., 2.], None)
    with pytest.raises(ValueError, match='both total CPU'):
        efficiency_bootstrap([0., 1.], None, [0., 1.], None, cost_a=1.)
    with pytest.raises(ValueError, match='positive'):
        efficiency_bootstrap([0., 1.], None, [0., 1.], None, cost_a=1., cost_b=0.)
