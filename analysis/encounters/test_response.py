"""Independent finite-state checks of the single-insertion response."""

import numpy as np
import pytest

from analysis.encounters.response import (
    finite_state_response, randomized_replacement_response,
)


PD = np.array([[.8, .2], [.3, .7]])
PJ = np.array([[.4, .6], [.1, .9]])
P0 = np.array([.7, .3])
RD = np.array([.2, .8])
RJ = np.array([.3, .5])
TERMINAL = np.array([-.4, .3])


def direct_value(transitions, rewards, initial=P0, terminal=TERMINAL):
    """Forward distribution propagation, independent of backward derivative."""
    probability = initial.copy()
    observable = 0.
    for transition, reward in zip(transitions, rewards):
        observable += probability @ reward
        probability = probability @ transition
    return observable + probability @ terminal


def test_exact_response_equals_sum_of_single_replacements():
    n = 7
    expected = finite_state_response(PD, PJ, P0, RD, RJ, n, TERMINAL)
    baseline = direct_value([PD] * n, [RD] * n)
    differences = []
    for k in range(n):
        transitions, rewards = [PD] * n, [RD] * n
        transitions[k], rewards[k] = PJ, RJ
        differences.append(direct_value(transitions, rewards) - baseline)
    assert expected['baseline'] == pytest.approx(baseline, abs=1e-14)
    assert expected['derivative'] == pytest.approx(sum(differences), abs=1e-14)
    # Uniformly selecting an insertion and multiplying by n is the full sum.
    assert n * np.mean(differences) == pytest.approx(expected['derivative'])


def test_exact_response_is_derivative_not_endpoint_change():
    n = 11
    result = finite_state_response(PD, PJ, P0, RD, RJ, n, TERMINAL)
    baseline = direct_value([PD] * n, [RD] * n)
    h = 1e-5
    slopes = []
    for a in [h, h / 2]:
        value = direct_value([(1-a)*PD + a*PJ] * n,
                             [(1-a)*RD + a*RJ] * n)
        slopes.append((value-baseline)/a)
    richardson = 2*slopes[1] - slopes[0]
    assert richardson == pytest.approx(result['derivative'], abs=1e-8)
    endpoint_change = direct_value([PJ] * n, [RJ] * n) - baseline
    assert abs(endpoint_change - result['derivative']) > .1


def test_randomized_replacement_recovers_independent_exact_response():
    n, count = 8, 10000

    def callback(matrix, reward):
        def step(states, time, dt, rng):
            indices = states[:, 0].astype(int)
            draw = rng.random(len(states))
            new = (draw > matrix[indices, 0]).astype(float)[:, None]
            return new, reward[indices]
        return step

    initial = np.zeros((count, 1))
    # Deterministic initial state avoids uncertainty in the input law.
    expected = finite_state_response(PD, PJ, [1., 0.], RD, RJ, n, TERMINAL)
    sample = randomized_replacement_response(
        callback(PD, RD), callback(PJ, RJ), initial, n, .2, seed=80913,
        terminal_reward=lambda states, time: TERMINAL[states[:, 0].astype(int)])
    assert abs(sample['raw_mean']-expected['derivative']) < 6*sample['raw_standard_error']
    assert np.array_equal(sample['raw_response'],
                          n*(sample['replaced_observable']-sample['baseline_observable']))
    assert sample['normalized_mean'] == pytest.approx(sample['raw_mean']/(n*.2))


def test_state_dependent_draw_counts_do_not_change_other_paths():
    def variable_draws(states, time, dt, rng):
        assert states.shape == (1, 2)
        state = states.copy()
        # The first coordinate controls how many RNG values this path uses.
        draws = rng.normal(size=int(state[0, 0]))
        state[0, 1] += draws[-1]
        return state, np.array([dt*state[0, 1]])

    def replacement(states, time, dt, rng):
        state, reward = variable_draws(states, time, dt, rng)
        state[0, 1] += .4
        return state, reward

    initial = np.array([[1., 0.], [2., 0.], [3., 0.], [4., 0.]])
    kwargs = dict(n_steps=5, dt=.1, seed=43)
    before = randomized_replacement_response(variable_draws, replacement, initial, **kwargs)
    changed = initial.copy()
    changed[0, 0] = 200
    after = randomized_replacement_response(variable_draws, replacement, changed, **kwargs)
    # Altering row zero's draw count cannot shift any other path's noise,
    # including the paired post-insertion continuations used for uncertainty.
    for key in ['baseline_observable', 'replaced_observable', 'raw_response']:
        assert np.array_equal(before[key][1:], after[key][1:])
    null = randomized_replacement_response(variable_draws, variable_draws, initial, **kwargs)
    assert np.array_equal(null['raw_response'], np.zeros(4))


def test_absorbing_reward_and_input_validation():
    def absorbing(states, time, dt, rng):
        return states, dt*states[:, 0]

    initial = np.array([[0.], [1.], [0.], [1.]])
    result = randomized_replacement_response(absorbing, absorbing, initial, 5, .1, seed=2)
    assert np.array_equal(result['baseline_observable'], [.0, .5, .0, .5])
    assert np.array_equal(result['raw_response'], np.zeros(4))
    with pytest.raises(ValueError, match='positive integer'):
        randomized_replacement_response(absorbing, absorbing, initial, 0, .1, seed=2)
    with pytest.raises(ValueError, match='valid integer'):
        randomized_replacement_response(absorbing, absorbing, initial, 5, .1, seed=2,
                                        replacement_steps=[0, 1, 2, 5])
    with pytest.raises(ValueError, match='row stochastic'):
        finite_state_response(PD*2, PJ, P0, RD, RJ, 3)
