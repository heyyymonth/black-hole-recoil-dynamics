"""Independent moment controls for the complete-transition response identity.

Constant-rate deterministic impulses have known compound-Poisson cumulants.
Polynomial terminal rewards isolate the stochastic identity without relying on
cluster profiles or residence geometry. Degree six exhibits a nonzero finite-a
remainder even when the response derivative is implemented correctly.
"""
from math import comb, erf, exp, sqrt

import numpy as np
import pytest

from analysis.response.transport import run_mixture, run_weighted_response


class ConstantImpulse:
    def __init__(self, rate, impulse):
        self.rate, self.impulse = rate, impulse

    def moments(self, states):
        n = len(states)
        drift = np.zeros((n, 3)); drift[:, 0] = self.rate*self.impulse
        raw_second = np.zeros((n, 3, 3)); raw_second[:, 0, 0] = self.rate*self.impulse**2
        return {'rate': np.full(n, self.rate), 'excluded': np.zeros(n),
                'A': drift, 'D': raw_second}

    def sample_events(self, states, rng):
        increments = np.zeros((len(states), 3)); increments[:, 0] = self.impulse
        return increments


def moments_from_cumulants(cumulants):
    """Bell recurrence for raw moments, independent of transport code."""
    order = len(cumulants)
    moments = np.zeros(order+1); moments[0] = 1
    for n in range(1, order+1):
        moments[n] = sum(comb(n-1, j-1)*cumulants[j-1]*moments[n-j]
                         for j in range(1, n+1))
    return moments


def moments_of_sum(first, second):
    return np.array([sum(comb(n, j)*first[j]*second[n-j] for j in range(n+1))
                     for n in range(len(first))])


def exact_mixture_moment(order, steps, per_step_mean_count, impulse, initial, a):
    jump = moments_from_cumulants([per_step_mean_count*impulse**j
                                  for j in range(1, order+1)])
    gaussian = moments_from_cumulants([per_step_mean_count*impulse,
                                      per_step_mean_count*impulse**2]+[0.]*(order-2))
    step_moments = (1-a)*gaussian+a*jump
    total = np.array([initial**j for j in range(order+1)])
    for _ in range(steps):
        total = moments_of_sum(total, step_moments)
    return total[order]


def exact_insertion_derivative(order, steps, lam, impulse, initial):
    # Total mean/variance are unchanged by replacing exactly one transition.
    common = [initial+steps*lam*impulse, steps*lam*impulse**2]
    baseline = moments_from_cumulants(common+[0.]*(order-2))[order]
    one_replacement = moments_from_cumulants(
        common+[lam*impulse**j for j in range(3, order+1)])[order]
    return steps*(one_replacement-baseline)


@pytest.mark.parametrize('order', [4, 6])
def test_multistep_weighted_response_matches_exact_polynomial_derivative(order):
    n, steps, horizon, rate, impulse, initial = 50000, 4, 1., 1.2, .35, .2
    dt = horizon/steps; lam = rate*dt
    state = np.zeros((n, 6)); state[:, 3] = initial
    result = run_weighted_response(ConstantImpulse(rate, impulse), np.zeros_like,
        state, horizon, dt, [2026110201, order], spatial_stop=1e6, velocity_stop=1e6)
    branch = {key: result[key]['states'][:, 3]**order for key in
              ('jump_zero', 'jump_positive', 'diffusion_plus', 'diffusion_minus')}
    sample = steps*(result['jump_zero_weight']*branch['jump_zero']
                    +result['jump_positive_weight']*branch['jump_positive']
                    -.5*(branch['diffusion_plus']+branch['diffusion_minus']))
    expected = exact_insertion_derivative(order, steps, lam, impulse, initial)
    se = sample.std(ddof=1)/np.sqrt(n)
    assert abs(sample.mean()-expected) < 6*se+1e-12
    for key in branch:
        assert result[key]['alive'].all()
    baseline = exact_mixture_moment(order, steps, lam, impulse, initial, 0.)
    endpoint = exact_mixture_moment(order, steps, lam, impulse, initial, 1.)
    if order == 4:
        assert endpoint-baseline == pytest.approx(expected, abs=1e-12)
    else:
        remainder = 10*steps*(steps-1)*lam**2*impulse**6
        assert endpoint-baseline-expected == pytest.approx(remainder, abs=1e-12)
        assert remainder > 0


def test_sixth_reward_exact_mixture_has_quadratic_remainder():
    steps, lam, impulse, initial = 4, .3, .35, .2
    baseline = exact_mixture_moment(6, steps, lam, impulse, initial, 0.)
    derivative = exact_insertion_derivative(6, steps, lam, impulse, initial)
    quadratic = 10*steps*(steps-1)*lam**2*impulse**6
    for a in [.001, .05, .25, .5, 1.]:
        actual = exact_mixture_moment(6, steps, lam, impulse, initial, a)
        assert actual-baseline-a*derivative == pytest.approx(quadratic*a*a, abs=1e-12)


def test_multistep_intermediate_mixture_matches_exact_sixth_moment():
    n, steps, horizon, rate, impulse, initial, a = 80000, 4, 1., 1.2, .35, .2, .5
    dt = horizon/steps
    states = np.zeros((n, 6)); states[:, 3] = initial
    result = run_mixture(ConstantImpulse(rate, impulse), np.zeros_like,
        states, horizon, dt, 2026110202, a, spatial_stop=1e6, velocity_stop=1e6)
    sample = result['states'][:, 3]**6
    exact = exact_mixture_moment(6, steps, rate*dt, impulse, initial, a)
    assert abs(sample.mean()-exact) < 6*sample.std(ddof=1)/np.sqrt(n)
    assert result['alive'].all()


def test_weighted_response_includes_exact_midstep_stopping_reward():
    # With one frozen step, every path resides for half a horizon before the
    # impulse. A velocity-domain stop then removes the remaining half. All
    # survivors stay inside the generously sized residence region, so the
    # response is -1/2 times the jump-minus-Gaussian stopping probability.
    n, rate, impulse, limit = 50000, .8, .7, 1.15
    states = np.zeros((n,6))
    result = run_weighted_response(ConstantImpulse(rate,impulse),np.zeros_like,
        states,1.,1.,2026110203,spatial_stop=100.,velocity_stop=limit,
        region_radius=10.)
    def cdf(x):
        return .5*(1+erf(x/sqrt(2)))
    jump_stop = 1-exp(-rate)*(1+rate)  # N >= 2 since 1*c < limit < 2*c.
    mean, sigma = rate*impulse,sqrt(rate)*impulse
    gaussian_stop = cdf((-limit-mean)/sigma)+1-cdf((limit-mean)/sigma)
    expected = -.5*(jump_stop-gaussian_stop)
    reward = {name:result[name]['normalized_residence'] for name in
              ('jump_zero','jump_positive','diffusion_plus','diffusion_minus')}
    sample = (result['jump_zero_weight']*reward['jump_zero']
              +result['jump_positive_weight']*reward['jump_positive']
              -.5*(reward['diffusion_plus']+reward['diffusion_minus']))
    assert abs(sample.mean()-expected) < 6*sample.std(ddof=1)/sqrt(n)
    for name,values in reward.items():
        stopped = result[name]['stop_cause']=='velocity_domain'
        np.testing.assert_array_equal(values[stopped],.5)
        np.testing.assert_array_equal(values[~stopped],1.)
