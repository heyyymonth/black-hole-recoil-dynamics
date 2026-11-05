"""Marginal-law, private-namespace and complete-transition precision controls."""
from math import erf, exp, sqrt

import numpy as np
import pytest

from analysis.matched.transport import _gaussian_increment as original_increment
from analysis.response import transport as original
from analysis.precision.coupling import (
    make_transport, symmetric_gaussian_increment, symmetric_root,
)


class StateLaw:
    """Compound Gaussian impulses with exact first/raw-second moments."""
    mean = np.array([.14, -.07, .11])
    sigma = np.array([.3, .2, .4])

    def moments(self, states):
        rate = 2*(1+np.linalg.norm(states[:, :3], axis=1))
        raw_second = np.diag(self.sigma**2)+np.outer(self.mean, self.mean)
        return {'rate': rate, 'excluded': .3+np.abs(states[:, 0]),
                'A': rate[:, None]*self.mean,
                'D': rate[:, None, None]*raw_second}

    def sample_events(self, states, rng):
        return self.mean+self.sigma*rng.standard_normal((len(states), 3))


def compare_nested(first, second):
    assert first.keys() == second.keys()
    for name in first:
        if isinstance(first[name], dict):
            compare_nested(first[name], second[name])
        elif isinstance(first[name], np.ndarray):
            np.testing.assert_array_equal(first[name], second[name], err_msg=name)
        else:
            assert first[name] == second[name], name


def test_eigen_private_solver_bitwise_matches_frozen_full_weighted_run():
    states = np.random.default_rng(31).normal(scale=.28, size=(512, 6))
    states[0, 0] = 1.2
    states[1, 3] = 2.
    settings = dict(spatial_stop=1.2, velocity_stop=1.4, region_radius=.4,
                    density=lambda x: np.exp(-np.sum(x*x, axis=1)))
    old = original.run_weighted_response(StateLaw(), lambda x: -x, states,
                                         1.03, .1, [20261105, 12], **settings)
    new = make_transport('eigen').run_weighted_response(
        StateLaw(), lambda x: -x, states, 1.03, .1, [20261105, 12], **settings)
    compare_nested(old, new)
    stops = set(np.concatenate([old[name]['stop_cause'] for name in
        ('jump_zero', 'jump_positive', 'diffusion_plus', 'diffusion_minus')]))
    assert {'horizon', 'spatial_domain', 'velocity_domain', 'excluded_channel'} <= stops
    assert np.any(old['positive_branch_count'] > 0)
    assert np.any(old['diffusion_plus']['exposure'] > 0)


def test_private_namespaces_do_not_mutate_original_or_other_arm():
    old_advance = original.advance_step
    old_gaussian = original._gaussian_increment
    eigen = make_transport('eigen')
    symmetric = make_transport('symmetric')
    assert original.advance_step is old_advance
    assert original._gaussian_increment is old_gaussian is original_increment
    assert eigen.advance_step.__globals__ is eigen.__dict__
    assert symmetric.advance_step.__globals__ is symmetric.__dict__
    assert eigen._gaussian_increment is old_gaussian
    assert symmetric._gaussian_increment is symmetric_gaussian_increment
    symmetric._gaussian_increment = None
    assert eigen._gaussian_increment is old_gaussian
    assert original._gaussian_increment is old_gaussian
    assert make_transport('symmetric')._gaussian_increment is symmetric_gaussian_increment


@pytest.mark.parametrize('D', [
    np.diag([0., 0., 0.]), np.diag([0., .3, 1.2]),
    np.array([[2., .3, .1], [.3, 1., .2], [.1, .2, 3.]]),
])
def test_exact_normal_cubature_preserves_drift_and_covariance(D):
    # Six nodes integrate every first and second standard-normal moment exactly.
    noise = np.concatenate((sqrt(3)*np.eye(3), -sqrt(3)*np.eye(3)))
    n, dt = len(noise), .137
    drift = np.array([.2, -.3, .1])
    args = (np.tile(drift, (n, 1)), np.tile(D, (n, 1, 1)),
            np.full(n, dt), noise)
    for increment in (original_increment, symmetric_gaussian_increment):
        values = increment(*args)
        np.testing.assert_allclose(values.mean(axis=0), dt*drift, atol=1e-15)
        centered = values-dt*drift
        np.testing.assert_allclose(centered.T@centered/n, dt*D,
                                   atol=1e-14, rtol=1e-13)
    values = symmetric_gaussian_increment(*args)
    opposite = symmetric_gaussian_increment(*args[:-1], -noise)
    np.testing.assert_allclose((values+opposite)/2,
                               np.tile(dt*drift, (n, 1)), atol=1e-15)


def test_repeated_eigenvalue_basis_and_signs_do_not_change_symmetric_root(monkeypatch):
    # Rotation within a repeated-eigenvalue subspace and arbitrary column signs
    # represent the same covariance; its unique PSD root cannot depend on them.
    angle = .73
    Q = np.array([[np.cos(angle), -np.sin(angle), 0.],
                  [np.sin(angle), np.cos(angle), 0.], [0., 0., -1.]])
    base = np.array([[1., 2., 3.], [-2., 3., 1.], [3., -1., 2.]])
    orientation, _ = np.linalg.qr(base)
    D = (orientation @ np.diag([.4, .4, 2.]) @ orientation.T)[None]
    expected = symmetric_root(D)
    eigh = np.linalg.eigh

    def changed_basis(matrix):
        values, vectors = eigh(matrix)
        # The mathematical covariance has a repeated pair; rounding may split
        # that pair at machine precision, immaterial to this tolerance check.
        return values, vectors@Q

    with monkeypatch.context() as patch:
        patch.setattr(np.linalg, 'eigh', changed_basis)
        changed = symmetric_root(D)
    np.testing.assert_allclose(changed, expected, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(expected@expected, D, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(expected, expected.transpose(0, 2, 1), atol=1e-15)


def test_covariance_crossing_has_continuous_symmetric_root():
    distances = []
    for eps in (1e-3, 1e-6, 1e-9):
        pair = np.array([np.diag([1-eps, 1+eps, 2.]),
                         np.diag([1+eps, 1-eps, 2.])])
        roots = symmetric_root(pair)
        distance = np.linalg.norm(roots[1]-roots[0])
        distances.append(distance)
        assert distance < 2*eps
    assert distances[2] < distances[1] < distances[0]


def test_roundoff_psd_clipping_matches_frozen_tolerance():
    D = np.diag([-5e-13, .5, 1.])[None]
    root = symmetric_root(D)
    np.testing.assert_array_equal(root[0], np.diag([0., sqrt(.5), 1.]))
    args = (np.zeros((1, 3)), D, np.ones(1), np.ones((1, 3)))
    assert np.isfinite(original_increment(*args)).all()
    assert np.isfinite(symmetric_gaussian_increment(*args)).all()
    invalid = np.diag([-2e-12, .5, 1.])[None]
    for increment in (original_increment, symmetric_gaussian_increment):
        with pytest.raises(ValueError, match='positive semidefinite'):
            increment(args[0], invalid, args[2], args[3])


@pytest.mark.parametrize('D', [np.eye(3), np.zeros((2, 2, 2)),
    np.full((1, 3, 3), np.nan), np.array([[[1., .1, 0.], [0., 1., 0.], [0., 0., 1.]]])])
def test_invalid_tensor_rejected(D):
    with pytest.raises(ValueError):
        symmetric_root(D)


@pytest.mark.parametrize('field', ['negative_time', 'bad_noise', 'nan_drift'])
def test_invalid_increment_inputs_rejected(field):
    A, D = np.zeros((2, 3)), np.tile(np.eye(3), (2, 1, 1))
    time, normal = np.ones(2), np.ones((2, 3))
    if field == 'negative_time': time[0] = -1.
    elif field == 'bad_noise': normal = np.ones((2, 2))
    else: A[0, 0] = np.nan
    with pytest.raises(ValueError):
        symmetric_gaussian_increment(A, D, time, normal)


def test_unknown_factor_rejected():
    with pytest.raises(ValueError, match='factor'):
        make_transport('unregistered')


def test_symmetric_complete_stopped_response_matches_exact_one_step_expectation():
    # A single deterministic x impulse allows an analytical Poisson/Gaussian
    # stop comparison, while each trajectory includes its half-step residence.
    class Impulse:
        rate, kick = .8, .7

        def moments(self, states):
            n = len(states)
            A, D = np.zeros((n, 3)), np.zeros((n, 3, 3))
            A[:, 0] = self.rate*self.kick
            D[:, 0, 0] = self.rate*self.kick**2
            return {'rate': np.full(n, self.rate), 'excluded': np.zeros(n), 'A': A, 'D': D}

        def sample_events(self, states, rng):
            out = np.zeros((len(states), 3)); out[:, 0] = self.kick
            return out

    law, n, limit = Impulse(), 30000, 1.15
    result = make_transport('symmetric').run_weighted_response(
        law, np.zeros_like, np.zeros((n, 6)), 1., 1., [20261105, 91],
        spatial_stop=100., velocity_stop=limit, region_radius=10.)
    cdf = lambda x: .5*(1+erf(x/sqrt(2)))
    mean, sigma = law.rate*law.kick, sqrt(law.rate)*law.kick
    jump_stop = 1-exp(-law.rate)*(1+law.rate)
    gaussian_stop = cdf((-limit-mean)/sigma)+1-cdf((limit-mean)/sigma)
    expected = -.5*(jump_stop-gaussian_stop)
    score = (result['jump_zero_weight']*result['jump_zero']['normalized_residence']
             +result['jump_positive_weight']*result['jump_positive']['normalized_residence']
             -.5*(result['diffusion_plus']['normalized_residence']
                  +result['diffusion_minus']['normalized_residence']))
    assert abs(score.mean()-expected) < 6*score.std(ddof=1)/sqrt(n)
    for name in ('jump_zero', 'jump_positive', 'diffusion_plus', 'diffusion_minus'):
        stopped = result[name]['stop_cause'] == 'velocity_domain'
        np.testing.assert_array_equal(result[name]['normalized_residence'][stopped], .5)
        np.testing.assert_array_equal(result[name]['normalized_residence'][~stopped], 1.)
