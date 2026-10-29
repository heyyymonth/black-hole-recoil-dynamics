"""Independent finite-angle conservation, quadrature and sampling checks."""

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from physics import (DiscreteBath, PC_PER_KMS_MYR, closest_approach_annulus, impact_moments,
                     local_moments, scatter_impulse, velocity_ball)


def test_elastic_two_body_energy_and_momentum():
    rng = np.random.default_rng(8391)
    remnant, field, G = 60., 14., .004302
    v, u = rng.normal(size=(2, 128, 3))*8
    delta = scatter_impulse(v-u, field, remnant, G,
                             np.geomspace(1.e-6, 2., 128), rng.uniform(0, 2*np.pi, 128))
    v_after, u_after = v+delta, u-remnant/field*delta
    np.testing.assert_allclose(remnant*v_after+field*u_after, remnant*v+field*u, atol=4.e-13)
    initial = remnant*np.sum(v*v, axis=1)+field*np.sum(u*u, axis=1)
    final = remnant*np.sum(v_after*v_after, axis=1)+field*np.sum(u_after*u_after, axis=1)
    np.testing.assert_allclose(final, initial, rtol=3.e-15)
    # Scalar relative velocity also broadcasts against a vector of impacts.
    angles = np.array([0., 1., 2.])
    a = scatter_impulse([1., 2., 3.], field, remnant, G, [.001, .01, .1], angles)
    b = np.array([scatter_impulse([1., 2., 3.], field, remnant, G, impact, phi)
                  for impact, phi in zip([.001, .01, .1], angles)])
    np.testing.assert_allclose(a, b)


def test_analytic_moments_against_independent_b_azimuth_quadrature():
    w = np.array([1.3, -.7, 2.1])
    G, mass, remnant = .004302, 14., 60.
    x, weights = leggauss(180)
    phi = np.arange(256)*2*np.pi/256
    for lower, upper in [(0., .002), (.003, .08), (.2, 1.7), (.03, .030001)]:
        b = lower+(upper-lower)*(x+1)/2
        b, phases = np.meshgrid(b, phi, indexing="ij")
        impulses = scatter_impulse(np.broadcast_to(w, b.shape+(3,)), mass, remnant, G, b, phases)
        integration = np.linalg.norm(w)*b*(upper-lower)/2*weights[:, None]*(2*np.pi/len(phi))
        expected = impact_moments(w, mass, remnant, G, lower, upper)
        np.testing.assert_allclose(integration.sum(), expected["rate"], rtol=3.e-12)
        np.testing.assert_allclose(np.einsum("ij,ija->a", integration, impulses), expected["A"], rtol=3.e-10, atol=1.e-16)
        np.testing.assert_allclose(np.einsum("ij,ija,ijb->ab", integration, impulses, impulses), expected["D"], rtol=3.e-10, atol=1.e-16)


def test_zero_speed_empty_domain_and_slow_speed_limits():
    zero = impact_moments([0., 0., 0.], 14., 60., .004302, 0., .1)
    assert zero["rate"] == 0
    assert not np.any(zero["A"]) and not np.any(zero["D"])
    np.testing.assert_array_equal(scatter_impulse([0., 0., 0.], 14., 60., .004302, 0., 1.), np.zeros(3))
    empty = impact_moments([1., 2., 3.], 14., 60., .004302, .1, .1)
    assert empty["rate"] == 0 and not np.any(empty["D"])
    for speed in [1.e-8, 1.e-5, .01]:
        result = impact_moments([0., 0., speed], 14., 60., .004302, 0., .1)
        assert np.isfinite(result["D"]).all()
        assert np.linalg.eigvalsh(result["D"]).min() >= 0
        # At b/b90 -> 0 the deflection reverses relative velocity.
        if speed < 1.e-4:
            np.testing.assert_allclose(result["A"][2]/result["rate"], -2*14/74*speed, rtol=1.e-12)
            np.testing.assert_allclose(result["D"][2, 2]/result["rate"], (2*14/74*speed)**2, rtol=1.e-12)


def test_tensor_covariance_and_positive_diffusion():
    w = np.array([2.1, .7, -1.3])
    angle = .731
    R = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
    a = impact_moments(w, 14., 60., .004302, .001, .1)
    b = impact_moments(R@w, 14., 60., .004302, .001, .1)
    np.testing.assert_allclose(b["A"], R@a["A"], rtol=3.e-15)
    np.testing.assert_allclose(b["D"], R@a["D"]@R.T, rtol=5.e-15, atol=1.e-18)
    assert np.linalg.eigvalsh(a["D"]).min() > 0


def test_weak_deflection_limit_recovers_inverse_square_moments():
    w = np.array([2., -.5, 1.3])
    speed = np.linalg.norm(w)
    field, remnant, G = 14., 60., .004302
    a = G*(field+remnant)/speed**2
    result = impact_moments(w, field, remnant, G, 1000*a, 10000*a)
    logarithm = np.log(10.)
    expected_A = -4*np.pi*G**2*field*(field+remnant)*logarithm*w/speed**3
    unit = w/speed
    expected_D = 4*np.pi*G**2*field**2*logarithm/speed*(np.eye(3)-np.outer(unit, unit))
    np.testing.assert_allclose(result["A"], expected_A, rtol=3.e-7)
    np.testing.assert_allclose(result["D"], expected_D, rtol=1.e-6)


def test_quadrature_volume_and_isotropic_bath():
    cutoff, n = 5., 7.
    u, weights = velocity_ball(cutoff, 12)
    np.testing.assert_allclose(weights.sum(), 4*np.pi*cutoff**3/3, rtol=2.e-15)
    np.testing.assert_allclose(weights@u, 0., atol=5.e-13)
    number_df = lambda u: np.full((1, len(u)), n/(4*np.pi*cutoff**3/3))
    r = local_moments(number_df, [0, 0, 0], cutoff, [14.], 60., .004302, .001, .1, 20)
    np.testing.assert_allclose(r["A"], 0., atol=1.e-14)
    np.testing.assert_allclose(r["D"], np.eye(3)*np.trace(r["D"])/3, atol=2.e-15)
    # A uniform ball has mean speed 3U/4, giving an independent total rate.
    np.testing.assert_allclose(r["rate"], n*np.pi*(.1**2-.001**2)*3*cutoff/4/PC_PER_KMS_MYR, rtol=3.e-15)
    cold_nodes, cold_weights = velocity_ball(20., 36, radial_scale=.5)
    cold_density = np.exp(-np.sum(cold_nodes**2, axis=1)/(.5**2*2))/(2*np.pi*.5**2)**1.5
    np.testing.assert_allclose(cold_weights@cold_density, 1., rtol=2.e-8)


def test_sampled_impulses_have_matched_raw_moments():
    u = np.array([[0., 0., 0.], [1., -.5, .3], [-.2, .8, -.7]])
    number = np.array([[3., 1., 2.], [1., 2., .5]])
    bath = DiscreteBath(u, number, [7., 19.], 60., .004302, .001, .025)
    velocity = np.array([2., -.3, .6])
    expected = bath.moments(velocity)
    draws = bath.sample_events(velocity, np.random.default_rng(23191), 180000)
    mean, second = draws.mean(axis=0), np.einsum("na,nb->ab", draws, draws)/len(draws)
    mean_se = draws.std(axis=0, ddof=1)/np.sqrt(len(draws))
    products = draws[:, :, None]*draws[:, None, :]
    second_se = products.std(axis=0, ddof=1)/np.sqrt(len(draws))
    assert np.all(abs(mean-expected["A"]/expected["rate"]) < 6*mean_se)
    assert np.all(abs(second-expected["D"]/expected["rate"]) < 6*second_se)


def test_invalid_domain_is_rejected_and_no_event_for_empty_bath():
    with pytest.raises(ValueError):
        impact_moments([0., 0., 1.], 14., 60., .004302, .2, .1)
    bath = DiscreteBath(np.zeros((1, 3)), np.zeros((1, 1)), [14.], 60., .004302, 0., .1)
    assert bath.moments([0., 0., 1.])["rate"] == 0
    with pytest.raises(ValueError, match="zero-rate"):
        bath.sample_events([0., 0., 1.], np.random.default_rng(1))


def test_closest_approach_exclusion_and_batch_consistency():
    speeds = np.array([0., 1.e-6, 1., 10.])
    lo, hi, empty = closest_approach_annulus(speeds, 14., 60., .004302, .001, .03)
    assert np.all(empty[:2]) and not empty[3]
    assert np.all(lo[empty] == hi[empty])
    for k in np.flatnonzero(~empty):
        np.testing.assert_allclose(lo[k]**2, .001**2+2*.004302*74*.001/speeds[k]**2)
    # The same speed-dependent bounds are integrated for batches and events.
    lower = lambda w, m: closest_approach_annulus(w, m, 60., .004302, .0001, .1)[0]
    u = np.array([[0., 0., 0.], [1., -.5, .3], [-.2, .8, -.7]])
    bath = DiscreteBath(u, np.array([[3., 1., 2.], [1., 2., .5]]),
                        [7., 19.], 60., .004302, lower, .1)
    velocities = np.array([[0., 0., 0.], [1., 2., 3.], [-1., -.2, .4]])
    batched = bath.moments_many(velocities, batch_size=2)
    for k, v in enumerate(velocities):
        scalar = bath.moments(v)
        for key in scalar:
            np.testing.assert_allclose(batched[key][k], scalar[key], rtol=5.e-15, atol=1.e-16)
    np.testing.assert_allclose(bath.rates_many(velocities), batched["rate"], rtol=5.e-15)
