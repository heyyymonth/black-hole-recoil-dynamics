"""Independent analytical limits for collisionless trajectory occupation."""
import numpy as np
import pytest

from analysis.residence.orbits import integrate_orbit, host_acceleration, host_scales


def test_harmonic_radial_crossings_and_residence():
    orbit = integrate_orbit(lambda x: -x, [2., 0., 0.], [0., 0., 0.],
                            2*np.pi, stop_radius=3., max_step=.03)
    expected = np.array([np.pi/3, 2*np.pi/3, 4*np.pi/3, 5*np.pi/3])
    np.testing.assert_allclose([e['time'] for e in orbit.crossings], expected,
                               atol=3e-10, rtol=0)
    assert [e['direction'] for e in orbit.crossings] == ['entry','exit','entry','exit']
    assert orbit.residence == pytest.approx(2*np.pi/3, abs=1e-9)
    t = np.linspace(0., 2*np.pi, 101)
    y = orbit.solution(t)
    np.testing.assert_allclose(y[0], 2*np.cos(t), atol=2e-10)
    np.testing.assert_allclose(y[3], -2*np.sin(t), atol=2e-10)
    assert orbit.stop_cause == 'horizon'
    assert len(orbit.turning_times) == 3


def test_kepler_circular_orbit_conserves_energy_and_angular_momentum():
    orbit = integrate_orbit(lambda x: -x/np.linalg.norm(x)**3,
                            [1.,0.,0.], [0.,1.,0.], 4*np.pi,
                            region_radius=1.1, stop_radius=3., max_step=.04)
    y = orbit.solution(np.linspace(0., orbit.end_time, 201))
    x, v = y[:3].T, y[3:].T
    np.testing.assert_allclose(.5*np.sum(v*v,axis=1)-1/np.linalg.norm(x,axis=1),
                               -.5, atol=2e-10, rtol=0)
    np.testing.assert_allclose(np.cross(x,v), np.tile([0.,0.,1.],(len(x),1)),
                               atol=2e-10, rtol=0)
    assert orbit.residence == pytest.approx(4*np.pi)
    assert not orbit.crossings
    assert not orbit.turning_times


@pytest.mark.parametrize('horizon', [2*np.pi-1e-8, np.nextafter(2*np.pi,0.),
                                    2*np.pi, np.nextafter(2*np.pi,np.inf), 2*np.pi+1e-8])
def test_region_tangencies_are_not_crossings(horizon):
    orbit = integrate_orbit(lambda x: -x, [1.,0.,0.], [0.,0.,0.],
                            horizon, stop_radius=2., max_step=.03)
    assert not orbit.crossings
    assert orbit.residence == pytest.approx(horizon, abs=1e-9)
    assert any(abs(t-np.pi) < 1e-9 for t in orbit.region_touches)
    assert all(abs(t/np.pi-round(t/np.pi)) < 1e-9 for t in orbit.region_touches)


def test_outer_grazing_orbit_never_enters_region():
    orbit = integrate_orbit(lambda x: -x, [1.,0.,0.], [0.,2.,0.],
                            2*np.pi, stop_radius=3., max_step=.03)
    assert not orbit.crossings
    assert orbit.residence == 0.


def test_stationary_and_circular_boundary_have_zero_strict_residence():
    stationary = integrate_orbit(lambda x: np.zeros(3), [1.,0.,0.], [0.,0.,0.],
                                 2., stop_radius=3.)
    circular = integrate_orbit(lambda x: -x, [1.,0.,0.], [0.,1.,0.],
                               2*np.pi, stop_radius=3., max_step=.03)
    for orbit in (stationary, circular):
        assert orbit.residence == 0
        assert not orbit.crossings
        assert not orbit.turning_times


def test_stop_is_validity_censoring_and_initial_outside_stops_immediately():
    orbit = integrate_orbit(lambda x: np.zeros(3), [0.,0.,0.], [1.,0.,0.],
                            5., stop_radius=2., max_step=.1)
    assert orbit.end_time == pytest.approx(2., abs=1e-10)
    assert orbit.residence == pytest.approx(1., abs=1e-10)
    assert orbit.stop_cause == 'validity_boundary'
    for radius in (2.,3.):
        early = integrate_orbit(lambda x: np.zeros(3), [radius,0.,0.], [0.,0.,0.],
                                5., stop_radius=2.)
        assert early.end_time == early.residence == 0
        assert early.stop_cause == 'validity_boundary'
        np.testing.assert_array_equal(early.solution(0), [radius,0.,0.,0.,0.,0.])


def test_zero_velocity_inside_and_boundary_launch_inward():
    inside = integrate_orbit(lambda x: np.zeros(3), [.2,0.,0.], [0.,0.,0.],
                             3., stop_radius=2.)
    assert inside.residence == 3
    inward = integrate_orbit(lambda x: np.zeros(3), [1.,0.,0.], [-1.,0.,0.],
                             1., stop_radius=2.)
    assert inward.residence == pytest.approx(1.)
    assert not inward.crossings


def test_host_scaling_sign_and_regular_center():
    class HarmonicHost:
        rt, G, M = 10., 1., 1000.
        @staticmethod
        def psi_spline(r, derivative=0):
            return -4*r if derivative else 200-2*r*r
    host = HarmonicHost()
    scales = host_scales(host, 2.)
    assert scales['velocity_kms'] == 4.
    force = host_acceleration(host, 2., 4.)
    np.testing.assert_allclose(force(np.array([.5,.25,0.])), [-.5,-.25,0.])
    np.testing.assert_array_equal(force(np.zeros(3)), np.zeros(3))


@pytest.mark.parametrize('bad', [0., -1., np.nan])
def test_invalid_horizon_rejected(bad):
    with pytest.raises(ValueError):
        integrate_orbit(lambda x: -x, [1.,0.,0.], [0.,0.,0.], bad, stop_radius=2.)
