"""Compare radial quadrature with exact Cartesian harmonic motion."""
import numpy as np
import pytest

from analysis.residence.radial_reference import radial_residence


def harmonic_occupation(x, v, horizon, radius):
    """Exact crossing times from |x cos t + v sin t|^2, with omega=1."""
    a = (x@x+v@v)/2
    b, c = (x@x-v@v)/2, x@v
    amplitude = np.hypot(b, c)
    threshold = (radius**2-a)/amplitude
    if threshold >= 1:
        return horizon
    if threshold <= -1:
        return 0.
    phase, alpha = np.arctan2(c, b), np.arccos(threshold)
    crossings = [0., horizon]
    for k in range(-3, int(horizon/np.pi)+4):
        for sign in [-1, 1]:
            t = (phase+sign*alpha+2*np.pi*k)/2
            if 0 < t < horizon:
                crossings.append(t)
    crossings.sort()
    total = 0.
    for left, right in zip(crossings[:-1], crossings[1:]):
        t = (left+right)/2
        if np.linalg.norm(x*np.cos(t)+v*np.sin(t)) < radius:
            total += right-left
    return total


@pytest.mark.parametrize('velocity', [[np.sqrt(3.75),0,0],
                                      [-np.sqrt(3.75),0,0],
                                      [0,np.sqrt(3.75),0], [.4,1.7,.2]])
@pytest.mark.parametrize('horizon', [.07, .8, 2.3, 10., 31.2])
def test_finite_phase_matches_exact_harmonic_motion(velocity, horizon):
    x, v = np.array([.5,0,0]), np.array(velocity)
    result = radial_residence(lambda r: 20-r*r/2, x, v, horizon, order=96)
    assert result['radial_period'] == pytest.approx(np.pi, abs=2e-8)
    assert result['residence'] == pytest.approx(harmonic_occupation(x,v,horizon,1.), abs=3e-8)


def test_turning_point_and_full_empty_regions():
    x, v = np.array([2.,0,0]), np.zeros(3)
    result = radial_residence(lambda r: 20-r*r/2, x,v,7.)
    assert result['initial_phase'] == pytest.approx(np.pi/2, abs=1e-8)
    assert result['residence'] == pytest.approx(harmonic_occupation(x,v,7.,1.), abs=2e-8)
    assert radial_residence(lambda r: 20-r*r/2,x,v,7.,region_radius=3.)['residence'] == pytest.approx(7.)
    assert radial_residence(lambda r: 20-r*r/2,x,v,7.,region_radius=0.)['residence'] == 0.


def test_radial_center_launch_and_invalid_orbits():
    result = radial_residence(lambda r: 20-r*r/2,[0,0,0],[2,0,0],1.)
    assert result['initial_phase'] == 0.
    assert result['residence'] == pytest.approx(np.pi/6, abs=2e-8)
    with pytest.raises(ValueError, match='outer turning'):
        radial_residence(lambda r: 20-r*r/2,[.5,0,0],[4,0,0],1.,outer_radius=2.)
    with pytest.raises(ValueError, match='circular'):
        radial_residence(lambda r: 20-r*r/2,[1,0,0],[0,1,0],1.)
