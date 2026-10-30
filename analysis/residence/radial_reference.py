"""Independent finite-horizon occupation quadrature for a spherical orbit.

The potential convention is Phi=-psi, and lengths, velocities and times must
use one consistent dimensionless (or physical natural-unit) system. Motion is
conservative in a time-independent spherical potential, with no collisions,
external forces or stopping boundaries. A single noncircular bound radial
interval is required. A radial trajectory passes smoothly through a finite
central potential; its radius then begins the next outward radial phase.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq, minimize_scalar


def radial_residence(psi, x0, v0, horizon, region_radius=1.0,
                     outer_radius=4.0, order=96):
    """Return turning radii, radial period and actual finite-time residence.

    ``outer_radius`` must exceed the initial radius and outer turning point,
    with negative squared radial speed there. ``psi`` accepts a scalar radius.
    Initial radial direction and phase are retained, including launches at a
    turning point. This is not a cycle-average approximation. Circular orbits
    are deliberately rejected rather than handled as ill-conditioned double
    roots. Occupation is in r < region_radius; isolated equality times have
    zero measure for the noncircular orbits considered here.
    """
    x0, v0 = np.asarray(x0, dtype=float), np.asarray(v0, dtype=float)
    if (x0.shape != (3,) or v0.shape != (3,) or
            not np.isfinite(x0).all() or not np.isfinite(v0).all()):
        raise ValueError('initial position and velocity must be finite three-vectors')
    if (not np.isfinite([horizon, region_radius, outer_radius]).all() or
            horizon < 0 or region_radius < 0 or outer_radius <= 0):
        raise ValueError('finite nonnegative horizon/region and positive outer bracket required')
    if not isinstance(order, (int, np.integer)) or order < 16:
        raise ValueError('quadrature order must be an integer >= 16')
    r0 = float(np.linalg.norm(x0))
    if outer_radius <= r0:
        raise ValueError('outer_radius must exceed the initial radius')

    def potential(radius):
        value = float(psi(float(radius)))
        if not np.isfinite(value):
            raise ValueError('psi must be finite on the sampled radial interval')
        return value

    energy = float(np.dot(v0, v0)/2-potential(r0))
    angular_momentum = float(np.linalg.norm(np.cross(x0, v0)))
    ell2 = angular_momentum**2
    radial_velocity = float(np.dot(x0, v0)/r0) if r0 > 0 else float(np.linalg.norm(v0))

    def speed_squared(radius):
        if radius == 0 and ell2 > 0:
            return -np.inf
        centrifugal = ell2/radius**2 if ell2 else 0.
        return 2*(energy+potential(radius))-centrifugal

    if speed_squared(outer_radius) >= 0:
        raise ValueError('outer_radius does not bracket the bound outer turning point')
    root_xtol = max(np.finfo(float).tiny, 2e-14*outer_radius)
    if ell2 == 0:
        pericenter = 0.
        if speed_squared(0.) <= 0:
            raise ValueError('stationary or unsupported central orbit')
        apocenter = brentq(speed_squared, 0., outer_radius,
                          xtol=root_xtol, rtol=8*np.finfo(float).eps)
    else:
        lower = max(r0*1e-12, np.finfo(float).tiny**.25)
        if speed_squared(lower) >= 0:
            raise ValueError('unable to bracket the inner turning point')
        peak = minimize_scalar(lambda r: -speed_squared(r),
                               bounds=(lower, outer_radius), method='bounded',
                               options={'xatol': root_xtol})
        scale = max(float(np.dot(v0, v0)), abs(energy+potential(r0)), 1.)
        if not peak.success or speed_squared(peak.x) <= 1e-12*scale:
            raise ValueError('circular or unresolved noncircular radial interval')
        pericenter = brentq(speed_squared, lower, peak.x,
                           xtol=root_xtol, rtol=8*np.finfo(float).eps)
        apocenter = brentq(speed_squared, peak.x, outer_radius,
                          xtol=root_xtol, rtol=8*np.finfo(float).eps)
    if not pericenter-10*root_xtol <= r0 <= apocenter+10*root_xtol:
        raise ValueError('initial state is outside the solved radial interval')

    delta = apocenter-pericenter
    gauss, weights = leggauss(order)
    psi_apo = potential(apocenter)
    psi_peri = potential(pericenter) if ell2 else None

    def stable_speed_squared(radius):
        # Subtract the nearest turning-point energy. This makes the root
        # exactly zero in the integrand without subtracting two O(E) values.
        if ell2 and radius-pericenter < apocenter-radius:
            reference, psi_reference = pericenter, psi_peri
        else:
            reference, psi_reference = apocenter, psi_apo
        value = 2*(potential(radius)-psi_reference)
        if ell2:
            value += ell2*(radius-reference)*(radius+reference)/(reference*radius)**2
        return value

    def time_to(radius):
        if radius <= pericenter:
            return 0.
        theta_end = np.arcsin(np.sqrt(np.clip((radius-pericenter)/delta, 0., 1.)))
        theta = (gauss+1)*theta_end/2
        radial = pericenter+delta*np.sin(theta)**2
        jacobian = 2*delta*np.sin(theta)*np.cos(theta)
        radial_speed2 = np.array([stable_speed_squared(r) for r in radial])
        if np.any(radial_speed2 <= 0) or not np.isfinite(radial_speed2).all():
            raise ValueError('nonpositive radial speed within quadrature interval')
        return float(theta_end/2*np.sum(weights*jacobian/np.sqrt(radial_speed2)))

    half_period = time_to(apocenter)
    radial_period = 2*half_period
    radius_tolerance = 20*root_xtol
    if abs(r0-pericenter) <= radius_tolerance:
        initial_phase = 0.
    elif abs(r0-apocenter) <= radius_tolerance:
        initial_phase = half_period
    else:
        initial_phase = time_to(r0)
        if radial_velocity < 0:
            initial_phase = radial_period-initial_phase
    if region_radius <= pericenter:
        region_outward_time = 0.
    elif region_radius >= apocenter:
        region_outward_time = half_period
    else:
        region_outward_time = time_to(region_radius)

    def accumulated_occupation(phase):
        cycles = np.floor(phase/radial_period)
        remainder = phase-cycles*radial_period
        return (cycles*2*region_outward_time
                + min(remainder, region_outward_time)
                + max(remainder-(radial_period-region_outward_time), 0.))

    residence = float(accumulated_occupation(initial_phase+horizon)
                      - accumulated_occupation(initial_phase))
    # Bound only sub-ulp accumulation roundoff, never integrate an invalid path.
    if residence < -1e-12*max(horizon, 1.) or residence > horizon+1e-12*max(horizon, 1.):
        raise ArithmeticError('occupation lies outside the finite horizon')
    residence = float(np.clip(residence, 0., horizon))
    return {'pericenter': float(pericenter), 'apocenter': float(apocenter),
            'radial_period': float(radial_period), 'half_period': float(half_period),
            'initial_phase': float(initial_phase),
            'region_outward_time': float(region_outward_time),
            'residence': residence,
            'residence_fraction': residence/horizon if horizon else 0.,
            'energy': energy, 'angular_momentum': angular_momentum,
            'quadrature_order': int(order),
            'scope': 'conservative noncircular bound spherical orbit; no stopping boundaries'}
