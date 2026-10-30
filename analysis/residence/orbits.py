"""Collisionless reference trajectories and strict-region occupation.

Coordinates, velocities and times are dimensionless and mutually consistent.
The integrator supplies a deterministic zero-collision control, not dynamical
friction, a stochastic residence prediction, or an irreversible escape rule.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp


PC_PER_KMS_MYR = 0.9777922216807892


@dataclass
class OrbitResult:
    solution: Callable
    end_time: float
    stop_cause: str
    residence: float
    crossings: list
    turning_times: list
    region_touches: list
    time_grid: np.ndarray
    states: np.ndarray
    function_evaluations: int


def host_scales(host, radius):
    """Return r, v and t scales from the same adopted spherical potential."""
    if not np.isfinite(radius) or not 0 < radius < host.rt:
        raise ValueError('reference radius must lie strictly inside the host')
    v2 = -radius * float(host.psi_spline(radius, 1))
    if not np.isfinite(v2) or v2 <= 0:
        raise ValueError('reference radius requires an attractive host force')
    velocity = np.sqrt(v2)
    return {'radius_pc': float(radius), 'velocity_kms': float(velocity),
            'time_Myr': float(PC_PER_KMS_MYR * radius / velocity)}


def host_acceleration(host, radius_scale, velocity_scale):
    """Dimensionless smooth host gravity, with the center evaluated regularly.

    The relative potential is psi = constant - Phi, so psi' is negative.
    All physical species already contribute to this potential. No extra BH
    force or remnant self-force is added. Outside the DF support the spherical
    mass continues as a point mass, for solver stages near a support stop.
    """
    if (not np.isfinite([radius_scale, velocity_scale]).all()
            or radius_scale <= 0 or velocity_scale <= 0):
        raise ValueError('positive finite radius and velocity scales required')

    def acceleration(position):
        position = np.asarray(position, dtype=float)
        radius = float(np.linalg.norm(position))
        if radius == 0:
            return np.zeros(3)
        physical_radius = radius_scale * radius
        derivative = (float(host.psi_spline(physical_radius, 1))
                      if physical_radius < host.rt else
                      -host.G * host.M / physical_radius**2)
        return radius_scale / velocity_scale**2 * derivative * position / radius

    return acceleration


def _unique_times(values, horizon):
    tolerance = 1e-10 * max(1., horizon)
    output = []
    for value in sorted(float(t) for t in values):
        if not output or value - output[-1] > tolerance:
            output.append(value)
    return output


def integrate_orbit(acceleration, x0, v0, horizon, *, region_radius=1.,
                    stop_radius, rtol=1e-10, atol=1e-12, max_step=None,
                    boundary_tolerance=1e-9):
    """Integrate dx/dt=v, dv/dt=a(x), and accumulate time with |x|<R.

    Region crossings are roots located by the ODE solver; radial turning
    roots additionally identify contacts without a membership change. The
    occupation is the sum of complete intervals between those roots. An
    initial point at or outside the validity boundary stops at t=0. Region
    boundary points are outside; an orbit staying on that boundary contributes
    zero residence. ``boundary_tolerance`` resolves numerically stationary
    contacts only; it is recorded by the calling protocol and must be refined.

    As for all event solvers, max_step must resolve the shortest orbit segment;
    separate step/tolerance refinement is required to check missed crossings.
    """
    x0, v0 = np.asarray(x0, dtype=float), np.asarray(v0, dtype=float)
    if (x0.shape != (3,) or v0.shape != (3,)
            or not np.isfinite(np.r_[x0, v0]).all()):
        raise ValueError('initial position and velocity must be finite three-vectors')
    parameters = [horizon, region_radius, stop_radius, rtol, atol, boundary_tolerance]
    if not np.isfinite(parameters).all() or min(parameters) <= 0:
        raise ValueError('horizon, radii and tolerances must be positive and finite')
    if max_step is None:
        max_step = horizon / 256
    if not np.isfinite(max_step) or max_step <= 0:
        raise ValueError('max_step must be positive and finite')
    state0 = np.r_[x0, v0]
    if np.linalg.norm(x0) >= stop_radius:
        def constant(times):
            times = np.asarray(times)
            return np.broadcast_to(state0.reshape((6,) + (1,) * times.ndim),
                                   (6,) + times.shape).copy()
        return OrbitResult(constant, 0., 'validity_boundary', 0., [], [], [],
                           np.array([0.]), state0[:, None], 0)

    def rhs(t, state):
        force = np.asarray(acceleration(state[:3]), dtype=float)
        if force.shape != (3,) or not np.isfinite(force).all():
            raise ValueError('acceleration must return a finite three-vector')
        return np.r_[state[3:], force]

    def region(t, state):
        return float(state[:3] @ state[:3] - region_radius**2)

    def stop(t, state):
        return float(state[:3] @ state[:3] - stop_radius**2)

    stop.terminal, stop.direction = True, 1

    def turning(t, state):
        return float(state[:3] @ state[3:])

    result = solve_ivp(rhs, (0., horizon), state0, method='DOP853',
                       dense_output=True, events=(region, stop, turning),
                       rtol=rtol, atol=atol, max_step=max_step)
    if not result.success:
        raise RuntimeError(f'orbit integration failed: {result.message}')
    end = float(result.t[-1])
    time_tolerance = 1e-10 * max(1., horizon)
    turns = _unique_times([t for t in result.t_events[2]
                          if time_tolerance < t < end-time_tolerance], horizon)
    # Zero-speed or circular paths have x.v=0 throughout. They have no isolated
    # radial turning event; retain only roots with a resolved sign change.
    isolated_turns = []
    for i, time in enumerate(turns):
        left = (turns[i-1] if i else 0.)
        right = (turns[i+1] if i+1 < len(turns) else end)
        before = result.sol((left+time)/2)
        after = result.sol((time+right)/2)
        scale = max(1., np.linalg.norm(before[:3])*np.linalg.norm(before[3:]),
                    np.linalg.norm(after[:3])*np.linalg.norm(after[3:]))
        qb, qa = before[:3] @ before[3:], after[:3] @ after[3:]
        if qb*qa < 0 and min(abs(qb), abs(qa)) > boundary_tolerance*scale:
            isolated_turns.append(time)
    roots = _unique_times([t for t in result.t_events[0]
                          if time_tolerance < t < end-time_tolerance], horizon)
    touches = [t for t in isolated_turns
               if abs(np.linalg.norm(result.sol(t)[:3])-region_radius)
               <= boundary_tolerance*max(1., region_radius)]
    boundaries = _unique_times([0., *roots, *touches, end], horizon)
    sides = []
    residence = 0.
    radius_tolerance = boundary_tolerance * max(1., region_radius)
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        radius = np.linalg.norm(result.sol((left+right)/2)[:3])
        inside = bool(radius < region_radius-radius_tolerance)
        sides.append(inside)
        if inside:
            residence += right-left
    crossings = []
    for i, time in enumerate(boundaries[1:-1]):
        if sides[i] != sides[i+1]:
            crossings.append({'time': float(time),
                              'direction': 'entry' if sides[i+1] else 'exit'})
    # A touch is not an entry or exit. Do not count an initial boundary launch
    # as a crossing; its direction is supplied by subsequent interior motion.
    touches = [t for t in touches if not any(abs(t-c['time']) <= time_tolerance
                                            for c in crossings)]
    return OrbitResult(result.sol, end,
                       'validity_boundary' if len(result.t_events[1]) else 'horizon',
                       float(residence), crossings, isolated_turns, touches,
                       result.t, result.y, result.nfev)
