"""Radial orbits in a static Plummer potential, in units b and sqrt(b^3/GM).

Production calculations use energy quadrature, not time-stepped trajectories.
The period is center-to-center: density has period P, signed position 2P.
No dynamical friction, encounters, tides or host evolution is included.
"""
from functools import cached_property
import math
from scipy.integrate import quad
from scipy.optimize import brentq


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def core_exit_threshold(core_radius):
    """Largest central kick fraction whose apocenter is inside radius c b."""
    c = _positive(core_radius, "core radius")
    s = math.hypot(1, c)
    # Rationalized form avoids subtracting nearly equal numbers at small c.
    return (c / s) * math.sqrt(s / (s + 1))


class PlummerOrbit:
    """A bound radial kick with 0 <= eta < 1; eta=0 is stationary.

    Energy-quadrature accuracy is tested through eta=0.9999. Radii and horizons
    are dimensionless. The special stationary object has no return period.
    Parameters are read-only to keep cached orbit integrals consistent.
    """

    def __init__(self, eta, epsabs=1e-11, epsrel=1e-11):
        eta = float(eta)
        if not math.isfinite(eta) or not 0 <= eta < 1:
            raise ValueError("eta must be finite and in [0, 1)")
        self._eta = eta
        self._epsabs = _positive(epsabs, "absolute tolerance")
        self._epsrel = _positive(epsrel, "relative tolerance")

    @property
    def eta(self):
        return self._eta

    @property
    def apocenter(self):
        e = self.eta
        return e * math.sqrt(2 - e*e) / (1 - e*e)

    def _dt_dtheta(self, theta):
        sa = 1 / (1 - self.eta*self.eta)
        s = math.hypot(1, self.apocenter * math.sin(theta))
        return math.sqrt(s * sa * (s + sa) / 2)

    def _outgoing_integral(self, theta, weighted=False):
        if theta <= 0:
            return 0.0
        xa = self.apocenter
        def integrand(angle):
            dt = self._dt_dtheta(angle)
            if weighted:
                return dt / math.hypot(1, xa * math.sin(angle))**5
            return dt
        # Resolve the narrow high-density part of near-escape excursions.
        split = [0.0]
        for radius in [.01, .1, 1., 10., 100.]:
            if radius < xa:
                angle = math.asin(radius / xa)
                if angle < theta:
                    split.append(angle)
        split.append(theta)
        return math.fsum(quad(integrand, a, b, epsabs=self._epsabs / (len(split)-1),
                             epsrel=self._epsrel, limit=150)[0]
                         for a, b in zip(split[:-1], split[1:]))

    @cached_property
    def period(self):
        if self.eta == 0:
            return None
        return 2 * self._outgoing_integral(math.pi / 2)

    @cached_property
    def _cycle_exposure(self):
        if self.eta == 0:
            raise ValueError("A stationary object has no orbital cycle")
        return 2 * self._outgoing_integral(math.pi / 2, weighted=True)

    def cycle_exposure(self):
        return self._cycle_exposure

    def _core_exit_time(self, core_radius):
        c = _positive(core_radius, "core radius")
        if self.eta == 0:
            raise ValueError("A stationary object has no orbital cycle")
        if c >= self.apocenter:
            return self.period / 2
        return self._outgoing_integral(math.asin(c / self.apocenter))

    def cycle_core_time(self, core_radius):
        return 2 * self._core_exit_time(core_radius)

    @staticmethod
    def _horizon(horizon):
        h = float(horizon)
        if not math.isfinite(h) or h < 0:
            raise ValueError("horizon must be finite and nonnegative")
        return h

    def _partial_outgoing_exposure(self, time):
        if time <= 0:
            return 0.0
        if time >= self.period / 2:
            return self.cycle_exposure() / 2
        theta = brentq(lambda a: self._outgoing_integral(a) - time,
                       0., math.pi / 2, xtol=5e-15, rtol=1e-14)
        return self._outgoing_integral(theta, weighted=True)

    def exposure(self, horizon):
        """Integral of normalized Plummer density along the orbit, to H."""
        h = self._horizon(horizon)
        if h == 0 or self.eta == 0:
            return h
        n = math.floor(h / self.period)
        remainder = h - n*self.period
        if remainder <= self.period / 2:
            partial = self._partial_outgoing_exposure(remainder)
        else:
            partial = self.cycle_exposure() - self._partial_outgoing_exposure(self.period-remainder)
        return n*self.cycle_exposure() + partial

    def core_time(self, horizon, core_radius):
        """Exact integrated time inside the operational sphere |x| < c."""
        h = self._horizon(horizon)
        c = _positive(core_radius, "core radius")
        if h == 0 or self.eta == 0 or c >= self.apocenter:
            return h
        tc = self._core_exit_time(c)
        n = math.floor(h / self.period)
        remainder = h - n*self.period
        partial = min(remainder, tc) + max(0., remainder-(self.period-tc))
        return n*2*tc + partial
