"""Exact frozen spherical-shell monopoles and collisionless radial returns.

G=1. Radii and masses may use any mutually consistent units. If radii are
in b and masses in M, times are in sqrt(b**3/(G*M)). Each supplied object
is spread over an infinitesimally thin spherical shell; this is not an
interpolation to a smooth mass profile. Angular coordinates and velocities
do not enter this frozen monopole approximation.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad


@lru_cache(maxsize=7)
def _unit_gauss(order: int):
    nodes, weights = leggauss(order)
    return (nodes + 1) / 2, weights / 2


class ShellPotential:
    """Frozen monopole of positive-mass shells at positive radii.

    Unsorted inputs are sorted and equal-radius shells are combined. Potential
    is continuous. ``acceleration`` at a shell uses the inward-side limit;
    its jump is physical for the stated shell model. A shell boundary has zero
    measure in a return-time integral. The innermost region is force-free.
    """

    def __init__(self, radii, masses):
        r = np.asarray(radii, dtype=float)
        m = np.asarray(masses, dtype=float)
        if (r.ndim != 1 or m.ndim != 1 or r.shape != m.shape or r.size == 0
                or not np.isfinite(r).all() or not np.isfinite(m).all()
                or np.any(r <= 0) or np.any(m <= 0)):
            raise ValueError("Shell radii and masses must be equal nonempty 1-D arrays of finite positive values")
        order = np.argsort(r, kind="stable")
        r, m = r[order], m[order]
        radii_unique, first = np.unique(r, return_index=True)
        masses_combined = np.add.reduceat(m, first)
        self.radii = radii_unique
        self.masses = masses_combined
        self._enclosed = np.concatenate(([0.0], np.cumsum(self.masses)))
        self._outer = np.concatenate((np.cumsum((self.masses / self.radii)[::-1])[::-1], [0.0]))
        self.total_mass = float(self._enclosed[-1])
        self.central_potential = -float(self._outer[0])
        if not math.isfinite(self.total_mass) or not math.isfinite(self.central_potential):
            raise ValueError("Shell totals must be finite in the chosen units")
        self.escape_speed = math.sqrt(-2 * self.central_potential)
        self._edge_potential = -(self._enclosed[:-1] / self.radii + self._outer[:-1])
        for array in (self.radii, self.masses, self._enclosed, self._outer, self._edge_potential):
            array.setflags(write=False)

    @staticmethod
    def _radius(radius):
        r = np.asarray(radius, dtype=float)
        if not np.isfinite(r).all() or np.any(r < 0):
            raise ValueError("Evaluation radii must be finite and nonnegative")
        return r

    def potential(self, radius):
        """Potential with zero at infinity; accepts a scalar or an array."""
        r = self._radius(radius)
        inside = np.searchsorted(self.radii, r, side="left")
        inner_term = np.divide(self._enclosed[inside], r, out=np.zeros_like(r), where=r > 0)
        result = -(inner_term + self._outer[inside])
        return float(result) if result.ndim == 0 else result

    def acceleration(self, radius):
        """Inward radial acceleration; inward-side limit at exact shell radii."""
        r = self._radius(radius)
        inside = np.searchsorted(self.radii, r, side="left")
        result = np.divide(-self._enclosed[inside], r * r, out=np.zeros_like(r), where=r > 0)
        return float(result) if result.ndim == 0 else result

    def _kick_energy(self, kick_speed: float):
        v = float(kick_speed)
        if not math.isfinite(v) or v < 0 or v >= self.escape_speed:
            raise ValueError("Kick speed must be finite, nonnegative, and below central escape speed")
        eta = v / self.escape_speed
        # Factoring the difference avoids subtraction of nearly equal energies.
        return v, self.central_potential * (1 - eta) * (1 + eta)

    def apocenter(self, kick_speed: float) -> float:
        """Outer turning radius for a bound kick; zero for a stationary object."""
        v, energy = self._kick_energy(kick_speed)
        if v == 0:
            return 0.0
        inside = int(np.searchsorted(self._edge_potential, energy, side="right"))
        inside = max(1, inside)
        denominator = -(energy + self._outer[inside])
        result = float(self._enclosed[inside] / denominator)
        # At speeds below the resolvable energy increment the outer excursion
        # rounds to the first shell. Its dominant central free flight is still
        # well-defined and is computed from the supplied kick, not this energy.
        return max(float(self.radii[0]), result)

    def first_return(self, kick_speed: float, *, epsabs: float = 1e-10,
                     epsrel: float = 1e-10) -> float:
        """First center-to-center time, without drag or encounters.

        The central free-flight interval is exact. Outside it, the integral is
        split at every shell. A quadratic change of variable regularizes each
        upper endpoint; energy gaps are accumulated from the turning point to
        avoid subtracting nearly equal potentials. Vectorized Gaussian rules
        are checked against increasing orders. Adaptive quadrature is a
        fallback for any segments whose convergence remains insufficient.

        A zero kick has no return event and raises ValueError. Unlike a smooth
        cored potential, a finite shell model has a force-free central region:
        its small-positive-kick return time grows as 2*r_min/kick_speed.
        """
        v, _ = self._kick_energy(kick_speed)
        if v == 0:
            raise ValueError("A stationary object has no first-return event")
        if (not math.isfinite(epsabs) or not math.isfinite(epsrel)
                or epsabs <= 0 or epsrel <= 0):
            raise ValueError("Quadrature tolerances must be finite and positive")
        ra = self.apocenter(v)
        central_time = float(self.radii[0] / v)
        if ra <= self.radii[0]:
            return 2 * central_time

        edges = np.concatenate((self.radii[self.radii < ra], [ra]))
        lower, upper = edges[:-1], edges[1:]
        width = upper - lower
        enclosed = self._enclosed[np.searchsorted(self.radii, lower, side="right")]
        interval_drop = enclosed * width / (lower * upper)
        # gap_upper[i] = Phi(ra)-Phi(upper[i]); construct directly from later
        # segments rather than subtracting a large segment from a cumulative sum.
        gap_upper = np.zeros_like(width)
        if width.size > 1:
            gap_upper[:-1] = np.cumsum(interval_drop[:0:-1])[::-1]

        def segment_rules(order):
            u, weights = _unit_gauss(order)
            displacement = width[:, None] * u[None, :] ** 2
            r = upper[:, None] - displacement
            gap = gap_upper[:, None] + enclosed[:, None] * displacement / (r * upper[:, None])
            values = 2 * width[:, None] * u[None, :] / np.sqrt(2 * gap)
            return values @ weights

        previous = segment_rules(16)
        for order in (32, 64, 128, 256):
            current = segment_rules(order)
            differences = np.abs(current - previous)
            tolerance = max(epsabs / 2, epsrel * (central_time + float(np.sum(current))))
            if float(np.sum(differences)) <= tolerance:
                return 2 * (central_time + float(np.sum(current)))
            previous = current

        # Only unresolved segments require Python-level adaptive integration.
        unresolved = differences > tolerance / len(width)
        for i in np.flatnonzero(unresolved):
            def integrand(u):
                displacement = width[i] * u * u
                r = upper[i] - displacement
                if gap_upper[i] == 0:
                    return math.sqrt(2 * width[i] * r * upper[i] / enclosed[i])
                gap = gap_upper[i] + enclosed[i] * displacement / (r * upper[i])
                return 2 * width[i] * u / math.sqrt(2 * gap)

            current[i] = quad(integrand, 0, 1, epsabs=epsabs / (2 * len(width)),
                              epsrel=epsrel, limit=200)[0]
        return 2 * (central_time + float(np.sum(current)))
