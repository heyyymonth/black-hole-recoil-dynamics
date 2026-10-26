"""Static boundary diagnostics. Units: pc, km/s, Msun; no trajectory evolution."""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

class EnclosedHost:
    def __init__(self, host, stride=1):
        self.host = host
        indices = np.unique(np.r_[np.arange(0, len(host.r), stride), len(host.r)-1])
        # Interpolate M/r^3, whose finite central limit follows from the density.
        # Interpolating M directly distorts the cubic core scaling on a sparse mesh.
        mass = host.mcj[host.physical].sum(axis=0)
        mean = np.empty_like(host.r)
        mean[0] = 4*np.pi/3*host.rhoj[host.physical, 0].sum()
        mean[1:] = mass[1:]/host.r[1:]**3
        self.mean = PchipInterpolator(host.r[indices], mean[indices])
        self.G, self.M, self.rt = host.G, host.M, host.rt

    def enclosed(self, r):
        r = np.asarray(r)
        if not np.isfinite(r).all() or np.any(r < 0):
            raise ValueError('Radius must be finite and nonnegative')
        return np.where(r >= self.rt, self.M, r**3*self.mean(np.minimum(r, self.rt)))

    def potential(self, r):
        r = np.asarray(r)
        inside = -self.host.psi(r) - self.G*self.M/self.rt
        return np.where(r <= self.rt, inside, -self.G*self.M/np.maximum(r, self.rt))


def differential_acceleration(Rvec, rvec, vc):
    """Exact Galactic difference in translating, NONROTATING axes."""
    Rvec, rvec = np.asarray(Rvec), np.asarray(rvec)
    other = Rvec + rvec
    return vc**2*(Rvec/np.sum(Rvec**2) - other/np.sum(other**2, axis=-1, keepdims=True))


def tidal_tensor(Rvec, vc):
    Rvec = np.asarray(Rvec)
    R = np.linalg.norm(Rvec)
    n = Rvec/R
    return vc**2/R**2*(2*np.outer(n, n)-np.eye(3))


def max_force_ratio(r, enclosed, G, R, vc):
    """Exact all-direction bound for r<R, normalized to spherical cluster force."""
    r = np.asarray(r)
    if np.any((r <= 0) | (r >= R)):
        raise ValueError('Require 0 < r < Galactic distance')
    return vc**2*r**3/(G*enclosed(r)*R*(R-r))


def root(function, lower, upper):
    if function(lower)*function(upper) >= 0:
        raise ValueError('No bracketed outer boundary within the declared domain')
    return brentq(function, lower, upper, xtol=1e-10, rtol=1e-13)


def jacobi_radius(enclosed, G, lam, lower, upper):
    if lam <= 0:
        raise ValueError('A positive extensive effective eigenvalue is required')
    return root(lambda r: G*enclosed(r)/r**3-lam, lower, upper)


def diagnostics(host, R, vc, tolerances):
    lower, upper = 1e-3*host.host.rh, min(host.rt, .2*R)
    tau = vc**2/R**2
    rj = jacobi_radius(host.enclosed, host.G, 2*tau, lower, upper)
    result = {
        'R_pc': R, 'vc_kms': vc, 'tidal_norm_kms2_pc2': tau,
        'effective_extensive_eigenvalue_kms2_pc2': 2*tau,
        'rJ_point_mass_pc': (host.G*host.M/(2*tau))**(1/3),
        'rJ_extended_pc': rj, 'rJ_over_rh': rj/host.host.rh,
        'mass_outside_rJ_fraction': 1-float(host.enclosed(rj))/host.M,
        'root_relative_residual': abs(host.G*host.enclosed(rj)/rj**3/(2*tau)-1),
        'eta_at_rh': float(max_force_ratio(host.host.rh, host.enclosed, host.G, R, vc)),
        'eta_at_3rh': float(max_force_ratio(3*host.host.rh, host.enclosed, host.G, R, vc)),
        'hill_central_opening_speed_kms': np.sqrt(2*(host.potential(rj)-tau*rj*rj-host.potential(0))),
    }
    for side, sign in [('inward', -1), ('outward', 1)]:
        # Exact circular rotating-frame radial balance, including centrifugal force.
        balance = lambda r: host.G*host.enclosed(r)/r**3 - vc**2/(R*(R+sign*r))-tau
        rs = root(balance, lower, upper)
        result[f'exact_{side}_saddle_pc'] = rs
        result[f'exact_{side}_relative_residual'] = abs(balance(rs)/(2*tau))
    for eps in tolerances:
        rv = root(lambda r: max_force_ratio(r, host.enclosed, host.G, R, vc)-eps, lower, upper)
        label = f'{eps:g}'
        result[f'r_valid_{label}_pc'] = rv
        result[f'valid_{label}_relative_residual'] = abs(max_force_ratio(rv, host.enclosed, host.G, R, vc)/eps-1)
    return {k: float(v) for k, v in result.items()}
