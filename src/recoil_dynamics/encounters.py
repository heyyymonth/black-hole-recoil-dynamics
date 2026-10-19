"""Velocity-weighted close-passage integrals for an imposed Plummer tracer bath.

The test orbit is never perturbed. Velocities use sqrt(G M / b), time uses
sqrt(b^3 / (G M)). These are encounter optical-depth ingredients, not outcomes.
"""
from functools import cached_property
import math

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import beta, betainc

from .plummer import PlummerOrbit


def _finite(value, name, positive=False):
    value = float(value)
    if not math.isfinite(value) or (value <= 0 if positive else value < 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return value


_BETA_NORM = beta(1.5, 4.5)


def _speed_moment(power, fraction, upper=False):
    """Integral of y**power against the normalized Plummer speed density."""
    a = (power + 3) / 2
    # Complement identity is stable near unity and works with pinned SciPy 1.10.
    incomplete = (betainc(4.5, a, (1-fraction)*(1+fraction)) if upper
                  else betainc(a, 4.5, fraction*fraction))
    return beta(a, 4.5) / _BETA_NORM * incomplete


def relative_speed_moments(speed, escape_speed):
    """Return <w>, <1/w> for isotropic Plummer target velocities.

    For y=u/escape_speed, p(y)=512/(7*pi)*y²*(1-y²)^(7/2).
    Angular integration gives <w>_angle=max(v,u)+min(v,u)²/(3*max(v,u))
    and <1/w>_angle=1/max(v,u). Remaining speed integrals use incomplete beta
    functions split at u=v. The singular single point u=v=0 is integrable.
    """
    v = _finite(speed, "speed")
    ve = _finite(escape_speed, "escape speed", positive=True)
    s = v / ve
    if s == 0:
        return ve * 1024/(693*math.pi), 512/(63*math.pi*ve)
    if s >= 1:
        return v + ve*ve/(12*v), 1/v
    below0 = _speed_moment(0, s)
    below2 = _speed_moment(2, s)
    above1 = _speed_moment(1, s, upper=True)
    above_inverse = _speed_moment(-1, s, upper=True)
    mean = s*below0 + below2/(3*s) + above1 + s*s*above_inverse/3
    inverse = below0/s + above_inverse
    return ve*mean, inverse/ve


class EncounterOrbit:
    """Density-weighted relative-speed moments along a ballistic orbit.

    integrals(H) returns J_g = integral rho/rho0 <w> dt and
    J_f = integral rho/rho0 <1/w> dt. No physical host or tidal cutoff is
    imposed here; the experiment runner declares those separately.
    """

    def __init__(self, eta, epsabs=1e-10, epsrel=1e-10):
        self.orbit = PlummerOrbit(eta, epsabs=epsabs, epsrel=epsrel)
        self._epsabs = float(epsabs)
        self._epsrel = float(epsrel)

    def _state(self, theta):
        xa = self.orbit.apocenter
        x = xa*math.sin(theta)
        s = math.hypot(1, x)
        sa = 1/(1-self.orbit.eta**2)
        # Rationalize 1/s - 1/sa: no cancellation at the turning point.
        v = math.sqrt(2/(s*sa*(s+sa))) * xa * abs(math.cos(theta))
        dt_dtheta = math.sqrt(s*sa*(s+sa)/2)
        return s, v, dt_dtheta

    def _outgoing(self, theta, component):
        if theta <= 0:
            return 0.0

        def integrand(angle):
            s, v, dt = self._state(angle)
            if component is None:
                return dt
            moments = relative_speed_moments(v, math.sqrt(2/s))
            return dt * s**-5 * moments[component]

        edges = [0.0]
        for radius in [.01, .1, 1., 10., 100.]:
            if radius < self.orbit.apocenter:
                angle = math.asin(radius/self.orbit.apocenter)
                if angle < theta:
                    edges.append(angle)
        edges.append(theta)
        return math.fsum(quad(integrand, a, b,
                             epsabs=self._epsabs/(len(edges)-1),
                             epsrel=self._epsrel, limit=150)[0]
                         for a, b in zip(edges[:-1], edges[1:]))

    @cached_property
    def _cycle(self):
        return tuple(2*self._outgoing(math.pi/2, k) for k in (0, 1))

    def integrals(self, horizon):
        h = _finite(horizon, "horizon")
        if h == 0:
            return 0., 0.
        if self.orbit.eta == 0:
            return tuple(h*m for m in relative_speed_moments(0, math.sqrt(2)))
        period = self.orbit.period
        n = math.floor(h/period)
        remainder = h-n*period
        outgoing_time = min(remainder, period-remainder)
        if outgoing_time == 0:
            theta = 0.0
        elif outgoing_time >= period/2:
            theta = math.pi/2
        else:
            theta = brentq(lambda angle: self._outgoing(angle, None)-outgoing_time,
                           0., math.pi/2, xtol=5e-15, rtol=1e-14)
        partial = tuple(self._outgoing(theta, k) for k in (0, 1))
        if remainder > period/2:
            partial = tuple(full-part for full, part in zip(self._cycle, partial))
        return tuple(n*full+part for full, part in zip(self._cycle, partial))


def physical_scales(config):
    """Convert the declared illustrative host into pc, km/s, Myr scales."""
    names = ("host_mass_msun", "plummer_scale_pc", "remnant_mass_msun",
             "target_binary_mass_msun", "target_binary_count",
             "target_binary_semimajor_axis_au", "passage_radius_au",
             "horizon_over_tdyn", "apocenter_limit_over_half_mass_radius")
    values = {name: _finite(config[name], name, positive=True) for name in names}
    constants = {name: _finite(value, name, positive=True)
                 for name, value in config['illustrative_constants'].items()}
    mass = values['host_mass_msun']
    b = values['plummer_scale_pc']
    radius = values['passage_radius_au']*constants['au_in_km']/constants['pc_in_km']
    velocity = math.sqrt(constants['G_pc_kms2_per_msun']*mass/b)
    tdyn = b*constants['pc_in_km']/velocity/constants['seconds_per_myr']
    n0 = 3*values['target_binary_count']/(4*math.pi*b**3)
    half_mass_radius = 1/math.sqrt(2**(2/3)-1)
    focusing = 2*(values['remnant_mass_msun']+values['target_binary_mass_msun'])/mass*b/radius
    return {
        'velocity_unit_kms': velocity,
        'central_escape_kms': math.sqrt(2)*velocity,
        'dynamical_time_myr': tdyn,
        'horizon_myr': values['horizon_over_tdyn']*tdyn,
        'half_mass_radius_over_b': half_mass_radius,
        'apocenter_limit_over_b': values['apocenter_limit_over_half_mass_radius']*half_mass_radius,
        'central_binary_density_pc_minus3': n0,
        'binary_mass_fraction': values['target_binary_count']*values['target_binary_mass_msun']/mass,
        'passage_radius_over_binary_semimajor_axis': values['passage_radius_au']/values['target_binary_semimajor_axis_au'],
        'focusing_coefficient': focusing,
        # V*t_d=b, so this factor converts BOTH dimensionless integrals once.
        'optical_depth_prefactor': n0*math.pi*radius**2*b,
    }


def experiment_rows(config):
    """Four declared cases; reject rather than silently use inadmissible kicks."""
    scales = physical_scales(config)
    h = float(config['horizon_over_tdyn'])
    coefficient = scales['focusing_coefficient']
    prefactor = scales['optical_depth_prefactor']
    g0, f0 = relative_speed_moments(0, math.sqrt(2))
    tau0 = prefactor*h*(g0+coefficient*f0)
    tolerance = _finite(config['proxy_relative_error_tolerance'], "proxy tolerance", positive=True)
    rows = []
    for eta in config['eta']:
        enc = EncounterOrbit(eta, epsabs=config['epsabs'], epsrel=config['epsrel'])
        orb = enc.orbit
        if orb.apocenter > scales['apocenter_limit_over_b']:
            raise ValueError(f"eta={eta}: apocenter exceeds declared admission radius")
        jg, jf = enc.integrals(h)
        geom = prefactor*jg
        focus = prefactor*coefficient*jf
        tau = geom+focus
        density = orb.exposure(h)/h
        ratio = tau/tau0
        error = density/ratio-1
        rows.append({
            'eta': eta,
            'kick_kms': float(eta)*scales['central_escape_kms'],
            'apocenter_over_b': orb.apocenter,
            'horizon_over_tdyn': h,
            'horizon_myr': scales['horizon_myr'],
            'density_ratio': density,
            'geometric_integral': jg,
            'focusing_integral': jf,
            'geometric_ratio': jg/(h*g0),
            'focusing_ratio': jf/(h*f0),
            'tau_geometric': geom,
            'tau_focusing': focus,
            'tau_total': tau,
            'tau_stationary': tau0,
            'encounter_ratio': ratio,
            'density_proxy_tau': density*tau0,
            'proxy_relative_error': error,
            'focus_fraction': focus/tau,
            'exceeds_tolerance': int(abs(error) > tolerance),
        })
    return scales, rows
