"""Finite-angle elastic encounters and their matched first two moments.

The rate is F(u) |v-u| d^3u b db dphi. Impact parameters refer to the
asymptotic relative orbit; the background is local, frozen and independently
resampled. These functions do not supply encounter cutoffs, binaries,
relativistic capture, correlations, or a self-consistent finite population.

Lengths are pc, velocities km/s and masses Msun. ``impact_moments`` returns
coefficients per unit number density in natural time pc/(km/s). Bath methods
return rates and velocity moments per Myr. D is the RAW second jump moment,
the coefficient of one half the velocity Hessian in the backward generator.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss

PC_PER_KMS_MYR = 0.9777922216807892


def _vectors(relative_velocity):
    w = np.asarray(relative_velocity, dtype=float)
    if w.ndim < 1 or w.shape[-1] != 3 or not np.isfinite(w).all():
        raise ValueError("relative velocities must be finite vectors of length three")
    speed = np.linalg.norm(w, axis=-1)
    direction = np.divide(w, speed[..., None], out=np.zeros_like(w),
                          where=speed[..., None] > 0)
    return w, speed, direction


def _parameters(speed, field_mass, remnant_mass, G, b_min, b_max):
    field_mass, remnant_mass, b_min, b_max = np.broadcast_arrays(
        np.asarray(field_mass, dtype=float), np.asarray(remnant_mass, dtype=float),
        np.asarray(b_min, dtype=float), np.asarray(b_max, dtype=float), speed)[:4]
    if (not np.isfinite(G) or G <= 0 or
            not all(np.isfinite(a).all() for a in (field_mass, remnant_mass, b_min, b_max)) or
            np.any(field_mass <= 0) or np.any(remnant_mass <= 0) or
            np.any(b_min < 0) or np.any(b_max < b_min)):
        raise ValueError("positive masses/G and finite 0 <= b_min <= b_max required")
    return field_mass, remnant_mass, b_min, b_max


def scatter_impulse(relative_velocity, field_mass, remnant_mass, G,
                    impact_parameter, azimuth):
    """Remnant velocity increment for an exact Newtonian elastic deflection.

    Inputs broadcast over the vector leading dimensions. The arbitrary
    perpendicular basis only chooses the origin of a uniform impact azimuth.
    At zero relative speed the increment is exactly zero. Closest approaches
    outside the Newtonian single-object model must be excluded by the caller.
    """
    _, speed, direction = _vectors(relative_velocity)
    mass, remnant, b, _ = _parameters(speed, field_mass, remnant_mass, G,
                                     impact_parameter, impact_parameter)
    speed, mass, remnant, b, phi = np.broadcast_arrays(
        speed, mass, remnant, b, np.asarray(azimuth, dtype=float))
    direction = np.broadcast_to(direction, speed.shape+(3,))
    if not np.isfinite(phi).all():
        raise ValueError("azimuth must be finite")
    # Work with G M and b w^2 rather than b90 = GM/w^2: no zero-speed division.
    gm = G * (mass + remnant)
    hyp = np.hypot(gm, b * speed**2)
    a_ratio, b_ratio = gm / hyp, b * speed**2 / hyp
    cos_minus_one = -2 * a_ratio**2
    sine = 2 * a_ratio * b_ratio
    axis = np.eye(3)[np.argmin(np.abs(direction), axis=-1)]
    tangent = np.cross(direction, axis)
    norm = np.linalg.norm(tangent, axis=-1)
    tangent = np.divide(tangent, norm[..., None], out=np.zeros_like(tangent),
                        where=norm[..., None] > 0)
    other = np.cross(direction, tangent)
    azimuth_direction = np.cos(phi)[..., None] * tangent + np.sin(phi)[..., None] * other
    fraction = mass / (mass + remnant)
    return (fraction * speed)[..., None] * (
        cos_minus_one[..., None] * direction + sine[..., None] * azimuth_direction)


def _log_functions(z):
    """log(1+z)/z and log(1+z)/z - 1/(1+z), stable at small z."""
    z = np.asarray(z, dtype=float)
    small = z < 1.e-4
    safe = np.where(small, 1., z)
    first = np.log1p(safe) / safe
    second = first - 1 / (1 + safe)
    # The second expression starts at order z, so direct subtraction is poor
    # for very slow relative motion (b/b90 -> 0).
    series_z = np.where(small, z, 0.)
    first_series = 1 - series_z/2 + series_z**2/3 - series_z**3/4 + series_z**4/5 - series_z**5/6
    second_series = series_z/2 - 2*series_z**2/3 + 3*series_z**3/4 - 4*series_z**4/5 + 5*series_z**5/6
    return np.where(small, first_series, first), np.where(small, second_series, second)


def impact_moments(relative_velocity, field_mass, remnant_mass, G, b_min, b_max):
    """Integrate rate, increment and raw increment-square over b and azimuth.

    Returns ``rate, A, D`` per unit number density and natural time. It is exact
    for a finite impact-parameter annulus, including large deflections. The
    velocity distribution still has to be integrated. b_min/b_max may depend
    on speed and species, provided identical bounds define sampled encounters.
    """
    _, speed, direction = _vectors(relative_velocity)
    mass, remnant, lower, upper = _parameters(speed, field_mass, remnant_mass, G,
                                             b_min, b_max)
    gm = G * (mass + remnant)
    eta = mass / (mass + remnant)
    area_difference = (upper - lower) * (upper + lower)
    y0 = (lower * speed**2 / gm)**2
    y1 = (upper * speed**2 / gm)**2
    # Avoid subtracting nearly equal squared bounds or large b90 terms.
    z = area_difference * (speed**2 / gm)**2 / (1 + y0)
    logarithm, transverse_remainder = _log_functions(z)
    integrated_log = area_difference / (1 + y0) * logarithm
    A_parallel = -2*np.pi * eta * speed**2 * integrated_log
    D_parallel = 4*np.pi * eta**2 * speed**3 * area_difference / ((1+y0)*(1+y1))
    D_transverse = (2*np.pi * eta**2 * speed**3 * area_difference / (1+y0)
                    * (transverse_remainder + y0 / ((1+y0)*(1+z))))
    nn = direction[..., :, None] * direction[..., None, :]
    diffusion = D_parallel[..., None, None]*nn + D_transverse[..., None, None]*(np.eye(3)-nn)
    return {"rate": np.pi * speed * area_difference,
            "A": A_parallel[..., None] * direction, "D": diffusion}


def closest_approach_annulus(speed, field_mass, remnant_mass, G,
                             minimum_pericenter, b_max):
    """Convert a physical pericenter floor to an allowed impact annulus.

    b_min^2 = rp_min^2 + 2 GM rp_min/w^2. If this is at or beyond b_max,
    return the zero-width annulus [b_max,b_max] and an explicit empty mask.
    This excludes a supplied channel; it does not model its physical outcome.
    The zero-speed rate is zero irrespective of the convention for its bounds.
    """
    speed = np.asarray(speed, dtype=float)
    if not np.isfinite(speed).all() or np.any(speed < 0):
        raise ValueError("relative speeds must be finite and nonnegative")
    mass, remnant, rp, upper = _parameters(speed, field_mass, remnant_mass,
                                          G, 0., b_max)
    rp = np.broadcast_to(np.asarray(minimum_pericenter, dtype=float), mass.shape)
    if not np.isfinite(rp).all() or np.any(rp < 0):
        raise ValueError("minimum pericenter must be finite and nonnegative")
    focusing = np.divide(np.sqrt(2*G*(mass+remnant)*rp), speed,
                         out=np.full_like(mass, np.inf), where=speed > 0)
    focusing = np.where(rp == 0, 0., focusing)
    lower = np.hypot(rp, focusing)
    empty = lower >= upper
    return np.minimum(lower, upper), upper, empty


def velocity_ball(cutoff, order=20, radial_scale=None):
    """Positive volume quadrature on |u| < cutoff; radial axis is Cartesian z.

    Opposite nodes are exact pairs. The radial and polar rules have ``order``
    Gauss-Legendre nodes; the azimuth has ``2*order`` equally spaced nodes.
    Optional ``radial_scale`` uses u = scale sinh(t) to resolve a cold species
    without truncating its DF below the common escape-speed support.
    """
    if not np.isfinite(cutoff) or cutoff <= 0 or not isinstance(order, (int, np.integer)) or order < 4:
        raise ValueError("positive cutoff and integer quadrature order >= 4 required")
    x, weights = leggauss(order)
    phi = np.arange(2*order)*np.pi/order
    mu, phi_grid = np.meshgrid(x, phi, indexing="ij")
    directions = np.stack((np.sqrt(1-mu*mu)*np.cos(phi_grid),
                           np.sqrt(1-mu*mu)*np.sin(phi_grid), mu), axis=-1).reshape(-1, 3)
    indices = np.arange(len(directions)).reshape(order, 2*order)
    opposite = np.roll(indices[::-1], order, axis=1).ravel()
    first = np.flatnonzero(np.arange(len(directions)) < opposite)
    directions[opposite[first]] = -directions[first]
    if radial_scale is None:
        radius = cutoff*(x+1)/2
        radial_weight = cutoff/2*weights*radius**2
    else:
        if not np.isfinite(radial_scale) or radial_scale <= 0:
            raise ValueError("radial_scale must be finite and positive")
        end = np.arcsinh(cutoff/radial_scale)
        t = end*(x+1)/2
        radius = radial_scale*np.sinh(t)
        radial_weight = end/2*weights*radial_scale*np.cosh(t)*radius**2
    angular_weight = np.repeat(weights, 2*order)*np.pi/order
    return ((directions[:, None, :]*radius[None, :, None]).reshape(-1, 3),
            (angular_weight[:, None]*radial_weight[None, :]).ravel())


class DiscreteBath:
    """Positive velocity quadrature used identically for moments and events.

    ``number_weights[j,k]`` is the number density represented by velocity node
    k of species j, in pc^-3. This is a numerical discrete approximation to a
    local DF; convergence in velocity quadrature must be checked separately.
    Bounds are scalars, species-by-node arrays, or callables ``b(speed,mass)``.
    Sampled events are independent draws from this bath, not persistent stars.
    """

    def __init__(self, velocities, number_weights, masses, remnant_mass, G,
                 b_min, b_max):
        velocities = np.asarray(velocities, dtype=float)
        number_weights = np.asarray(number_weights, dtype=float)
        masses = np.asarray(masses, dtype=float)
        if (velocities.ndim != 2 or velocities.shape[1] != 3 or
                masses.ndim != 1 or number_weights.shape != (len(masses), len(velocities)) or
                not np.isfinite(velocities).all() or not np.isfinite(number_weights).all() or
                np.any(number_weights < 0)):
            raise ValueError("finite velocity nodes and nonnegative species-by-node number weights required")
        _parameters(np.zeros(len(masses)), masses, remnant_mass, G, 0., 0.)
        self.velocities, self.number_weights, self.masses = velocities, number_weights, masses
        self.remnant_mass, self.G, self.b_min, self.b_max = remnant_mass, G, b_min, b_max

    def _state(self, velocity):
        velocity = np.asarray(velocity, dtype=float)
        if velocity.shape != (3,) or not np.isfinite(velocity).all():
            raise ValueError("velocity must be one finite three-vector")
        relative = velocity[None, :] - self.velocities
        speed = np.linalg.norm(relative, axis=-1)[None, :]
        masses = self.masses[:, None]
        shape = self.number_weights.shape
        low = self.b_min(speed, masses) if callable(self.b_min) else self.b_min
        high = self.b_max(speed, masses) if callable(self.b_max) else self.b_max
        low, high = np.broadcast_to(low, shape), np.broadcast_to(high, shape)
        # An empty physical domain must be supplied as equal bounds; no silent
        # clipping or reordering of invalid cutoffs changes the encounter law.
        _parameters(np.broadcast_to(speed, shape), masses, self.remnant_mass, self.G, low, high)
        return np.broadcast_to(relative, shape+(3,)), low, high

    def moments(self, velocity):
        relative, low, high = self._state(velocity)
        moments = impact_moments(relative, self.masses[:, None], self.remnant_mass,
                                  self.G, low, high)
        weights = self.number_weights / PC_PER_KMS_MYR
        rate_species = np.sum(weights*moments["rate"], axis=1)
        A_species = np.einsum("jk,jka->ja", weights, moments["A"])
        D_species = np.einsum("jk,jkab->jab", weights, moments["D"])
        return {"rate": float(rate_species.sum()), "A": A_species.sum(axis=0),
                "D": D_species.sum(axis=0), "rate_species": rate_species,
                "A_species": A_species, "D_species": D_species}

    def rates_many(self, velocities, batch_size=16):
        """Per-Myr total event rates without allocating moment tensors."""
        return self._many(velocities, batch_size, rates_only=True)["rate"]

    def moments_many(self, velocities, batch_size=16):
        """Chunked matched moments for states sharing the same local bath."""
        return self._many(velocities, batch_size, rates_only=False)

    def _many(self, velocities, batch_size, rates_only):
        velocities = np.asarray(velocities, dtype=float)
        if velocities.ndim != 2 or velocities.shape[1] != 3 or not np.isfinite(velocities).all():
            raise ValueError("velocities must be finite N by three vectors")
        if not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
            raise ValueError("batch_size must be a positive integer")
        count, species, nodes = len(velocities), len(self.masses), len(self.velocities)
        output = {"rate": np.zeros(count), "rate_species": np.zeros((count, species))}
        if not rates_only:
            output.update(A=np.zeros((count, 3)), D=np.zeros((count, 3, 3)),
                          A_species=np.zeros((count, species, 3)),
                          D_species=np.zeros((count, species, 3, 3)))
        weights = self.number_weights / PC_PER_KMS_MYR
        mass = self.masses[None, :, None]
        for start in range(0, count, batch_size):
            end = min(start+batch_size, count)
            relative = velocities[start:end, None, None, :] - self.velocities[None, None, :, :]
            speed = np.linalg.norm(relative, axis=-1)
            shape = (end-start, species, nodes)
            relative = np.broadcast_to(relative, shape+(3,))
            low = self.b_min(speed, mass) if callable(self.b_min) else self.b_min
            high = self.b_max(speed, mass) if callable(self.b_max) else self.b_max
            low, high = np.broadcast_to(low, shape), np.broadcast_to(high, shape)
            _parameters(np.broadcast_to(speed, shape), mass, self.remnant_mass, self.G, low, high)
            if rates_only:
                rate = np.pi * speed * (high-low)*(high+low)
            else:
                moments = impact_moments(relative, mass, self.remnant_mass, self.G, low, high)
                rate = moments["rate"]
                output["A_species"][start:end] = np.einsum("jk,bjka->bja", weights, moments["A"])
                output["D_species"][start:end] = np.einsum("jk,qjkab->qjab", weights, moments["D"])
                output["A"][start:end] = output["A_species"][start:end].sum(axis=1)
                output["D"][start:end] = output["D_species"][start:end].sum(axis=1)
            output["rate_species"][start:end] = np.einsum("jk,bjk->bj", weights, rate)
            output["rate"][start:end] = output["rate_species"][start:end].sum(axis=1)
        return output

    def sample_events(self, velocity, rng, count=1):
        """Sample conditional-on-event increments at a fixed state.

        Waiting times/evolution are the caller's responsibility. Sampling b^2
        uniformly includes the cross-section measure b db; velocity nodes are
        selected using flux, species abundance and the allowed annular area.
        """
        if not isinstance(count, (int, np.integer)) or count < 0:
            raise ValueError("count must be a nonnegative integer")
        relative, low, high = self._state(velocity)
        speeds = np.linalg.norm(relative, axis=-1)
        rates = self.number_weights * speeds * np.pi*(high-low)*(high+low)
        total = rates.sum()
        if total <= 0:
            if count == 0:
                return np.empty((0, 3))
            raise ValueError("cannot sample an event from a zero-rate bath")
        selected = rng.choice(rates.size, size=count, p=rates.ravel()/total)
        species, nodes = np.unravel_index(selected, rates.shape)
        lo, hi = low[species, nodes], high[species, nodes]
        impact = np.sqrt(lo*lo + rng.random(count)*(hi-lo)*(hi+lo))
        phi = rng.uniform(0, 2*np.pi, count)
        return scatter_impulse(relative[species, nodes], self.masses[species],
                               self.remnant_mass, self.G, impact, phi)


def local_moments(number_df, velocity, cutoff, masses, remnant_mass, G,
                  b_min, b_max, order=20, radial_scale=None):
    """Integrate a supplied number DF using a positive spherical quadrature."""
    nodes, weights = velocity_ball(cutoff, order, radial_scale=radial_scale)
    df = np.asarray(number_df(nodes), dtype=float)
    if df.shape != (len(masses), len(nodes)):
        raise ValueError("number_df must return species-by-velocity values")
    bath = DiscreteBath(nodes, df*weights[None, :], masses, remnant_mass, G, b_min, b_max)
    return bath.moments(velocity)
