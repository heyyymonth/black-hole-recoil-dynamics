"""Rosenbluth moments for bounded, possibly anisotropic NUMBER distributions.

Units: number DF in pc^-3 (km/s)^-3; masses in Msun; speeds in km/s;
returned A in (km/s)/Myr and D in (km/s)^2/Myr. Constant ln Lambda only.
The ray integral is centered on test velocity, so no inverse-power singularity
is evaluated numerically. Positive quadrature weights preserve PSD diffusion.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

PC_PER_KMS_MYR = 0.9777922216807892


def rosenbluth(number_df, velocity, cutoff, masses, remnant_mass, log_lambda,
               gravitational_constant, order=24):
    """number_df(u[N,3]) returns [species,N]; u is in cluster rest frame.

    The current validated API requires |velocity|<cutoff. Field support must
    be contained in the ball |u|<cutoff. Zero field streaming
    is not assumed by the integration. A single constant cutoff/log Lambda and
    a locally homogeneous weak-deflection closure are explicit inputs.
    """
    v = np.asarray(velocity, dtype=float)
    m = np.asarray(masses, dtype=float)
    if v.shape != (3,) or not np.isfinite(v).all():
        raise ValueError('velocity must be one finite Cartesian vector')
    scalars=np.asarray([cutoff,remnant_mass,log_lambda,gravitational_constant],dtype=float)
    if not np.isfinite(scalars).all() or np.any(scalars<=0):
        raise ValueError('positive physical inputs required')
    if m.ndim != 1 or np.any(m <= 0) or not np.isfinite(m).all() or order < 4:
        raise ValueError('positive species masses and quadrature order >=4 required')
    if np.linalg.norm(v)>=cutoff:
        raise ValueError('Exterior-velocity quadrature is not validated; an angular refinement scheme is required')
    mu, wmu = leggauss(order)
    phi = np.arange(2*order) * np.pi/order
    mu, phi = np.meshgrid(mu, phi, indexing='ij')
    directions = np.stack([np.sqrt(1-mu**2)*np.cos(phi),
                           np.sqrt(1-mu**2)*np.sin(phi), mu], axis=-1).reshape(-1,3)
    solid_weights = np.repeat(wmu, 2*order)*np.pi/order
    a = directions @ v
    discr = a*a + cutoff*cutoff - v@v
    good = discr > 0
    root = np.sqrt(np.maximum(discr,0))
    lo = np.maximum(0,a-root)
    hi = a+root
    good &= hi > lo
    directions, solid_weights = directions[good],solid_weights[good]
    lo, hi = lo[good],hi[good]
    x, wx = leggauss(order)
    s = lo[:,None]+(hi-lo)[:,None]*(x+1)/2
    weights = solid_weights[:,None]*(hi-lo)[:,None]*wx/2
    u = (v-directions[:,None,:]*s[:,:,None]).reshape(-1,3)
    f = np.asarray(number_df(u), dtype=float)
    if f.shape != (len(m),len(u)) or np.any(f<0) or not np.isfinite(f).all():
        raise ValueError('number DF must be finite, nonnegative, species by velocity')
    f = f.reshape(len(m),len(directions),order)
    # h = integral F s ds dOmega; grad h = -integral F n ds dOmega.
    h = np.einsum('jas,as,as->j',f,s,weights)
    radial0 = np.einsum('jas,as->ja',f,weights)
    radial1 = np.einsum('jas,as,as->ja',f,s,weights)
    grad_h = -radial0 @ directions
    transverse = np.eye(3)[None,:,:]-directions[:,:,None]*directions[:,None,:]
    hess_g = np.einsum('ja,abc->jbc',radial1,transverse)
    factor = 4*np.pi*gravitational_constant**2*log_lambda/PC_PER_KMS_MYR
    A_species = factor*(m*(remnant_mass+m))[:,None]*grad_h
    D_species = factor*m[:,None,None]**2*hess_g
    return {'A':A_species.sum(axis=0),'D':D_species.sum(axis=0),
            'A_species':A_species,'D_species':D_species,'h':h,
            'grad_h':grad_h,'hess_g':hess_g,
            'trace_expected':float(2*factor*np.sum(m*m*h))}
