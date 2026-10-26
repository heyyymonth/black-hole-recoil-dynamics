"""Galactic reference motion and host diagnostics; never evolves a remnant."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from galaxy import Galaxy
# Standard J2000 equatorial -> Galactic orthogonal matrix; X reflection matches catalog.
ICRS_TO_CATALOG = np.array([
    [.0548755604162154,.8734370902348850,.4838350155487132],
    [.4941094278755837,-.4448296299600112,.7469822444972189],
    [-.8676661490190047,-.1980763734312015,.4559837761750669]])
K = .001022712165045695  # km/s -> kpc/Myr


def initial_state(source, distance):
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError('Heliocentric distance must be finite and positive')
    ra,dec=np.deg2rad([source['RA_deg'],source['DEC_deg']])
    radial=np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
    east=np.array([-np.sin(ra),np.cos(ra),0.])
    north=np.array([-np.sin(dec)*np.cos(ra),-np.sin(dec)*np.sin(ra),np.cos(dec)])
    tangential_per_kpc=4.740470463533348*(source['pmra_cosdec_masyr']*east+source['pmdec_masyr']*north)
    delta=distance-source['distance_kpc']
    anchor=np.asarray(source['baseline_point_kpc_kms'],dtype=float)
    return anchor+delta*np.r_[ICRS_TO_CATALOG@radial,ICRS_TO_CATALOG@tangential_per_kpc]


def integrate(galaxy,y0,interval,max_step=.5,rtol=1e-10,atol=1e-11):
    def rhs(t,y):return K*np.r_[y[3:],galaxy.fields(y[:3])[1]]
    def radial(t,y):return y[:3]@y[3:]
    def plane(t,y):return y[2]
    result=solve_ivp(rhs,interval,y0,method='DOP853',dense_output=True,events=[radial,plane],rtol=rtol,atol=atol,max_step=max_step)
    if not result.success:raise RuntimeError(result.message)
    return result


def conserved(galaxy,states):
    x,v=states[:,:3],states[:,3:]
    energy=galaxy.fields(x)[0]+.5*np.sum(v*v,axis=1)
    lz=x[:,0]*v[:,1]-x[:,1]*v[:,0]
    return energy,lz


def angular_fields(x,v,acc):
    """Minimal radial-frame rotation and its physical derivative.

    omega: km/s/kpc; omega_dot: (km/s)^2/kpc^2, not per Myr.
    """
    s2=np.sum(x*x,axis=-1)
    omega=np.cross(x,v)/s2[...,None]
    derivative=np.cross(x,acc)/s2[...,None]-2*np.sum(x*v,axis=-1)[...,None]*omega/s2[...,None]
    return omega,derivative


def host_diagnostics(galaxy,host,states,epsilon=.01,derivative_step=.001):
    if epsilon != .01:
        raise ValueError('The named 1-percent diagnostic requires epsilon=0.01')
    if not np.isfinite(derivative_step) or derivative_step <= 0:
        raise ValueError('Positive finite tensor derivative step required')
    x,v=states[:,:3],states[:,3:]
    tensor=galaxy.fields(x)[2]
    norm=np.max(np.abs(np.linalg.eigvalsh(tensor)),axis=1)/1e6 # pc tensor
    omega,omega_dot=angular_fields(x,v,galaxy.fields(x)[1])
    effective=tensor+np.sum(omega**2,axis=1)[:,None,None]*np.eye(3)-omega[:,:,None]*omega[:,None,:]
    lam=np.linalg.eigvalsh(effective)[:,-1]/1e6
    shift=K*v*derivative_step
    plus=galaxy.fields(x+shift)[2];minus=galaxy.fields(x-shift)[2]
    tdot=(plus-minus)/(2*derivative_step)
    amplitude_rate=np.abs(np.linalg.norm(plus,ord=2,axis=(1,2))-np.linalg.norm(minus,ord=2,axis=(1,2)))/(2*derivative_step)/1e6/norm
    variation=np.linalg.norm(tdot,ord=2,axis=(1,2))/1e6/norm # Myr^-1
    radii=[];scales=[]
    for point,l in zip(x,lam):
        def ratio(r):return galaxy.hessian_ball_bound(point,r/1000)/1e6*r**3/(host.G*host.enclosed(r))
        lower=1e-3*host.host.rh;upper=min(host.rt,.1*np.linalg.norm(point)*1000)
        if ratio(lower)>=epsilon:raise ValueError('No resolved central force-accuracy domain')
        radii.append(brentq(lambda r:ratio(r)-epsilon,lower,upper,xtol=1e-8))
        if l<=0:raise ValueError('No extensive frozen-frame diagnostic in this state')
        scales.append(brentq(lambda r:host.G*host.enclosed(r)/r**3-l,lower,upper,xtol=1e-8))
    radii=np.asarray(radii);scales=np.asarray(scales)
    def tdyn(r):return np.sqrt(r**3/(host.G*host.enclosed(r)))/(K*1000)
    rh=host.host.rh
    return {'tidal_norm_kms2_pc2':norm,'tidal_variation_per_Myr':variation,
            'tidal_amplitude_variation_per_Myr':amplitude_rate,
            'euler_over_extensive_eigenvalue':np.linalg.norm(omega_dot,axis=1)/(lam*1e6),
            'r_sufficient_1pct_pc':radii,'frozen_frame_scale_pc':scales,
            'mass_beyond_frozen_scale_fraction':1-host.enclosed(scales)/host.M,
            'eta_linear_rh':norm*rh**3/(host.G*host.enclosed(rh)),
            'eta_linear_3rh':norm*(3*rh)**3/(host.G*host.enclosed(3*rh)),
            'eta_bound_rh':galaxy.hessian_ball_bound(x,rh/1000)/1e6*rh**3/(host.G*host.enclosed(rh)),
            'eta_bound_3rh':galaxy.hessian_ball_bound(x,3*rh/1000)/1e6*(3*rh)**3/(host.G*host.enclosed(3*rh)),
            'timescale_ratio_rh':tdyn(rh)*variation,
            'amplitude_timescale_ratio_rh':tdyn(rh)*amplitude_rate,
            'timescale_ratio_3rh':tdyn(3*rh)*variation,
            'timescale_ratio_frozen_scale':tdyn(scales)*variation}
