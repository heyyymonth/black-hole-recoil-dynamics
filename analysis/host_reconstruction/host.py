"""Numerical host interface shared by independent validation and collisions."""
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.special import gammainc
from numpy.polynomial.legendre import leggauss

def checked_radius(value):
    radius = np.asarray(value, dtype=float)
    if not np.isfinite(radius).all() or np.any(radius < 0):
        raise ValueError('radius must be finite and nonnegative')
    return radius


class PositiveProfile:
    """Interpolate positive densities/pressures in log space; enforce zero at rt.

    The final positive mesh point and asymptotic lowering exponent determine
    only the final, density-negligible interval. Intervals violating the nodal trend use a shape-preserving fallback.
    This avoids negative densities, artificial outer shells,
    and catastrophic relative errors for segregated species with long tails.
    """
    def __init__(self,r,y,power):
        self.rt=r[-1];self.last=r[-2];self.y_last=y[:,-2];self.power=power
        if np.any(y[:,:-1]<=0):raise ValueError('Positive interior profile required')
        self.knots=r[:-1]
        log_y=np.log(y[:,:-1])
        self.log_spline=CubicSpline(self.knots,log_y,axis=1)
        self.monotone_spline=PchipInterpolator(self.knots,log_y,axis=1)
        # A positive log-cubic can still invent enormous outer-density spikes
        # on irregular meshes. Inspect the derivative polynomial on EVERY
        # interval, including its interior extremum. Use a shape-preserving
        # interpolant throughout any interval that reverses its nodal trend.
        # Both interpolants agree at nodes, so the joined profile is continuous.
        a,b,c=self.log_spline.c[:3].transpose(0,2,1)
        width=np.diff(self.knots)[None,:]
        slope=np.diff(log_y,axis=1)/width
        vertex=np.divide(-b,3*a,out=np.zeros_like(a),where=a!=0)
        dv=np.where((vertex>0)&(vertex<width),3*a*vertex**2+2*b*vertex+c,c)
        d1=3*a*width**2+2*b*width+c
        direction=np.sign(slope)
        tolerance=1e-10*np.maximum(np.abs(slope),1/self.rt)
        self.shape_fallback=(np.minimum.reduce([direction*c,direction*d1,direction*dv]) < -tolerance)
        self.shape_fallback |= (direction==0)&(np.maximum.reduce([abs(c),abs(d1),abs(dv)])>tolerance)
    def __call__(self,r):
        r=checked_radius(r);clipped=np.minimum(r,self.last)
        interval=np.clip(np.searchsorted(self.knots,clipped,side='right')-1,0,len(self.knots)-2)
        logarithm=np.where(self.shape_fallback[:,interval],self.monotone_spline(clipped),self.log_spline(clipped))
        base=np.exp(logarithm)
        taper=np.maximum((self.rt-r)/(self.rt-self.last),0)**self.power
        return base*np.minimum(taper,1)

class Host:
    """Read-only numerical host in pc, km/s and Msun; DF is NUMBER normalized."""
    def __init__(self,arrays):
        self.data=arrays
        for key in ['r','mj','Mj','rhoj','mcj','s2j','raj','A','v2rj','v2tj','types']:
            setattr(self,key,arrays[key])
        self.G=float(arrays['G']);self.M=float(arrays['M']);self.rh=float(arrays['rh']);self.rt=float(arrays['rt'])
        self.psi_spline=CubicSpline(self.r,arrays['psi'])
        self.rho_spline=PositiveProfile(self.r,self.rhoj,float(arrays['g'])+1.5)
        self.pr_spline=PositiveProfile(self.r,self.rhoj*self.v2rj,float(arrays['g'])+2.5)
        self.pt_spline=PositiveProfile(self.r,self.rhoj*self.v2tj,float(arrays['g'])+2.5)
        self.g=float(arrays['g'])
        self.physical=self.types!='TR'
    def psi(self,r):
        r = checked_radius(r)
        return np.maximum(self.psi_spline(np.minimum(r,self.rt)),0)
    def number_df(self,r,u):
        r = checked_radius(r)
        u = np.asarray(u, dtype=float)
        if r.ndim != 0 or u.ndim != 2 or u.shape[1] != 3 or not np.isfinite(u).all():
            raise ValueError('number_df requires a scalar radius and finite N by 3 velocities')
        speed2=np.sum(u*u,axis=-1)
        E=(self.psi(r)-speed2[None,:]/2)/self.s2j[:,None]
        positive=np.maximum(E,0)
        energy=np.exp(positive)*gammainc(self.g,positive)
        J2=r*r*(u[:,0]**2+u[:,1]**2)
        mass_df=self.A[:,None]*energy*np.exp(-J2[None,:]/(2*self.raj[:,None]**2*self.s2j[:,None]))
        return mass_df/self.mj[:,None]
    def project(self,R,order):
        R = checked_radius(R)
        if R.ndim != 0 or R >= self.rt:
            raise ValueError('projected moments require 0 <= R < rt')
        x,w=leggauss(order);scale=.1*self.rh
        end=np.arcsinh(np.sqrt(self.rt*self.rt-R*R)/scale)
        t=(x+1)*end/2;z=scale*np.sinh(t);weight=w*end/2*scale*np.cosh(t)
        r=np.sqrt(R*R+z*z)
        rho=self.rho_spline(r);pr=self.pr_spline(r);pt=self.pt_spline(r)/2
        f=R*R/np.maximum(r*r,1e-300)
        sigma=2*(rho@weight)
        los=2*((pr*(1-f)+pt*f)@weight)/sigma
        pmr=2*((pr*f+pt*(1-f))@weight)/sigma
        pmt=2*(pt@weight)/sigma
        return np.array([sigma,los,pmr,pmt])
