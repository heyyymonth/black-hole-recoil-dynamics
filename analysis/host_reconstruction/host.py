"""Numerical host interface shared by independent validation and collisions."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import gammainc
from numpy.polynomial.legendre import leggauss

class PositiveProfile:
    """Interpolate positive densities/pressures in log space; enforce zero at rt.

    The final positive mesh point and asymptotic lowering exponent determine
    only the final, density-negligible interval. This avoids negative densities
    and catastrophic relative errors for segregated species with long tails.
    """
    def __init__(self,r,y,power):
        self.rt=r[-1];self.last=r[-2];self.y_last=y[:,-2];self.power=power
        if np.any(y[:,:-1]<=0):raise ValueError('Positive interior profile required')
        self.log_spline=CubicSpline(r[:-1],np.log(y[:,:-1]),axis=1)
    def __call__(self,r):
        r=np.asarray(r);base=np.exp(self.log_spline(np.minimum(r,self.last)))
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
        return np.maximum(self.psi_spline(np.clip(r,0,self.rt)),0)
    def number_df(self,r,u):
        u=np.asarray(u);speed2=np.sum(u*u,axis=-1)
        E=(self.psi(r)-speed2[None,:]/2)/self.s2j[:,None]
        positive=np.maximum(E,0)
        energy=np.exp(positive)*gammainc(self.g,positive)
        J2=r*r*(u[:,0]**2+u[:,1]**2)
        mass_df=self.A[:,None]*energy*np.exp(-J2[None,:]/(2*self.raj[:,None]**2*self.s2j[:,None]))
        return mass_df/self.mj[:,None]
    def project(self,R,order):
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
