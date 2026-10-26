"""Irrgang et al. (2013), Model I. Positions kpc, velocities km/s.

Potential: (km/s)^2; acceleration: (km/s)^2/kpc; tidal tensor: (km/s)^2/kpc^2.
The paper's potential unit is exactly 100 (km/s)^2; no rounded mass-unit conversion.
"""
import numpy as np

class Galaxy:
    def __init__(self, parameters):
        self.p = parameters
        positive = ['Mb_Mgal', 'Md_Mgal', 'Mh_Mgal', 'bb_kpc', 'ad_kpc',
                    'bd_kpc', 'ah_kpc', 'halo_cutoff_kpc', 'potential_unit_kms2']
        if any(not np.isfinite(parameters[k]) or parameters[k] <= 0 for k in positive):
            raise ValueError('Galactic model scales must be finite and positive')
        unit = parameters['potential_unit_kms2']
        self.cb = unit*parameters['Mb_Mgal']
        self.cd = unit*parameters['Md_Mgal']
        self.ch = unit*parameters['Mh_Mgal']
        self.bb = parameters['bb_kpc'];self.ad = parameters['ad_kpc'];self.bd = parameters['bd_kpc']
        self.ah = parameters['ah_kpc'];self.cutoff = parameters['halo_cutoff_kpc']
        if parameters['gamma'] != 2:raise ValueError('Only the published gamma=2 Model I is implemented')

    def fields(self, x):
        x = np.asarray(x,dtype=float)
        if x.shape[-1:] != (3,) or not np.isfinite(x).all():
            raise ValueError('Galactic positions must be finite Cartesian vectors')
        s = np.linalg.norm(x,axis=-1)
        if np.any(s<=0):raise ValueError('The halo cusp at the Galactic origin has no unique force direction')
        eye = np.eye(3);xx = x[..., :,None]*x[...,None,:]
        b = np.sqrt(s*s+self.bb**2)
        phi = -self.cb/b
        acc = -self.cb*x/b[...,None]**3
        tensor = self.cb*(3*xx/b[...,None,None]**5-eye/b[...,None,None]**3)
        z = x[...,2];zz = np.sqrt(z*z+self.bd**2);aa = self.ad+zz
        d = np.sqrt(x[...,0]**2+x[...,1]**2+aa*aa)
        q = x.copy();q[...,2] = aa*z/zz
        jj = np.broadcast_to(eye,x.shape[:-1]+(3,3)).copy();jj[...,2,2] += self.ad*self.bd**2/zz**3
        phi -= self.cd/d
        acc -= self.cd*q/d[...,None]**3
        tensor += self.cd*(3*q[..., :,None]*q[...,None,:]/d[...,None,None]**5-jj/d[...,None,None]**3)
        c = self.cutoff;a = self.ah
        kc = self.ch*c*c/(a*(a+c))
        inner = s<c
        ph = np.where(inner,self.ch/a*(np.log((a+s)/(a+c))-c/(a+c)),-kc/s)
        force = np.where(inner,self.ch/(a*(a+s)),kc/s**2)
        derivative = np.where(inner,-self.ch/(a*(a+s)**2),-2*kc/s**3)
        n = x/s[...,None]
        phi += ph
        acc -= force[...,None]*n
        tensor += -(force/s)[...,None,None]*eye+(force/s-derivative)[...,None,None]*n[..., :,None]*n[...,None,:]
        return phi,acc,tensor

    def hessian_ball_bound(self, x, radius_kpc):
        """Rigorous componentwise operator-norm bound throughout the spatial ball.

        Bulge: <=2 C/(s_min^2+b^2)^(3/2).
        Disk Hessian = C[J/d^3-3qq^T/d^5], J>=I, |q|<=d;
        its norm <= C max(2, Jzz_max)/d_min^3 by eigenvalue bounds.
        Halo gamma=2 norm = C/[a s(a+s)], decreasing with s inside cutoff.
        The summed bound is conservative; it is not a sampled angular maximum.
        """
        x = np.asarray(x,dtype=float);r = np.asarray(radius_kpc)
        if x.shape[-1:] != (3,) or not np.isfinite(x).all() or not np.isfinite(r).all():
            raise ValueError('Finite Cartesian centers and radii required')
        s = np.linalg.norm(x,axis=-1);smin = s-r
        if np.any(r<0) or np.any(smin<=0) or np.any(s+r>=self.cutoff):
            raise ValueError('Ball must avoid the origin and lie inside the halo cutoff')
        planar = np.hypot(x[...,0],x[...,1])
        zmin = np.maximum(np.abs(x[...,2])-r,0)
        zzmin = np.sqrt(zmin*zmin+self.bd**2)
        dmin = np.sqrt(np.maximum(planar-r,0)**2+(self.ad+zzmin)**2)
        jmax = 1+self.ad*self.bd**2/zzmin**3
        return (2*self.cb/(smin*smin+self.bb**2)**1.5
                +self.cd*np.maximum(2,jmax)/dmin**3
                +self.ch/(self.ah*smin*(self.ah+smin)))
