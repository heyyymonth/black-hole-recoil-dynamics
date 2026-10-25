"""Independent consistency checks against the exported numerical host."""
from pathlib import Path
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simpson
sys.path.insert(0,str(Path(__file__).parent))
from host import Host
ROOT=Path(__file__).resolve().parents[2]

def make():return Host(dict(np.load(ROOT/'results/host-reconstruction/host.npz',allow_pickle=False)))

def test_all_arrays_finite_and_mass_partition():
    h=make()
    for k,a in h.data.items():
        if k!='types':assert np.isfinite(a).all(),k
    np.testing.assert_allclose(h.Mj.sum(),h.M,rtol=1e-12)
    r=np.r_[0,np.geomspace(1e-7,h.rt,8192)]
    integrated=4*np.pi*simpson(h.rho_spline(r).sum(axis=0)*r*r,x=r)
    np.testing.assert_allclose(integrated,h.M,rtol=1e-6)


def test_number_distribution_integrates_to_mass_density():
    h=make();r=.1*h.rh;U=np.sqrt(2*h.psi(r));x,w=leggauss(72)
    v=(x+1)*U/2;vv,mu=np.meshgrid(v,x,indexing='ij')
    u=np.column_stack([(vv*np.sqrt(1-mu*mu)).ravel(),np.zeros(vv.size),(vv*mu).ravel()])
    weights=(2*np.pi*(w*U/2*v*v)[:,None]*w).ravel()
    rho=(h.number_df(r,u)@weights)*h.mj
    np.testing.assert_allclose(rho,h.rho_spline(r),rtol=1e-6)
    assert np.all(h.rho_spline(np.array([0,.1*h.rh,h.rh,3*h.rh]))>0)
    np.testing.assert_array_equal(h.rho_spline(h.rt),np.zeros(len(h.mj)))


def test_cutoff_boundary_and_velocity_symmetry():
    h=make();r=h.rh;U=np.sqrt(2*h.psi(r))
    u=np.array([[0,0,.2*U],[0,0,-.2*U],[.2*U,0,0],[0,.2*U,0],[0,0,1.01*U]])
    f=h.number_df(r,u)
    np.testing.assert_allclose(f[:,0],f[:,1]);np.testing.assert_allclose(f[:,2],f[:,3])
    np.testing.assert_array_equal(f[:,-1],np.zeros(len(h.mj)))


def test_central_diffusion_against_independent_radial_integral():
    """At r=v=0 the isotropic host reduces to sixteen 1-D integrals."""
    from scipy.integrate import quad
    import csv
    h=make();U=float(np.sqrt(2*h.psi(0)))
    hj=np.array([4*np.pi*quad(lambda speed:float(h.number_df(0,np.array([[speed,0,0]]))[j,0])*speed,0,U,epsabs=1e-10,epsrel=1e-11)[0] for j in range(len(h.mj))])
    expected=8*np.pi*h.G*h.G*np.sum(h.mj**2*hj)/.9777922216807892
    with (ROOT/'results/local-transport-kernels/states.csv').open() as f:
        central=next(csv.DictReader(f))
    np.testing.assert_allclose(float(central['trace_D_over_lnLambda_kms2_Myr']),expected,rtol=1e-9)
