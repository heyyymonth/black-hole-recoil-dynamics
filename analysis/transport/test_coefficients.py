"""Independent uniform-sphere analytic limits and covariance checks."""
import numpy as np
from coefficients import rosenbluth, PC_PER_KMS_MYR


def evaluate(v,order=32,df=None):
    n,U=7.,5.
    if df is None:
        df=lambda u: np.full((1,len(u)),3*n/(4*np.pi*U**3))
    return rosenbluth(df,v,U,[2.],60.,10.,.004302,order)


def test_uniform_sphere_rest_and_interior():
    n,U=7.,5.
    for speed in [0.,1.,3.]:
        v=np.array([speed,0,0]);r=evaluate(v)
        np.testing.assert_allclose(r['h'],n*(3*U**2-speed**2)/(2*U**3),rtol=1e-10)
        np.testing.assert_allclose(r['grad_h'][0],-n*v/U**3,atol=1e-13)
        expected=n/U*np.diag([1-3*speed**2/(5*U**2),1-speed**2/(5*U**2),1-speed**2/(5*U**2)])
        np.testing.assert_allclose(r['hess_g'][0],expected,rtol=1e-10,atol=1e-13)
        assert np.linalg.eigvalsh(r['D']).min()>0
        np.testing.assert_allclose(np.trace(r['D']),r['trace_expected'],rtol=1e-13)


def test_mass_prefactors_and_time_conversion():
    r=evaluate([1,0,0]);factor=4*np.pi*.004302**2*10/PC_PER_KMS_MYR
    np.testing.assert_allclose(r['A'],factor*2*62*r['grad_h'][0])
    np.testing.assert_allclose(r['D'],factor*4*r['hess_g'][0])


def test_anisotropic_refinement_and_rotation():
    v=np.array([1.2,.4,.8]);axis=np.array([0.,0.,1.])
    df=lambda u: np.exp(-np.sum(u*u,axis=1)/8-(u@axis)**2/3)[None,:]
    a=evaluate(v,24,df);b=evaluate(v,40,df)
    np.testing.assert_allclose(a['A'],b['A'],rtol=2e-6,atol=1e-10)
    np.testing.assert_allclose(a['D'],b['D'],rtol=2e-6,atol=1e-10)
    angle=.63;R=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
    c=evaluate(R@v,40,lambda u: df(u@R))
    np.testing.assert_allclose(c['A'],R@b['A'],rtol=2e-6,atol=1e-10)
    np.testing.assert_allclose(c['D'],R@b['D']@R.T,rtol=2e-6,atol=1e-10)


def test_exterior_velocity_is_not_silently_accepted():
    import pytest
    with pytest.raises(ValueError,match='Exterior-velocity'):
        evaluate([10.,0,0])
