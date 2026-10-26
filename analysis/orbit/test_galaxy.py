"""Independent force, tensor, bound, coordinate and orbit checks."""
from pathlib import Path
import json,sys
import numpy as np
import pytest
from numpy.testing import assert_allclose
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis/orbit'))
from galaxy import Galaxy
from orbits import ICRS_TO_CATALOG,initial_state,integrate,conserved,angular_fields,K
P=json.loads((ROOT/'analysis/orbit/protocol.json').read_text())


def test_force_is_potential_gradient_and_tensor_is_force_gradient():
    g=Galaxy(P['galaxy'])
    for x in [np.array([8.4,0.,0.]),np.array([7.59,-4.65,.71]),np.array([20.,5.,-4.])]:
        h=1e-4
        fd=np.array([-(g.fields(x+h*e)[0]-g.fields(x-h*e)[0])/(2*h) for e in np.eye(3)])
        dt=np.column_stack([(g.fields(x+h*e)[1]-g.fields(x-h*e)[1])/(2*h) for e in np.eye(3)])
        _,a,t=g.fields(x)
        assert_allclose(fd,a,rtol=1e-7,atol=1e-6)
        assert_allclose(dt,t,rtol=2e-7,atol=2e-5)
        assert_allclose(t,t.T,atol=1e-12)
    assert_allclose(np.sqrt(-8.4*g.fields([8.4,0,0])[1][0]),242.,atol=.1)


def test_disk_symmetry_and_halo_cutoff_matching():
    g=Galaxy(P['galaxy'])
    _,a,t=g.fields([10.,0.,0.]);assert_allclose(a[1:],0,atol=1e-14)
    assert_allclose(t,np.diag(np.diag(t)),atol=1e-14)
    c=g.cutoff;h=1e-5
    left=g.fields([c-h,0,0]);right=g.fields([c+h,0,0])
    assert_allclose(left[0],right[0],rtol=2e-7)
    assert_allclose(left[1],right[1],rtol=3e-7,atol=1e-12)
    # The density, hence second derivative, may jump at the halo cutoff.
    with pytest.raises(ValueError):g.hessian_ball_bound([199.,0.,0.],2)


def test_analytic_ball_bound_dominates_exact_forces_and_hessians():
    rng=np.random.default_rng(183)
    g=Galaxy(P['galaxy'])
    for center in [np.array([8.4,0.,0.]),np.array([7.59,-4.65,.71]),np.array([20.,5.,-4.])]:
        for radius in [.005,.03,.15]:
            n=rng.normal(size=(500,3));n/=np.linalg.norm(n,axis=1)[:,None]
            offset=radius*n
            bound=g.hessian_ball_bound(center,radius)
            exact=g.fields(center+offset)[1]-g.fields(center)[1]
            assert np.max(np.linalg.norm(exact,axis=1))<=bound*radius*(1+1e-12)
            inside=offset*rng.uniform(size=(500,1))
            assert np.max(np.linalg.norm(g.fields(center+inside)[2],ord=2,axis=(1,2)))<=bound*(1+1e-12)


def test_distance_transport_and_sky_basis():
    s=json.loads((ROOT/P['source']).read_text())
    assert_allclose(ICRS_TO_CATALOG@ICRS_TO_CATALOG.T,np.eye(3),atol=2e-14)
    assert_allclose(initial_state(s,s['distance_kpc']),s['baseline_point_kpc_kms'],atol=0)
    delta=initial_state(s,s['distance_kpc']+1)-initial_state(s,s['distance_kpc'])
    assert_allclose(np.linalg.norm(delta[:3]),1,rtol=1e-14)
    assert_allclose(delta[:3]@delta[3:],0,atol=1e-12)
    assert_allclose(np.linalg.norm(delta[3:]),4.740470463533348*np.hypot(s['pmra_cosdec_masyr'],s['pmdec_masyr']),rtol=1e-14)
    # Rounded Galactic l,b from the original table are independent orientation checks.
    n=delta[:3];l=np.degrees(np.arctan2(n[1],-n[0]))%360;b=np.degrees(np.arcsin(n[2]))
    assert_allclose([l,b],[277.229,8.640],atol=.00051)


def test_short_orbit_conservation_refinement_and_reversal():
    g=Galaxy(P['galaxy']);s=json.loads((ROOT/P['source']).read_text());y=np.asarray(s['baseline_point_kpc_kms'])
    sol=integrate(g,y,[0,40],max_step=.5,rtol=1e-10)
    t=np.linspace(0,40,81);track=sol.sol(t).T;e,l=conserved(g,track)
    assert np.max(np.abs(e/e[0]-1))<1e-9
    assert np.max(np.abs(l/l[0]-1))<1e-9
    back=integrate(g,track[-1],[40,0],max_step=.25,rtol=1e-12)
    assert_allclose(back.y[:,-1],y,rtol=0,atol=1e-7)


def test_rotating_frame_euler_derivative_uses_physical_units():
    g=Galaxy(P['galaxy']);s=json.loads((ROOT/P['source']).read_text());y=np.asarray(s['baseline_point_kpc_kms'])
    dt=.01
    plus=integrate(g,y,[0,dt],max_step=.005,rtol=1e-12).y[:,-1]
    minus=integrate(g,y,[0,-dt],max_step=.005,rtol=1e-12).y[:,-1]
    omega_plus=angular_fields(plus[:3],plus[3:],g.fields(plus[:3])[1])[0]
    omega_minus=angular_fields(minus[:3],minus[3:],g.fields(minus[:3])[1])[0]
    physical_derivative=(omega_plus-omega_minus)/(2*dt*K)
    expected=angular_fields(y[:3],y[3:],g.fields(y[:3])[1])[1]
    assert_allclose(physical_derivative,expected,rtol=3e-6,atol=2e-6)
