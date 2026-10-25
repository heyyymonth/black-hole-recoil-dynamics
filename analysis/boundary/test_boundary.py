"""Independent analytic limits and input-derived boundary consistency checks."""
from pathlib import Path
import json, sys
import numpy as np
import pytest
from numpy.testing import assert_allclose
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'analysis/host_reconstruction'))
sys.path.insert(0, str(ROOT/'analysis/boundary'))
from host import Host
from model import EnclosedHost, diagnostics, differential_acceleration, jacobi_radius, max_force_ratio, tidal_tensor


def test_point_mass_and_extended_plummer_saddles():
    # Independent closed-form force-balance radii, including the extended-mass correction.
    G, M, lam, a = .004302, 1.8e5, .0012, 10.
    point = (G*M/lam)**(1/3)
    assert_allclose(jacobi_radius(lambda r: M, G, lam, .001, 1000), point, rtol=1e-12)
    enclosed = lambda r: M*r**3/(r*r+a*a)**1.5
    expected = np.sqrt((G*M/lam)**(2/3)-a*a)
    assert_allclose(jacobi_radius(enclosed, G, lam, .001, 1000), expected, rtol=1e-12)
    with pytest.raises(ValueError):
        jacobi_radius(enclosed, G, -lam, .001, 1000)


def test_exact_force_bound_and_frame_tensor():
    R, vc, r = 8930., 220., 84.
    # Sample directions and verify the analytic angular maximum is attained inward.
    rng = np.random.default_rng(438)
    n = rng.normal(size=(1000,3));n /= np.linalg.norm(n,axis=1)[:,None]
    n = np.vstack([n,[-1.,0.,0.],[1.,0.,0.]])
    acc = differential_acceleration([R,0,0],r*n,vc)
    exact_max = vc**2*r/(R*(R-r))
    assert np.max(np.linalg.norm(acc,axis=1)) <= exact_max*(1+1e-12)
    assert_allclose(np.linalg.norm(acc[-2]),exact_max,rtol=1e-12)
    # Check the tensor by differences of exact Galactic forces; centrifugal force is absent.
    step = .01
    fd = np.column_stack([(differential_acceleration([R,0,0],step*e,vc)-differential_acceleration([R,0,0],-step*e,vc))/(2*step) for e in np.eye(3)])
    assert_allclose(fd,tidal_tensor([R,0,0],vc),rtol=2e-9,atol=1e-12)
    assert_allclose(np.linalg.eigvalsh(fd),vc**2/R**2*np.array([-1,-1,1]),rtol=2e-9)


def test_recorded_host_boundaries_and_interior_domain():
    host = Host(dict(np.load(ROOT/'results/host-reconstruction/host.npz',allow_pickle=False)))
    model = EnclosedHost(host)
    tiny = 1e-7*host.rh
    assert_allclose(model.enclosed(tiny)/tiny**3, 4*np.pi/3*host.rhoj[host.physical,0].sum(), rtol=1e-7)
    recorded = json.loads((ROOT/'results/outer-boundary/summary.json').read_text())['baseline']
    result = diagnostics(model,8930.,220.,[.01,.03,.1])
    for key,value in result.items():
        if 'residual' not in key:
            assert_allclose(value,recorded[key],rtol=2e-9)
    radius = result['r_valid_0.01_pc']
    grid = np.geomspace(1e-3*host.rh,radius,500)
    eta = max_force_ratio(grid,model.enclosed,host.G,8930.,220.)
    assert np.all(np.diff(eta)>0) and max(eta)<.01*(1+1e-9)
    assert result['exact_inward_saddle_pc'] < result['rJ_extended_pc'] < result['exact_outward_saddle_pc']
    assert result['rJ_extended_pc'] < result['rJ_point_mass_pc']
    # Potential and force use the SAME mass inventory at the opening radius.
    r = result['rJ_extended_pc']; h=.001
    gradient=(model.potential(r+h)-model.potential(r-h))/(2*h)
    assert_allclose(gradient,model.G*model.enclosed(r)/r**2,rtol=1e-4)
    with pytest.raises(ValueError):
        max_force_ratio(8930,model.enclosed,host.G,8930.,220.)
