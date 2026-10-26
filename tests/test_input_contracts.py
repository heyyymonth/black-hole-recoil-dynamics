"""Regression checks for invalid-domain handling and declared model normalization."""
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for folder in ['host_reconstruction', 'boundary', 'orbit']:
    sys.path.insert(0, str(ROOT / 'analysis' / folder))
from host import Host
from galaxy import Galaxy
from model import EnclosedHost
from orbits import host_diagnostics, initial_state

@pytest.fixture
def host():
    return Host(dict(np.load(ROOT/'results/host-reconstruction/host.npz', allow_pickle=False)))

@pytest.mark.parametrize('radius', [-1., np.nan, np.inf, -np.inf])
def test_invalid_radii_are_rejected_consistently(host, radius):
    for evaluator in [host.psi, host.rho_spline, host.pr_spline, host.pt_spline,
                      EnclosedHost(host).enclosed, EnclosedHost(host).potential]:
        with pytest.raises(ValueError):
            evaluator(radius)
    with pytest.raises(ValueError):
        host.number_df(radius, np.zeros((1, 3)))
    with pytest.raises(ValueError):
        host.project(radius, 16)


def test_exterior_distribution_and_undefined_projected_moments(host):
    np.testing.assert_array_equal(host.rho_spline(2*host.rt), 0)
    np.testing.assert_array_equal(host.number_df(2*host.rt, np.zeros((2,3))), 0)
    for radius in [host.rt, 2*host.rt]:
        with pytest.raises(ValueError):
            host.project(radius, 16)
    for velocity in [np.ones(3), [[np.nan, 0, 0]], np.ones((2,2))]:
        with pytest.raises(ValueError):
            host.number_df(0, velocity)


def test_galactic_normalization_is_taken_from_configuration():
    p=json.loads((ROOT/'analysis/orbit/protocol.json').read_text())['galaxy']
    doubled={**p, 'potential_unit_kms2':2*p['potential_unit_kms2']}
    a,b=Galaxy(p),Galaxy(doubled)
    x=[7.59,-4.65,.71]
    for expected,actual in zip(a.fields(x),b.fields(x)):
        np.testing.assert_allclose(actual,2*expected)
    np.testing.assert_allclose(b.hessian_ball_bound(x,.03),2*a.hessian_ball_bound(x,.03))
    for bad in [np.nan,0.,-1.]:
        with pytest.raises(ValueError):Galaxy({**p,'bd_kpc':bad})


def test_diagnostic_labels_cannot_silently_change_tolerance(host):
    p=json.loads((ROOT/'analysis/orbit/protocol.json').read_text())
    s=json.loads((ROOT/p['source']).read_text())
    g=Galaxy(p['galaxy']);state=np.array([s['baseline_point_kpc_kms']])
    with pytest.raises(ValueError):host_diagnostics(g,EnclosedHost(host),state,epsilon=.03)
    for distance in [0,-1,np.nan,np.inf]:
        with pytest.raises(ValueError):initial_state(s,distance)
    for x in [[np.nan,0,0],[1,2]]:
        with pytest.raises(ValueError):g.fields(x)


def test_collision_runner_rejects_mislabeled_normalization(monkeypatch):
    spec=importlib.util.spec_from_file_location('collision_runner', ROOT/'analysis/transport/run.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(ROOT/'analysis/transport'))
    spec.loader.exec_module(module)
    original=module.json.loads
    def changed(raw):
        value=original(raw)
        if value.get('experiment')=='local-transport-kernels-1':value['log_lambda']=2
        return value
    monkeypatch.setattr(module.json,'loads',changed)
    with pytest.raises(ValueError,match='per-unit'):
        module.run()


def test_missing_orbit_refinement_is_not_reported_as_zero_error(monkeypatch):
    spec=importlib.util.spec_from_file_location('orbit_verifier',ROOT/'analysis/orbit/verify.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    original=module.load
    def incomplete(name):
        value=original(name)
        if name=='results/orbit/summary.json':
            del value['cases'][0]['refinement_max_position_kpc']
        return value
    monkeypatch.setattr(module,'load',incomplete)
    with pytest.raises(AssertionError):module.validate()


def test_failed_collision_checks_raise_and_retain_failure_record(monkeypatch,tmp_path):
    import shutil
    spec=importlib.util.spec_from_file_location('failed_collision_runner',ROOT/'analysis/transport/run.py')
    module=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(ROOT/'analysis/transport'));spec.loader.exec_module(module)
    base=tmp_path/'results/host-reconstruction';base.mkdir(parents=True)
    for name in ['host.npz','summary.json']:
        shutil.copyfile(ROOT/'results/host-reconstruction'/name,base/name)
    monkeypatch.setattr(module,'ROOT',tmp_path)
    monkeypatch.setattr(module,'OUT',tmp_path/'results/local-transport-kernels')
    original=module.json.loads
    def force_failure(raw):
        value=original(raw)
        if value.get('experiment')=='local-transport-kernels-1':
            value['quadrature_orders']=[4,4]
            value['acceptance']['relative_vector_and_tensor_refinement']=-1
        return value
    monkeypatch.setattr(module.json,'loads',force_failure)
    with pytest.raises(AssertionError,match='Collision validation failed'):module.run()
    record=original((module.OUT/'summary.json').read_text())
    assert record['all_numerical_checks_passed'] is False
