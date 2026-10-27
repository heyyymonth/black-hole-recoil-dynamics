"""Physical limits, sampling semantics and evidence consistency for six hosts."""
from pathlib import Path
import importlib.util
import json
import sys
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]

def module(name):
    spec=importlib.util.spec_from_file_location('comparison_'+name,ROOT/'analysis/cluster_comparison'/f'{name}.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod


def test_joint_selection_is_weighted_reproducible_and_preserves_zero_weight():
    m=module('run');w,idx=m.selection(np.array([-np.inf,0.,np.log(3)]),0,4000)
    np.testing.assert_allclose(w,[0,.25,.75])
    np.testing.assert_array_equal(idx,m.selection(np.array([-np.inf,0.,np.log(3)]),0,4000)[1])
    assert not np.any(idx==0)
    assert abs(np.mean(idx==2)-.75)<.025
    with pytest.raises(ValueError):m.selection(np.array([-np.inf,-np.inf]),0)


def test_zero_speed_integral_uniform_ball_and_species_partition():
    m=module('run')
    class Uniform:
        G=.004302; mj=np.array([.5,10.]);physical=np.array([True,True])
        def psi(self,r):return 8.
        def number_df(self,r,u):return np.ones((2,len(u)))*np.array([2.,.01])[:,None]
    host=Uniform()
    expected=16*np.pi**2*host.G**2*host.mj**2*np.array([2.,.01])*16/.9777922216807892
    np.testing.assert_allclose(m.zero_speed(host,0,32),expected,rtol=1e-12)


def test_absent_bh_not_imputed_and_invalid_reservoir_rejected():
    m=module('posterior_inventory')
    r=m.summarize([0,2],[0,20],[np.nan,.3])
    assert r['zero_BH_fraction']==.5
    assert r['probe_to_BH_mass_ratio_given_present']['median']==3
    assert r['BH_half_mass_radius_pc_given_present']['median']==.3
    r=m.summarize([0,0],[0,0],[np.nan,123.])
    assert r['probe_to_BH_mass_ratio_given_present'] is None
    with pytest.raises(ValueError):m.summarize([1],[0],[.3])


def test_pinned_orbit_row_extraction_matches_ngc3201_anchor():
    m=module('environment')
    anchor=json.loads((ROOT/'data/orbit/ngc3201-source.json').read_text())
    rawpath=ROOT/anchor['raw_path']
    if rawpath.exists():raw=rawpath.read_text()
    else:raw=(ROOT/'data/cluster-comparison-environment/NGC3201-row.html').read_text()
    parsed,_=m.extract_source(raw,'NGC3201',anchor)
    for key in ['baseline_point_kpc_kms','distance_kpc','distance_error_kpc','pm_correlation','radial_velocity_kms']:
        assert parsed[key]==anchor[key]
    with pytest.raises(ValueError):m.extract_source(raw,'NGC9999',anchor)


def test_environment_sources_work_without_untracked_full_page(tmp_path,monkeypatch):
    import shutil
    m=module('environment')
    data=tmp_path/'data/cluster-comparison-environment'
    out=tmp_path/'results/cluster-comparison-environment'
    shutil.copytree(m.DATA,data)
    out.mkdir(parents=True)
    shutil.copy(m.OUT/'provenance.json',out/'provenance.json')
    anchor=tmp_path/'data/orbit/ngc3201-source.json'
    anchor.parent.mkdir(parents=True)
    shutil.copy(ROOT/'data/orbit/ngc3201-source.json',anchor)
    monkeypatch.setattr(m,'ROOT',tmp_path)
    monkeypatch.setattr(m,'DATA',data)
    monkeypatch.setattr(m,'OUT',out)
    sources=m.prepare_sources()
    assert len(sources)==6
    for cluster,source in sources.items():
        assert source==json.loads((data/f'{cluster}.json').read_text())
    row=data/'NGC6121-row.html'
    row.write_text(row.read_text().replace('6.42','6.43'))
    with pytest.raises(ValueError,match='extract changed'):m.prepare_sources()


def test_static_and_orbit_diagnostics_reproduce_original_ngc3201():
    m=module('environment')
    host=m.EnclosedHost(m.Host(dict(np.load(ROOT/'results/host-reconstruction/host.npz'))))
    old=json.loads((ROOT/'results/outer-boundary/summary.json').read_text())['baseline']
    new=m.static_diagnostics(host,old['R_pc'],old['vc_kms'],[.01,.03,.1])
    for key in ['rJ_extended_pc','r_valid_0.01_pc','exact_inward_saddle_pc','exact_outward_saddle_pc']:
        np.testing.assert_allclose(new[key],old[key],rtol=1e-8)
    # The point-mass continuation must allow an opening outside compact DF support.
    source=json.loads((ROOT/'data/orbit/ngc3201-source.json').read_text())
    p=json.loads((ROOT/'analysis/orbit/protocol.json').read_text());galaxy=m.Galaxy(p['galaxy'])
    from orbits import host_diagnostics as original
    state=m.initial_state(source,source['distance_kpc'])[None,:]
    old=original(galaxy,host,state);new=m.host_diagnostics(galaxy,host,state,.001)
    for key in old:np.testing.assert_allclose(new[key],old[key],rtol=1e-8)


def test_outer_profile_cannot_invent_a_shell_on_irregular_nodes():
    sys.path.insert(0,str(ROOT/'analysis/host_reconstruction'))
    from host import PositiveProfile
    r=np.array([0.,1.,10.,99.,99.9,99.99,100.])
    y=((100-r)**2)[None,:]
    p=PositiveProfile(r,y,2.)
    query=np.linspace(0,100,30001)
    interpolated=p(query)[0]
    assert np.all(np.diff(interpolated)<=1e-8)
    assert 0<=interpolated.min()<=interpolated.max()<=y.max()*(1+1e-12)
    np.testing.assert_allclose(p(r),y,atol=1e-10,rtol=1e-12)


def test_dense_even_df_has_zero_drift_but_streaming_df_does_not():
    sys.path.insert(0,str(ROOT/'analysis/transport'))
    from coefficients import rosenbluth
    even=lambda u: 1e9*np.exp(-np.sum(u*u,axis=1))[None,:]
    streaming=lambda u:even(u)*(1+.1*u[:,2])[None,:]
    a=rosenbluth(even,[0,0,0],3,[10],60,1,.004302,24)
    np.testing.assert_array_equal(a['A'],[0,0,0])
    b=rosenbluth(streaming,[0,0,0],3,[10],60,1,.004302,24)
    assert b['A'][2]>0
