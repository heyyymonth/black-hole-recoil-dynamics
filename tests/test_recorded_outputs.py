"""Independent checks of reported numbers, using stored inputs and raw sources."""
import csv
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def read_rows(name):
    with (ROOT/name).open() as stream:return list(csv.DictReader(stream))


def test_all_recorded_collision_tensors_and_species_accounting():
    states=read_rows('results/local-transport-kernels/states.csv')
    tensors=json.loads((ROOT/'results/local-transport-kernels/tensors.json').read_text())
    species=read_rows('results/local-transport-kernels/species.csv')
    for state,tensor in zip(states,tensors,strict=True):
        D=np.array(tensor['D_over_lnLambda']);A=np.array(tensor['A_over_lnLambda'])
        axis=np.array([0,0,1]) if state['orientation']=='radial' else np.array([1,0,0])
        np.testing.assert_allclose(D,D.T,atol=1e-14)
        assert np.linalg.eigvalsh(D).min()>=-1e-12
        np.testing.assert_allclose(np.linalg.eigvalsh(D).min(),float(state['minimum_D_eigenvalue']))
        np.testing.assert_allclose(axis@D@axis,float(state['D_parallel_over_lnLambda_kms2_Myr']))
        np.testing.assert_allclose((np.trace(D)-axis@D@axis)/2,float(state['D_perpendicular_mean_over_lnLambda_kms2_Myr']))
        np.testing.assert_allclose(-A@axis,float(state['drag_over_lnLambda_kms_Myr']))
        terms=[r for r in species if r['state']==state['state']]
        assert len(terms)==16 and len({r['bin'] for r in terms})==16
        np.testing.assert_allclose(sum(float(r['trace_D_over_lnLambda_kms2_Myr']) for r in terms),np.trace(D))


def test_host_paper_values_against_species_arrays():
    h=dict(np.load(ROOT/'results/host-reconstruction/host.npz',allow_pickle=False))
    s=json.loads((ROOT/'results/host-reconstruction/summary.json').read_text())['summary']
    bh=h['types']=='BH';physical=h['types']!='TR'
    np.testing.assert_allclose((h['Mj'][bh]/h['mj'][bh]).sum(),s['total_BH_expected_number'])
    np.testing.assert_allclose(h['Mj'][bh].sum(),s['total_BH_mass_Msun'])
    np.testing.assert_allclose(h['rhoj'][physical,0].sum(),s['central_density_Msun_pc3'])
    np.testing.assert_allclose(np.sqrt(2*h['psi'][0]),s['central_cutoff_kms'])
    np.testing.assert_allclose(np.sqrt(2*(h['psi'][0]+h['G']*h['M']/h['rt'])),s['central_escape_to_infinity_kms'])


def test_orbit_summaries_against_detected_events():
    summary=json.loads((ROOT/'results/orbit/summary.json').read_text())
    events=read_rows('results/orbit/events.csv')
    for case in summary['cases']:
        for kind in ['pericenter','apocenter']:
            selected=[e for e in events if e['case']==case['case'] and e['event']==kind]
            assert len(selected)==case[kind+'_count']
            np.testing.assert_allclose(np.mean([float(e['R_kpc']) for e in selected]),case[kind+'_mean_kpc'])
        lo,hi=case['pericenter_mean_kpc'],case['apocenter_mean_kpc']
        np.testing.assert_allclose((hi-lo)/(hi+lo),case['spherical_eccentricity_from_mean_extrema'])


def test_catalog_display_scales_and_all_hst_conversions():
    catalog=read_rows('data/observed-clusters/derived/catalog_scalars.csv')
    for row in catalog:
        mass,radius=float(row['mass_msun']),float(row['half_mass_3d_pc'])
        velocity=np.sqrt(4.30091727003628e-3*mass/radius)
        np.testing.assert_allclose(velocity,float(row['scale_v_kms']),rtol=2e-10)
        np.testing.assert_allclose(.9777922216807892*radius/velocity,float(row['scale_t_myr']),rtol=2e-10)
    rows=read_rows('data/observed-clusters/derived/hst_native_and_common_units.csv')
    assert len(rows)==321
    assert len({(r['cluster'],r['source'],r['bin_id']) for r in rows})==107
    for row in rows:
        D,angle,pm=map(float,[row['distance_kpc'],row['radius_arcsec'],row['sigma_masyr']])
        radius=D*1000*np.deg2rad(angle/3600)
        speed=4.740470463533348*D*pm
        np.testing.assert_allclose(radius,float(row['radius_pc_at_catalog_distance']),rtol=2e-10)
        np.testing.assert_allclose(speed,float(row['sigma_kms_at_catalog_distance']),rtol=2e-10)
        np.testing.assert_allclose(radius/float(row['half_mass_3d_pc']),float(row['projected_radius_over_rh']),rtol=2e-10)
        for side in ['upper','lower']:
            np.testing.assert_allclose(4.740470463533348*D*float(row[f'error_{side}_masyr']),float(row[f'measurement_error_{side}_kms']),rtol=2e-10)
