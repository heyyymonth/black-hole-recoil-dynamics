"""Verify orbit source extraction and recorded evidence; no integration or network."""
from pathlib import Path
import argparse,csv,hashlib,html,json,math,re
ROOT=Path(__file__).resolve().parents[2]
INPUTS=['analysis/orbit/'+n for n in ['galaxy.py','orbits.py','run.py','test_galaxy.py','protocol.json','verify.py']]+['data/orbit/ngc3201-source.json','data/orbit/ngc3201-row.html','results/host-reconstruction/host.npz','results/host-reconstruction/summary.json','results/calculations.json','analysis/boundary/model.py','analysis/host_reconstruction/host.py']
OUTPUTS=['results/orbit/'+n for n in ['summary.json','tracks.csv','events.csv']]
def hashes(paths):return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
def load(name):return json.loads((ROOT/name).read_text(),parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
def rows(name):
    with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))

def validate():
    p=load('analysis/orbit/protocol.json');s=load(p['source']);r=load(OUTPUTS[0])
    cells=[html.unescape(re.sub('<[^>]*>','',v)).strip() for v in re.findall(r'<td[^>]*>(.*?)(?=<td|$)',(ROOT/'data/orbit/ngc3201-row.html').read_text())]
    assert cells==s['row_cells'] and len(cells)==19
    values=[float(x.split('±')[0]) for x in cells[1:]]
    assert values[10:16]==s['baseline_point_kpc_kms']
    assert values[0:2]==[s['RA_deg'],s['DEC_deg']]
    assert values[4]==s['distance_kpc'] and values[6]==s['radial_velocity_kms']
    assert values[7:10]==[s['pmra_cosdec_masyr'],s['pmdec_masyr'],s['pm_correlation']]
    assert values[16:18]==[s['catalog_mean_pericenter_kpc'],s['catalog_mean_apocenter_kpc']]
    raw=ROOT/s['raw_path']
    if raw.exists():assert hashlib.sha256(raw.read_bytes()).hexdigest()==s['raw_sha256']
    assert r['revision']==p['revision']=='galactic-orbit-1' and r['all_checks_passed']
    assert r['galactic_orbits_integrated'] and not r['remnant_trajectories_performed'] and not r['host_re_equilibrated']
    assert r['protocol_sha256']==hashes(['analysis/orbit/protocol.json'])['analysis/orbit/protocol.json']
    assert [c['case'] for c in r['cases']]==p['cases']
    lim=p['acceptance']
    for c in r['cases']:
        assert c['energy_relative_drift']<lim['relative_energy_drift'] and c['Lz_relative_drift']<lim['relative_Lz_drift']
        assert c['pericenter_mean_kpc']<c['apocenter_mean_kpc'] and c['pericenter_count']>=2 and c['apocenter_count']>=2
        assert 0<c['r_sufficient_1pct_pc_min']<c['r_sufficient_1pct_pc_max']
        for key,tol in [('refinement_max_position_kpc','refined_position_kpc'),('refinement_max_velocity_kms','refined_velocity_kms'),('reversal_position_error_kpc','time_reversal_position_kpc')]:
            if c['case'] in ['catalog_reference', 'host_distance']:
                assert key in c and math.isfinite(c[key]) and 0 <= c[key] < lim[tol]
    for label in ['pericenter','apocenter']:assert r['checks']['catalog_'+label+'_fractional_difference']<lim['catalog_peri_apo_fractional_difference']
    assert r['checks']['sufficient_radius_cadence_relative_difference']<lim['bound_radius_relative_refinement']
    refined=r['checks']['refined_minimum_sufficient_radius_pc']
    assert 0 < r['checks']['proposed_fixed_isolated_control_pc'] == math.floor(refined)
    assert math.isclose(abs(refined/r['cases'][1]['r_sufficient_1pct_pc_min']-1), r['checks']['sufficient_radius_cadence_relative_difference'], rel_tol=1e-12)
    tracks=rows(OUTPUTS[1]);events=rows(OUTPUTS[2]);expected=round((p['integration']['time_Myr'][1]-p['integration']['time_Myr'][0])/p['integration']['output_step_Myr'])+1
    assert len(tracks)==expected*len(p['cases'])
    for c in r['cases']:
        case=[row for row in tracks if row['case']==c['case']]
        assert len(case)==expected and float(case[0]['t_Myr'])==0 and float(case[-1]['t_Myr'])==1000
        for name in ['r_sufficient_1pct_pc','frozen_frame_scale_pc','eta_bound_rh','timescale_ratio_rh']:
            vals=[float(row[name]) for row in case]
            assert all(math.isfinite(x) for x in vals)
            assert math.isclose(min(vals),c[name+'_min'],rel_tol=1e-12) and math.isclose(max(vals),c[name+'_max'],rel_tol=1e-12)
        assert all(float(row['eta_bound_rh'])>=float(row['eta_linear_rh']) for row in case)
    assert all(e['event'] in ['pericenter','apocenter','plane_up','plane_down'] for e in events)
    return {'revision':p['revision'],'scope':p['scope'],'input_sha256':hashes(INPUTS),'output_sha256':hashes(OUTPUTS)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');a=parser.parse_args()
    value=validate();out=ROOT/'results/orbit/provenance.json'
    if a.record:out.write_text(json.dumps(value,indent=2)+'\n')
    else:assert load(str(out.relative_to(ROOT)))==value,'Orbit evidence or input changed'
    print('Verified four Galactic reference orbits and host-support diagnostics; no remnant trajectories or host refit.')
if __name__=='__main__':main()
