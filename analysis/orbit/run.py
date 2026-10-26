"""Execute the registered Galactic-orbit and host-support audit."""
from pathlib import Path
import csv,hashlib,json,sys
import numpy as np
import scipy
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis/host_reconstruction'));sys.path.insert(0,str(ROOT/'analysis/boundary'))
from host import Host
from model import EnclosedHost
from galaxy import Galaxy
from orbits import initial_state,integrate,conserved,host_diagnostics

def write_csv(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)

def main():
    protocol=json.loads((ROOT/'analysis/orbit/protocol.json').read_text());source=json.loads((ROOT/protocol['source']).read_text())
    host=EnclosedHost(Host(dict(np.load(ROOT/protocol['host'],allow_pickle=False))))
    hs=json.loads((ROOT/protocol['host_summary']).read_text());galaxy=Galaxy(protocol['galaxy']);cfg=protocol['integration']
    distances=[source['distance_kpc'],hs['summary']['distance_kpc'],source['distance_kpc']-source['distance_error_kpc'],source['distance_kpc']+source['distance_error_kpc']]
    interval=cfg['time_Myr'];times=np.arange(interval[0],interval[1]+cfg['output_step_Myr']/2,cfg['output_step_Myr'])
    cases=[];events=[];profiles=[];checks={};out=ROOT/'results/orbit';out.mkdir(exist_ok=True)
    for name,distance in zip(protocol['cases'],distances):
        print('Integrating',name,flush=True)
        y0=initial_state(source,distance)
        sol=integrate(galaxy,y0,interval,cfg['max_step_Myr'],cfg['rtol'],cfg['atol'])
        states=sol.sol(times).T;energy,lz=conserved(galaxy,states)
        errE=float(np.max(np.abs((energy-energy[0])/energy[0])));errL=float(np.max(np.abs((lz-lz[0])/lz[0])))
        tidal=host_diagnostics(galaxy,host,states,protocol['force_accuracy']['epsilon'],protocol['tensor_derivative_step_Myr'])
        r=np.linalg.norm(states[:,:3],axis=1)
        case={'case':name,'distance_kpc':distance,'initial_state_kpc_kms':y0.tolist(),'energy_relative_drift':errE,'Lz_relative_drift':errL}
        extrema={'pericenter':[],'apocenter':[]}
        for t,y in zip(sol.t_events[0],sol.y_events[0]):
            acc=galaxy.fields(y[:3])[1];kind='pericenter' if np.sum(y[3:]**2)+y[:3]@acc>0 else 'apocenter'
            radius=float(np.linalg.norm(y[:3]));extrema[kind].append(radius)
            events.append({'case':name,'event':kind,'t_Myr':float(t),'R_kpc':radius,'Z_kpc':float(y[2])})
        for t,y in zip(sol.t_events[1],sol.y_events[1]):
            events.append({'case':name,'event':'plane_up' if y[5]>0 else 'plane_down','t_Myr':float(t),'R_kpc':float(np.linalg.norm(y[:3])),'Z_kpc':float(y[2])})
        for label,values in extrema.items():case[label+'_mean_kpc']=float(np.mean(values));case[label+'_count']=len(values)
        for key,values in tidal.items():
            case[key+'_min']=float(np.min(values));case[key+'_max']=float(np.max(values))
        case['maximum_abs_Z_kpc']=float(np.max(np.abs(states[:,2])))
        case['spherical_eccentricity_from_mean_extrema']=(case['apocenter_mean_kpc']-case['pericenter_mean_kpc'])/(case['apocenter_mean_kpc']+case['pericenter_mean_kpc'])
        if name in ['catalog_reference','host_distance']:
            # Independent step/tolerance refinement and backward integration.
            refined=integrate(galaxy,y0,interval,cfg['refinement_max_step_Myr'],cfg['refinement_rtol'],1e-13)
            refstates=refined.sol(times).T
            case['refinement_max_position_kpc']=float(np.max(np.linalg.norm(refstates[:,:3]-states[:,:3],axis=1)))
            case['refinement_max_velocity_kms']=float(np.max(np.linalg.norm(refstates[:,3:]-states[:,3:],axis=1)))
            reverse=integrate(galaxy,states[-1],interval[::-1],cfg['max_step_Myr'],cfg['rtol'],cfg['atol'])
            case['reversal_position_error_kpc']=float(np.linalg.norm(reverse.y[:3,-1]-y0[:3]))
        for i,t in enumerate(times):
            profiles.append({'case':name,'t_Myr':float(t),**dict(zip(['X_kpc','Y_kpc','Z_kpc','U_kms','V_kms','W_kms'],states[i])),**{k:float(v[i]) for k,v in tidal.items()}})
        cases.append(case)
    ref=cases[0]
    for label in ['pericenter','apocenter']:
        checks['catalog_'+label+'_fractional_difference']=abs(ref[label+'_mean_kpc']/source['catalog_mean_'+label+'_kpc']-1)
    # Cadence convergence of the most restrictive force domain and response indicators.
    sol=integrate(galaxy,initial_state(source,distances[1]),interval,cfg['refinement_max_step_Myr'],cfg['refinement_rtol'],1e-13)
    fine_times=np.arange(interval[0],interval[1]+cfg['output_step_Myr']/4,cfg['output_step_Myr']/2)
    fine=host_diagnostics(galaxy,host,sol.sol(fine_times).T,protocol['force_accuracy']['epsilon'],protocol['tensor_derivative_step_Myr']/2)
    checks['sufficient_radius_cadence_relative_difference']=abs(float(np.min(fine['r_sufficient_1pct_pc']))/cases[1]['r_sufficient_1pct_pc_min']-1)
    checks['refined_minimum_sufficient_radius_pc']=float(np.min(fine['r_sufficient_1pct_pc']))
    checks['proposed_fixed_isolated_control_pc']=float(np.floor(checks['refined_minimum_sufficient_radius_pc']))
    checks['timescale_ratio_rh_refinement_absolute_difference']=abs(float(np.max(fine['timescale_ratio_rh']))-cases[1]['timescale_ratio_rh_max'])
    checks['solar_circle_vc_at_8_4_kpc_kms']=float(np.sqrt(-8.4*galaxy.fields([8.4,0,0])[1][0]))
    limits=protocol['acceptance']
    passed=all(c['energy_relative_drift']<limits['relative_energy_drift'] and c['Lz_relative_drift']<limits['relative_Lz_drift'] for c in cases)
    passed=passed and all(c.get('refinement_max_position_kpc',0)<limits['refined_position_kpc'] and c.get('refinement_max_velocity_kms',0)<limits['refined_velocity_kms'] and c.get('reversal_position_error_kpc',0)<limits['time_reversal_position_kpc'] for c in cases)
    passed=passed and max(checks['catalog_pericenter_fractional_difference'],checks['catalog_apocenter_fractional_difference'])<limits['catalog_peri_apo_fractional_difference'] and checks['sufficient_radius_cadence_relative_difference']<limits['bound_radius_relative_refinement']
    result={'revision':protocol['revision'],'galactic_orbits_integrated':True,'remnant_trajectories_performed':False,'host_re_equilibrated':False,'cases':cases,'checks':checks,'all_checks_passed':bool(passed),'protocol_sha256':hashlib.sha256((ROOT/'analysis/orbit/protocol.json').read_bytes()).hexdigest(),'environment':{'python':sys.version.split()[0],'numpy':np.__version__,'scipy':scipy.__version__}}
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');write_csv(out/'tracks.csv',profiles);write_csv(out/'events.csv',events)
    print(json.dumps(result,indent=2),flush=True)
    if not passed:raise AssertionError('Orbit checks failed; result is not approved')
if __name__=='__main__':main()
