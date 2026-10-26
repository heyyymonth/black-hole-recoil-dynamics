"""Evaluate the anisotropic collision kernels of a numerically checked host."""
from pathlib import Path
import csv,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis/host_reconstruction'))
from host import Host
from coefficients import rosenbluth

OUT=ROOT/'results/local-transport-kernels'

def run():
    protocol_path=Path(__file__).with_name('protocol.json');p=json.loads(protocol_path.read_text())
    if p['log_lambda'] != 1:
        raise ValueError('This per-unit-lnLambda output requires log_lambda=1')
    if set(p['orientations']) != {'radial', 'tangential'}:
        raise ValueError('Expected declared radial and tangential orientations')
    base=ROOT/'results/host-reconstruction'
    host_summary=json.loads((base/'summary.json').read_text())
    if not host_summary['all_numerical_checks_passed']:raise ValueError('Host checks failed; collision evaluation blocked')
    host=Host(dict(np.load(base/'host.npz',allow_pickle=False)))
    rows=[];species=[];tensors=[]
    for fraction in p['r_over_rh']:
        r=host.rh*fraction;cutoff=float(np.sqrt(2*host.psi(r)))
        df=lambda u:host.number_df(r,u)[host.physical]
        for speed in p['v_over_local_cutoff']:
            for direction in (['radial'] if speed==0 else p['orientations']):
                unit=np.array([0.,0.,1.]) if direction=='radial' else np.array([1.,0.,0.])
                velocity=unit*speed*cutoff
                answers=[rosenbluth(df,velocity,cutoff,host.mj[host.physical],p['remnant_mass_Msun'],p['log_lambda'],host.G,order) for order in p['quadrature_orders']]
                coarse,fine=answers[-2:]
                Ascale=max(np.linalg.norm(fine['A']),1e-14)
                Aerror=float(np.linalg.norm(fine['A']-coarse['A'])/Ascale) if speed else 0.0
                zero_residual=float(np.linalg.norm(fine['A'])) if speed==0 else 0.0
                Derror=float(np.linalg.norm(fine['D']-coarse['D'])/np.linalg.norm(fine['D']))
                trace=float(np.trace(fine['D']));eigen=np.linalg.eigvalsh(fine['D'])
                identity=abs(trace/fine['trace_expected']-1)
                drag=-float(fine['A']@unit)
                # At zero velocity use the radial axis to define D components.
                parallel=float(unit@fine['D']@unit);perpendicular=(trace-parallel)/2
                j=len(rows)
                rows.append({'state':j,'r_over_rh':fraction,'r_pc':r,'orientation':direction,'v_over_cutoff':speed,'v_kms':speed*cutoff,'drag_over_lnLambda_kms_Myr':drag,'D_parallel_over_lnLambda_kms2_Myr':parallel,'D_perpendicular_mean_over_lnLambda_kms2_Myr':perpendicular,'trace_D_over_lnLambda_kms2_Myr':trace,'A_refinement_error':Aerror,'zero_speed_drift_residual_kms_Myr':zero_residual,'D_refinement_relative_error':Derror,'trace_identity_relative_error':identity,'minimum_D_eigenvalue':float(eigen[0])})
                tensors.append({'state':j,'A_over_lnLambda':fine['A'].tolist(),'D_over_lnLambda':fine['D'].tolist()})
                inds=np.where(host.physical)[0]
                for k,b in enumerate(inds):
                    species.append({'state':j,'bin':int(b),'type':str(host.types[b]),'expected_number_in_host':float(host.Mj[b]/host.mj[b]),'drag_over_lnLambda_kms_Myr':-float(fine['A_species'][k]@unit),'trace_D_over_lnLambda_kms2_Myr':float(np.trace(fine['D_species'][k])),'trace_D_fraction':float(np.trace(fine['D_species'][k])/trace)})
                print('state',j,'refinement',Aerror,Derror,flush=True)
    checks={'zero_speed_symmetry_pass':all(r['zero_speed_drift_residual_kms_Myr']<p['acceptance']['zero_speed_drift_absolute_kms_Myr'] for r in rows),'all_refinement_pass':all(max(r['A_refinement_error'],r['D_refinement_relative_error'])<p['acceptance']['relative_vector_and_tensor_refinement'] for r in rows),'trace_identity_pass':all(r['trace_identity_relative_error']<p['acceptance']['trace_identity_relative'] for r in rows),'positive_semidefinite_pass':all(r['minimum_D_eigenvalue']>=-p['acceptance']['negative_eigenvalue_tolerance'] for r in rows),'nonpositive_energy_drift_from_friction':all(r['drag_over_lnLambda_kms_Myr']>=-1e-12 for r in rows)}
    summary={'revision':'local-transport-kernels-1','states':len(rows),'remnant_mass_Msun':p['remnant_mass_Msun'],'normalization':'per unit constant ln Lambda','quadrature_orders':p['quadrature_orders'],'checks':checks,'all_numerical_checks_passed':all(checks.values()),'max_zero_speed_drift_residual_kms_Myr':max(r['zero_speed_drift_residual_kms_Myr'] for r in rows),'max_A_refinement_error':max(r['A_refinement_error'] for r in rows),'max_D_refinement_relative_error':max(r['D_refinement_relative_error'] for r in rows),'max_trace_identity_relative_error':max(r['trace_identity_relative_error'] for r in rows),'host_sha256':hashlib.sha256((base/'host.npz').read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256(protocol_path.read_bytes()).hexdigest(),'interpretation':'Local coefficients conditional on a fixed reconstructed host and weak-scattering closure, not calibrated friction times or remnant trajectories.'}
    OUT.mkdir(parents=True,exist_ok=True)
    for name,contents in [('states.csv',rows),('species.csv',species)]:
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(contents[0]),lineterminator="\n");w.writeheader();w.writerows(contents)
    (OUT/'tensors.json').write_text(json.dumps(tensors,indent=2,allow_nan=False)+'\n')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))
    if not summary['all_numerical_checks_passed']:
        raise AssertionError('Collision validation failed; outputs must not be adopted')

if __name__=='__main__':run()
