"""Persist host-force and collisionless controls for matched transport.

Run this module to evaluate controls; --verify only needs NumPy and the standard
library. Criteria below are declared before execution. A null trajectory check
validates deterministic transport, not the collisional residence prediction.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/matched-residence/controls.json'
CRITERIA={
    'force_relative_error':.001,
    'mass_quadrature_relative_change':1e-5,
    'normalized_residence_absolute_error':1e-4,
    'normalized_exposure_absolute_error':1e-4,
    'fine_error_no_larger_than_coarse':True,
}
SOURCES=['analysis/matched/controls.py','analysis/matched/run.py',
         'analysis/matched/transport.py','analysis/matched/kernel.py',
         'analysis/host_reconstruction/host.py','analysis/encounters/physics.py',
         'analysis/encounters/run.py']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checks(record):
    force=record['force']
    numerical=np.array(force['numerical_vectors'])
    reference=np.array(force['density_integral_vectors'])
    relative=np.linalg.norm(numerical-reference,axis=1)/np.linalg.norm(reference,axis=1)
    masses=np.array(force['enclosed_mass_by_order_Msun'])
    change=np.abs(masses[0]/masses[1]-1)
    passed={'origin_force_exact_zero':bool(np.all(np.array(force['origin_force'])==0)),
            'force_relative_error':bool(np.max(relative)<=CRITERIA['force_relative_error']),
            'mass_quadrature_relative_change':bool(np.max(change)<=CRITERIA['mass_quadrature_relative_change'])}
    grouped={}
    for row in record['null_trajectories']:
        name=f"{row['orientation']}-dt{row['dt']}"
        if row['stop_cause']!='horizon' or row['end_time']!=record['horizon']:
            passed[name+'-horizon']=False
        else:
            passed[name+'-horizon']=True
        errors={}
        for observable in ('normalized_residence','normalized_exposure'):
            errors[observable]=abs(row[observable]-row['reference_'+observable])
            passed[name+'-'+observable]=errors[observable]<=CRITERIA[observable+'_absolute_error']
        grouped.setdefault(row['orientation'],[]).append((row['dt'],errors))
    for orientation,rows in grouped.items():
        rows=sorted(rows)
        if len(rows)!=2:
            raise ValueError('Expected exactly two timestep controls per orientation')
        for observable in ('normalized_residence','normalized_exposure'):
            passed[orientation+'-'+observable+'-refinement']=rows[0][1][observable]<=rows[1][1][observable]
    return {'max_force_relative_error':float(np.max(relative)),
            'max_mass_quadrature_relative_change':float(np.max(change)),
            'checks':{key:bool(value) for key,value in passed.items()},
            'accepted':bool(all(passed.values()))}


def generate():
    sys.path.insert(0,str(ROOT))
    from numpy.polynomial.legendre import leggauss
    from analysis.matched.kernel import MatchedKernel
    from analysis.matched.run import host_functions,initial_state
    from analysis.matched.transport import run_ensemble

    protocol_path=ROOT/'analysis/matched/protocol.json'
    p=json.loads(protocol_path.read_text())
    kernel_path=ROOT/'results/matched-residence/kernel-coarse.npz'
    kernel=MatchedKernel(ROOT/p['host_path'],kernel_path,protocol_path)
    acceleration,density=host_functions(kernel)
    baseline_path=ROOT/'results/collisionless-residence/NGC5139.json'
    baseline=json.loads(baseline_path.read_text())
    references={row['summary']['orientation']:row['summary'] for row in baseline['records']}
    # Direct all-species mass integration is independent of psi_spline's
    # derivative and does not interpolate stored cumulative-mass arrays.
    radii=np.array([.01,.02,.05,.1,.25,.5,.75,1.,1.5,2.,3.,4.])
    positions=radii[:,None]*np.array([1.,2.,2.])[None,:]/3
    masses=[]
    for order in (64,128):
        nodes,weights=leggauss(order)
        radius=kernel.r_scale_pc*radii
        radial=radius[:,None]*(nodes+1)/2
        rho=kernel.host.rho_spline(radial.ravel()).sum(axis=0).reshape(radial.shape)
        masses.append(2*np.pi*radius*np.sum(weights*rho*radial**2,axis=1))
    force_size=kernel.host.G*masses[-1]/(kernel.r_scale_pc*radii**2*kernel.v_scale_kms**2)
    reference_force=-force_size[:,None]*positions/radii[:,None]
    force_record={'radii_over_density_radius':radii.tolist(),
                  'positions':positions.tolist(),'quadrature_orders':[64,128],
                  'enclosed_mass_by_order_Msun':np.array(masses).tolist(),
                  'numerical_vectors':acceleration(positions).tolist(),
                  'density_integral_vectors':reference_force.tolist(),
                  'origin_force':acceleration(np.zeros((1,3)))[0].tolist(),
                  'scope':'Full host density integrated independently; dimensionless acceleration'}

    class NullKernel:
        def moments(self,states):
            n=len(states)
            return {'rate':np.zeros(n),'excluded':np.zeros(n),
                    'A':np.zeros((n,3)),'D':np.zeros((n,3,3))}
        def sample_events(self,states,rng):
            return np.zeros((len(states),3))

    null_rows=[]
    for orientation in p['launch']['orientations']:
        ref=references[orientation]
        for dt in p['transport']['dt_over_t_dyn']:
            result=run_ensemble(NullKernel(),acceleration,initial_state(kernel,p,orientation)[None,:],
                                p['horizon_over_t_dyn'],dt,2026110103,'diffusion',
                                spatial_stop=4.,velocity_stop=4.,density=density,
                                region_radius=p['region_radius_over_density_radius'])
            null_rows.append({'orientation':orientation,'dt':dt,
                              'normalized_residence':float(result['normalized_residence'][0]),
                              'reference_normalized_residence':ref['residence_fraction'],
                              'normalized_exposure':float(result['normalized_exposure'][0]),
                              'reference_normalized_exposure':ref['normalized_BH_density_exposure'],
                              'end_time':float(result['end_time'][0]),
                              'stop_cause':str(result['stop_cause'][0]),
                              'terminal_state':result['states'][0].tolist()})
    inputs=['analysis/matched/protocol.json',p['host_path'],p['encounter_protocol'],
            p['baseline_protocol'],'analysis/host_reconstruction/protocol.json',
            str(kernel_path.relative_to(ROOT)),str(kernel_path.with_suffix('.json').relative_to(ROOT)),
            str(baseline_path.relative_to(ROOT))]
    record={'revision':'matched-transport-controls-1','cluster':p['cluster'],
            'criteria':CRITERIA,
            'criteria_basis':'Force tolerance inherited from host validation; mass quadrature 1e-5; null normalized observable tolerance 1e-4, 2 percent of the registered 0.005 mean numerical budget; each fine-step error must not exceed its coarse error. Declared in source before evaluation.',
            'horizon':p['horizon_over_t_dyn'],'numpy_version':np.__version__,
            'force':force_record,'null_trajectories':null_rows,
            'source_sha256':{name:digest(ROOT/name) for name in SOURCES},
            'input_sha256':{name:digest(ROOT/name) for name in inputs},
            'interpretation':'Deterministic and force controls only; not collisional timestep or kernel-grid convergence.'}
    record.update(checks(record))
    OUT.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    if not record['accepted']:
        raise ValueError('Control criteria failed; all results retained in controls.json')
    return record


def verify():
    record=json.loads(OUT.read_text())
    if record['criteria']!=CRITERIA:
        raise ValueError('Criteria mismatch')
    for group in ('source_sha256','input_sha256'):
        for name,expected in record[group].items():
            if digest(ROOT/name)!=expected:
                raise ValueError('Stale control dependency: '+name)
    recalculated=checks(record)
    for key,value in recalculated.items():
        equal = (np.isclose(record[key],value,rtol=2e-11,atol=2e-13)
                 if isinstance(value,float) else record[key]==value)
        if not equal:
            raise ValueError('Stored control result differs from recomputation: '+key)
    if not recalculated['accepted']:
        raise ValueError('Control criteria not accepted')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=verify() if args.verify else generate()
    print(json.dumps({key:result[key] for key in ('accepted','max_force_relative_error',
                                                 'max_mass_quadrature_relative_change')},indent=2))
