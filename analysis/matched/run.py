"""Execute registered independent ensembles; do not change controls from outcomes."""
from pathlib import Path
import argparse
import itertools
import json
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analysis.matched.kernel import MatchedKernel,digest
from analysis.matched.transport import run_ensemble

PROTOCOL=ROOT/'analysis/matched/protocol.json'
OUT=ROOT/'results/matched-residence'
CODE=['analysis/matched/run.py','analysis/matched/transport.py','analysis/matched/kernel.py',
      'analysis/host_reconstruction/host.py','analysis/encounters/physics.py','analysis/encounters/run.py']
ARRAYS=('states','residence','normalized_residence','exposure','normalized_exposure','end_time',
        'stop_cause','alive','elastic_events','evaluated_steps')


def settings(p):
    return list(itertools.product(p['transport']['grids'],p['transport']['dt_over_t_dyn'],
                                  p['launch']['orientations'],p['transport']['methods']))


def setting_name(grid,dt,orientation,method):
    return f'{grid}-dt{dt:.3f}-{orientation}-{method}'


def host_functions(kernel):
    h=kernel.host;rscale=kernel.r_scale_pc;vscale=kernel.v_scale_kms
    bh=h.types=='BH';rho0=float(h.rho_spline(0)[bh].sum())
    def acceleration(x):
        r=np.linalg.norm(x,axis=1)
        derivative=h.psi_spline(r*rscale,1)
        scale=np.divide(rscale*derivative,vscale*vscale*r,out=np.zeros_like(r),where=r>0)
        return x*scale[:,None]
    def density(x):
        return h.rho_spline(np.linalg.norm(x,axis=1)*rscale)[bh].sum(axis=0)/rho0
    return acceleration,density


def initial_state(kernel,p,orientation):
    r0=p['launch']['radius_over_density_radius']
    ra=p['launch']['radial_apocenter_over_density_radius']
    speed=np.sqrt(2*(kernel.host.psi(r0*kernel.r_scale_pc)
                     -kernel.host.psi(ra*kernel.r_scale_pc)))/kernel.v_scale_kms
    axis=np.array([1.,0.,0.]) if orientation=='radial_outward' else np.array([0.,1.,0.])
    return np.r_[r0,0.,0.,speed*axis]


def identity(grid,p):
    return {'protocol_sha256':digest(PROTOCOL),
            'host_sha256':digest(ROOT/p['host_path']),
            'encounter_protocol_sha256':digest(ROOT/p['encounter_protocol']),
            'kernel_sha256':digest(OUT/f'kernel-{grid}.npz'),
            'kernel_record_sha256':digest(OUT/f'kernel-{grid}.json'),
            'implementation_sha256':{name:digest(ROOT/name) for name in CODE}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',choices=['coarse','fine'],required=True)
    parser.add_argument('--method',choices=['diffusion','jump'])
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    p=json.loads(PROTOCOL.read_text())
    metadata=json.loads((OUT/f'kernel-{args.grid}.json').read_text())
    if not (metadata['completed'] and metadata['all_vertex_quadrature_passed']
            and metadata['all_holdout_quadrature_passed']):
        raise ValueError('Completed local quadrature and holdout validation required')
    k=MatchedKernel(ROOT/p['host_path'],OUT/f'kernel-{args.grid}.npz',PROTOCOL)
    force,density=host_functions(k)
    # Inherited isolated-force limit is orders of magnitude beyond this table.
    baseline=json.loads((ROOT/'results/collisionless-residence/NGC5139.json').read_text())
    validity=min(r['summary']['stop_radius_pc'] for r in baseline['records'])/k.r_scale_pc
    spatial=min(k.axes[0][-1],validity);velocity=k.axes[1][-1]
    n=p['transport']['paths_per_setting'];horizon=p['horizon_over_t_dyn']
    for index,(grid,dt,orientation,method) in enumerate(settings(p)):
        if grid!=args.grid or (args.method is not None and args.method!=method):continue
        name=setting_name(grid,dt,orientation,method);path=OUT/f'{name}.npz'
        ids=identity(grid,p)
        if args.resume and path.exists():
            previous=json.loads(path.with_suffix('.json').read_text())
            if all(previous.get(key)==value for key,value in ids.items()) and previous['samples_sha256']==digest(path):
                print('verified resume',name,flush=True);continue
            raise ValueError('Stale output cannot be resumed')
        seed=[p['transport']['seed'],index]
        initial=initial_state(k,p,orientation)
        print('starting',name,'paths',n,flush=True);start=time.monotonic()
        result=run_ensemble(k,force,np.repeat(initial[None,:],n,axis=0),horizon,dt,seed,method,
                             spatial_stop=spatial,velocity_stop=velocity,density=density,
                             region_radius=p['region_radius_over_density_radius'])
        np.savez_compressed(path,**{key:result[key] for key in ARRAYS})
        record={**ids,'revision':p['revision'],'cluster':p['cluster'],'grid':grid,'dt':dt,
                'orientation':orientation,'method':method,'paths':n,'horizon':horizon,
                'seed':seed,'seed_scheme':'NumPy SeedSequence([registered_seed, fixed_product_index])',
                'numpy_version':np.__version__,'initial_state':initial.tolist(),
                'scales':metadata['scales'],'spatial_stop':float(spatial),'velocity_stop':float(velocity),
                'inherited_force_domain_over_density_radius':float(validity),
                'units':'r_rho, v_dyn, t_dyn; normalized residence/exposure divided by full horizon',
                'samples_sha256':digest(path),'scheme':result['scheme'],
                'wall_seconds':time.monotonic()-start}
        path.with_suffix('.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        print('finished',name,'seconds',round(record['wall_seconds'],2),
              'mean residence',float(result['normalized_residence'].mean()),
              'stops',dict(zip(*np.unique(result['stop_cause'],return_counts=True))),flush=True)


if __name__=='__main__':main()
