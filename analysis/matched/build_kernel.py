"""Build and check positive finite-angle rate-measure interpolation tables."""
from pathlib import Path
import argparse
import json
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from analysis.matched.kernel import (MatchedKernel,Host,canonical_velocity,
    center_rotation,digest,host_scales,refine_vertex)

PROTOCOL=ROOT/'analysis/matched/protocol.json'
OUT=ROOT/'results/matched-residence'
IMPLEMENTATION=['analysis/matched/kernel.py','analysis/matched/build_kernel.py',
                'analysis/host_reconstruction/host.py','analysis/encounters/physics.py',
                'analysis/encounters/run.py']


def axes_for(protocol,grid):
    axes=[np.asarray(protocol['kernel']['coarse_axes'][key],dtype=float)
          for key in ('radius','speed','pitch')]
    if grid=='fine':
        axes=[np.sort(np.r_[axis,(axis[1:]+axis[:-1])/2]) for axis in axes]
    elif grid!='coarse':
        raise ValueError('Only the two registered grids are supported')
    return axes


def identity(protocol):
    return {'protocol_sha256':digest(PROTOCOL),'host_sha256':digest(ROOT/protocol['host_path']),
            'encounter_protocol_sha256':digest(ROOT/protocol['encounter_protocol']),
            'implementation_sha256':{name:digest(ROOT/name) for name in IMPLEMENTATION}}


def state_at(coordinates):
    r,s,pitch=coordinates
    return np.r_[0.,0.,r,canonical_velocity(s,pitch)]


def finite_summary(answer):
    return {key:np.asarray(answer[key]).tolist() for key in ('rate','excluded','A','D')}


def save_checkpoint(path,data,metadata):
    temporary=path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**data)
    temporary.replace(path)
    metadata['table_sha256']=digest(path)
    temp_json=path.with_suffix('.tmp.json')
    temp_json.write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    temp_json.replace(path.with_suffix('.json'))


def build(grid,resume=False):
    protocol=json.loads(PROTOCOL.read_text())
    encounter=json.loads((ROOT/protocol['encounter_protocol']).read_text())
    host_path=ROOT/protocol['host_path']
    host_record=json.loads(host_path.with_name('summary.json').read_text())
    if not host_record['all_numerical_checks_passed']:
        raise ValueError('The adopted host must first pass its numerical checks')
    host=Host(dict(np.load(host_path,allow_pickle=False)))
    scales=host_scales(host)
    axes=axes_for(protocol,grid)
    shape=tuple(map(len,axes));count=int(np.prod(shape))
    factors={'rate':scales[2],'excluded':scales[2],'A':scales[2]/scales[1],
             'D':scales[2]/scales[1]**2}
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/f'kernel-{grid}.npz'
    keys=identity(protocol)
    data={'radius':axes[0],'speed':axes[1],'pitch':axes[2],
          'rate':np.zeros(shape),'excluded':np.zeros(shape),
          'A':np.zeros(shape+(3,)),'D':np.zeros(shape+(3,3)),
          'order':np.zeros(shape,dtype=int),'completed':np.zeros(shape,dtype=bool),
          'accepted':np.zeros(shape,dtype=bool),'scales':np.asarray(scales)}
    records=[None]*count
    if resume and path.exists():
        previous=json.loads(path.with_suffix('.json').read_text())
        if any(previous.get(key)!=value for key,value in keys.items()):
            raise ValueError('Checkpoint does not match current source/protocol/host')
        if previous['table_sha256']!=digest(path):
            raise ValueError('Checkpoint table hash mismatch')
        data=dict(np.load(path,allow_pickle=False));records=previous['vertices']
        if any(not np.array_equal(data[key],axis) for key,axis in zip(('radius','speed','pitch'),axes)):
            raise ValueError('Checkpoint grid differs')
    metadata={**keys,'revision':protocol['revision'],'grid':grid,'shape':list(shape),
              'vertex_count':count,'units':'r_rho, v_dyn, t_dyn; rates and moments per dimensionless time',
              'scales':{'r_scale_pc':scales[0],'v_scale_kms':scales[1],'time_scale_Myr':scales[2]},
              'interpolation':protocol['kernel']['interpolation'],
              'symmetry':'Rate measures and moments averaged over the same physical symmetries; center velocity orientation is rotated from radial canonical law',
              'vertices':records,'all_vertex_quadrature_passed':False,'completed':False}
    # Reuse the rotationally identical central and zero-speed calculations.
    base_cache={}
    for flat in range(count):
        index=np.unravel_index(flat,shape)
        if data['completed'][index]:
            continue
        coordinates=tuple(float(axes[i][j]) for i,j in enumerate(index))
        r,s,pitch=coordinates
        base=(r,s,0. if r==0 or s==0 else pitch)
        if base not in base_cache:
            base_cache[base]=refine_vertex(host,base,scales,protocol,encounter)
        base_answer,base_record=base_cache[base]
        answer=dict(base_answer)
        record={**base_record,'coordinates':list(coordinates)}
        if r==0 and s>0:
            rotation=center_rotation(pitch)
            answer['A']=rotation@base_answer['A']
            answer['D']=rotation@base_answer['D']@rotation.T
        for key,factor in factors.items():
            data[key][index]=np.asarray(answer[key])*factor
        data['order'][index]=record['order']
        data['accepted'][index]=record['accepted']
        data['completed'][index]=True
        records[flat]=record
        if (flat+1)%25==0 or flat+1==count:
            metadata['accepted_vertices']=int(data['accepted'].sum())
            metadata['completed_vertices']=int(data['completed'].sum())
            save_checkpoint(path,data,metadata)
            print(grid,metadata['completed_vertices'],'/',count,
                  'accepted',metadata['accepted_vertices'],flush=True)
    metadata['all_vertex_quadrature_passed']=bool(data['accepted'].all())
    metadata['completed_vertices']=count
    metadata['accepted_vertices']=int(data['accepted'].sum())
    save_checkpoint(path,data,metadata)
    if not metadata['all_vertex_quadrature_passed']:
        print('UNRESOLVED vertex quadrature; no interpolation validation or trajectories',flush=True)
        return False
    # Use the same positive interpolation as a real trajectory, without
    # constructing its event CDF until an actual event is sampled.
    kernel=MatchedKernel(host_path,path,PROTOCOL)
    rng=np.random.default_rng(protocol['kernel']['holdout_seed'])
    points=rng.uniform([0,0,0],[axes[0][-1],axes[1][-1],1.],
                       size=(protocol['kernel']['holdout_points'],3))
    # Half sample the region traversed by the registered reference launches.
    points[:len(points)//2,:2]*=.5
    controls=np.array([[0,0,0],[0,1,0],[.5,0,0],[.5,1,0],
                       [.5,1,1e-4],[.5,1,1-1e-4],[.5,1,1],[1,1,.5]])
    points=np.vstack([points,controls])
    predictions=kernel.moments(np.array([state_at(point) for point in points]))
    holdouts=[]
    for i,point in enumerate(points):
        answer,record=refine_vertex(host,point,scales,protocol,encounter)
        direct={key:np.asarray(answer[key])*factor for key,factor in factors.items()}
        prediction={key:predictions[key][i] for key in factors}
        sigma=answer['sigma']/scales[1]
        scale=max(np.linalg.norm(direct['A']),np.trace(direct['D'])/sigma)
        errors={'rate_relative_error':float(abs(prediction['rate']/direct['rate']-1)),
                'excluded_relative_error':float(abs(prediction['excluded']/direct['excluded']-1)),
                'scaled_drift_error':float(np.linalg.norm(prediction['A']-direct['A'])/scale),
                'diffusion_relative_error':float(np.linalg.norm(prediction['D']-direct['D'])/np.linalg.norm(direct['D']))}
        holdouts.append({'coordinates':point.tolist(),'direct':finite_summary(direct),
                         'interpolated':finite_summary(prediction),'errors':errors,
                         'direct_quadrature':record})
    metadata['holdout_points']=holdouts
    metadata['all_holdout_quadrature_passed']=all(row['direct_quadrature']['accepted'] for row in holdouts)
    metadata['holdout_error_maxima']={key:max(row['errors'][key] for row in holdouts)
                                      for key in holdouts[0]['errors']}
    metadata['holdout_interpretation']=protocol['kernel']['validation_interpretation']
    metadata['completed']=True
    save_checkpoint(path,data,metadata)
    print(grid,'complete',json.dumps(metadata['holdout_error_maxima']),flush=True)
    return metadata['all_holdout_quadrature_passed']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',choices=['coarse','fine'],required=True)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not build(args.grid,args.resume):
        raise SystemExit(1)
