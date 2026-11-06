"""Execute fixed independent symmetric-root predictions before fresh validation."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse
import gc
import json
import os
import shutil
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analysis.response import run as reference
from analysis.response.disk_cache import DiskMatchedKernel
from analysis.precision.coupling import make_transport

OUT=ROOT/'results/retest'
PROTOCOL=ROOT/'analysis/retest/protocol.json'
read=reference.read_json
digest=reference.digest
CODE=('analysis/retest/run.py','analysis/retest/analyze.py','analysis/retest/statistics.py',
      'analysis/precision/coupling.py','analysis/response/disk_cache.py')


def jobs(p,stage):
    index=0
    for si,setting in enumerate(p['numerical_settings']):
        for li,launch in enumerate(p['launches']):
            if stage=='prediction':
                variants=[None];method='weighted'
                n=p['prediction']['primary_quartets'] if setting['role']=='primary' else p['prediction']['control_quartets']
            else:
                variants=p['validation']['primary_a_values'] if setting['role']=='primary' else p['validation']['control_a_values']
                method='paired_mixture';n=p['validation']['pairs_per_setting']
            for vi,a in enumerate(variants):
                name=reference.job_name(stage,setting,launch,method,a)
                yield {'name':name,'stage':stage,'setting_index':si,'launch_index':li,
                    'variant_index':vi,'setting':setting,'launch':launch,'method':method,'a':a,
                    'paths':n,'seed':[p['seed'],50 if stage=='prediction' else 60,si,li,vi],
                    'worker':index%p['workers']}
                index+=1


def identity(p):
    rp=read(ROOT/p['reference_protocol']);rid=reference.identity(rp)
    if rid!=read(ROOT/p['reference_design'])['identity']:
        raise ValueError('Frozen reference inputs changed')
    precision=read(ROOT/p['precision_design'])
    for path,sha in precision['identity']['source_sha256'].items():
        if digest(ROOT/path)!=sha:raise ValueError('Frozen precision source changed: '+path)
    names=[str(PROTOCOL.relative_to(ROOT)),p['reference_design'],p['precision_design'],
           p['development_evidence'],'results/response/summary.json','results/retest/preflight.txt',*CODE]
    return {'reference_identity':rid,'source_sha256':{name:digest(ROOT/name) for name in names}}


def ensure_design(p,ids):
    core={'revision':p['revision'],'scope':p['scope'],'identity':ids,
          'jobs':list(jobs(p,'prediction'))+list(jobs(p,'validation')),
          'stage_order':['prediction','prediction-register','validation']}
    path=OUT/'execution-design.json'
    if path.exists():
        old=read(path)
        if any(old.get(k)!=v for k,v in core.items()):raise ValueError('Changed fixed design')
    else:
        if list(OUT.glob('prediction-*.npz')) or list(OUT.glob('validation-*.npz')):
            raise ValueError('Outcomes predate design')
        reference.write_new_json(path,{**core,'created_utc':reference.utc_now()})
    return digest(path)


def checked(job,ids,dependencies):
    path=OUT/(job['name']+'.npz');r=read(path.with_suffix('.json'))
    if not (r['complete'] and r['identity']==ids and r['job']==job
            and r['stage_dependencies']==dependencies and r['samples_sha256']==digest(path)):
        raise ValueError('Changed/incomplete outcome '+job['name'])
    return r


def lock_predictions(p,ids,design):
    deps={'execution_design_sha256':design};inventory=[]
    for job in jobs(p,'prediction'):
        r=checked(job,ids,deps);path=OUT/(job['name']+'.npz')
        inventory.append({'name':job['name'],'arrays':str(path.relative_to(ROOT)),
            'arrays_sha256':digest(path),'record':str(path.with_suffix('.json').relative_to(ROOT)),
            'record_sha256':digest(path.with_suffix('.json')),'summary':r['summary']})
    core={'revision':p['revision'],'identity':ids,**deps,'prediction_inventory':inventory,
          'scope':'All predictions frozen before independent validation; earlier outcomes excluded'}
    path=OUT/'prediction-register.json'
    if path.exists():
        old=read(path)
        if any(old.get(k)!=v for k,v in core.items()):raise ValueError('Changed prediction register')
    else:
        if list(OUT.glob('validation-*.npz')):raise ValueError('Validation predates register')
        reference.write_new_json(path,{**core,'created_utc':reference.utc_now()})
    return digest(path)


def prepare_caches(p):
    for worker in range(p['workers']):
        for grid in ('coarse','fine'):
            source=ROOT/p['cache']['source_root']/grid
            target=ROOT/p['cache']['private_root']/f'worker{worker}'/grid
            if target.exists():continue
            target.parent.mkdir(parents=True,exist_ok=True)
            if source.exists():shutil.copytree(source,target,copy_function=os.link)
            else:target.mkdir()


def worker_main(p,ids,dependencies,stage,worker,resume):
    rp=read(ROOT/p['reference_protocol']);matched=read(ROOT/rp['matched_protocol'])
    transport=make_transport('symmetric');kernel=None;grid=None
    for job in jobs(p,stage):
        if job['worker']!=worker:continue
        path=OUT/(job['name']+'.npz')
        if path.exists() or path.with_suffix('.json').exists():
            if not resume:raise FileExistsError('Existing result '+job['name'])
            checked(job,ids,dependencies);print('verified resume',job['name'],flush=True);continue
        if identity(p)!=ids:raise ValueError('Sources changed before execution')
        setting=job['setting']
        if setting['grid']!=grid:
            kernel=None;gc.collect();grid=setting['grid'];table=ROOT/rp['kernel_tables'][grid]
            meta=read(table.with_suffix('.json'))
            if not all(meta[k] for k in ('completed','all_vertex_quadrature_passed','all_holdout_quadrature_passed')):
                raise ValueError('Unaccepted kernel table')
            cache=p['cache']
            kernel=DiskMatchedKernel(ROOT/rp['host_path'],table,ROOT/rp['matched_protocol'],
                disk_cache_dir=ROOT/cache['private_root']/f'worker{worker}'/grid,
                namespace=digest(ROOT/cache['namespace_record']),grid=grid,
                cache_bytes=cache['retained_array_bytes'],cache_vertices=cache['vertex_cap'],
                disk_max_bytes=cache['disk_max_bytes_per_grid'])
        force,density=reference.host_functions(kernel)
        inherited=read(ROOT/'results/collisionless-residence/NGC5139.json')
        domain=min(r['summary']['stop_radius_pc'] for r in inherited['records'])/kernel.r_scale_pc
        controls={'spatial_stop':float(min(kernel.axes[0][-1],domain)),
            'velocity_stop':float(kernel.axes[1][-1]),'density':density,
            'region_radius':matched['region_radius_over_density_radius']}
        initial=reference.initial_state(kernel,matched,job['launch'])
        args=(kernel,force,np.repeat(initial[None,:],job['paths'],axis=0),p['horizon'],setting['dt'],job['seed'])
        before=kernel.disk_statistics();started=reference.utc_now();wall=time.perf_counter();cpu=time.process_time()
        print('starting',job['name'],'worker',worker,'paths',job['paths'],flush=True)
        if stage=='prediction':result=transport.run_weighted_response(*args,**controls)
        else:result=transport.run_paired_mixture(*args,mixture_probability=job['a'],**controls)
        wall=time.perf_counter()-wall;cpu=time.process_time()-cpu
        arrays,branches=reference.flatten_result(result,job['method'],p['horizon'],setting['dt'])
        if identity(p)!=ids:raise ValueError('Sources changed during execution')
        score='response_sample' if stage=='prediction' else 'paired_difference'
        cache=kernel.disk_statistics();cache['cache_path']=str(Path(cache['cache_path']).relative_to(ROOT))
        r={'complete':True,'revision':p['revision'],'identity':ids,'job':job,'stage_dependencies':dependencies,
            'started_utc':started,'completed_utc':reference.utc_now(),'wall_seconds':wall,'cpu_seconds':cpu,
            'horizon':p['horizon'],'dt':setting['dt'],'n_steps':result.get('n_steps',int(np.ceil(p['horizon']/setting['dt']))),
            'initial_state':initial.tolist(),'spatial_stop':controls['spatial_stop'],'velocity_stop':controls['velocity_stop'],
            'inherited_force_domain':float(domain),'region_radius':controls['region_radius'],
            'scales':{'density_radius_pc':kernel.r_scale_pc,'velocity_kms':kernel.v_scale_kms,'time_Myr':kernel.time_scale_Myr},
            'coupling':'Symmetric PSD root; '+result['coupling'],'branches':list(branches),
            'summary':reference.sample_summary(arrays[score]),
            'numerical_missing_bound':reference.sample_summary(arrays['numerical_missing_bound']),
            'stops':{b:dict(zip(*[x.tolist() for x in np.unique(result[b]['stop_cause'],return_counts=True)])) for b in branches},
            'cache':cache,'cache_before':{k:v for k,v in before.items() if k!='cache_path'},
            'array_inventory':{k:{'shape':list(v.shape),'dtype':str(v.dtype)} for k,v in arrays.items()}}
        temporary=path.with_suffix('.npz.tmp')
        with temporary.open('wb') as f:np.savez_compressed(f,**arrays)
        temporary.replace(path);r['samples_sha256']=digest(path)
        reference.write_new_json(path.with_suffix('.json'),r)
        print('finished',job['name'],'cpu',round(cpu,2),'wall',round(wall,2),'mean',r['summary']['mean'],'SE',r['summary']['standard_error'],flush=True)
    return worker


def execute_stage(p,ids,deps,stage,resume):
    with ProcessPoolExecutor(max_workers=p['workers']) as pool:
        futures=[pool.submit(worker_main,p,ids,deps,stage,w,resume) for w in range(p['workers'])]
        for future in futures:future.result()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['prediction','validation','all'],default='all')
    parser.add_argument('--resume',action='store_true');args=parser.parse_args()
    p=read(PROTOCOL);OUT.mkdir(parents=True,exist_ok=True);ids=identity(p);design=ensure_design(p,ids)
    prepare_caches(p);deps={'execution_design_sha256':design}
    if args.stage in ('prediction','all'):execute_stage(p,ids,deps,'prediction',args.resume)
    register=lock_predictions(p,ids,design)
    if args.stage in ('validation','all'):
        execute_stage(p,ids,{**deps,'prediction_register_sha256':register},'validation',args.resume)
    print('completed independent response stage',args.stage,flush=True)


if __name__=='__main__':main()
