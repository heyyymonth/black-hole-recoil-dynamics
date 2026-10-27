"""Repeat the registered local experiments across six conditional cluster hosts.

References use the original maximum-quadrature-weight rule. A separate fixed-seed
weighted joint-row pilot tests sensitivity; it is not a converged posterior study.
Every attempted row is retained, including failures and legitimate zero-BH hosts.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis/host_reconstruction'))
from host import Host

CLUSTERS=['NGC0104','NGC1851','NGC3201','NGC5139','NGC6121','NGC6397']
OUT=ROOT/'results/cluster-comparison'
ANALYSIS=ROOT/'analysis/cluster_comparison'
SEED=20261026
DRAWS=8
CORE_RADII=[0.,.001,.003,.01,.03,.1,.3,1.,3.]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def csv_save(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def selection(logw, cluster_index, count=DRAWS):
    if np.isnan(logw).any() or np.isposinf(logw).any() or not np.isfinite(logw).any():
        raise ValueError('Invalid nested quadrature log weights')
    w=np.exp(logw-np.max(logw));w/=w.sum()
    rng=np.random.default_rng(np.random.SeedSequence([SEED,cluster_index]))
    return w,rng.choice(len(w),size=count,replace=True,p=w)

def prepare():
    import h5py
    sources=json.loads((ROOT/'data/observed-clusters/sources.json').read_text())
    pinned={x['file']:x for x in sources['files']}
    base=json.loads((ROOT/'analysis/host_reconstruction/protocol.json').read_text())
    plan={'revision':'six-cluster-local-comparison-1','clusters':CLUSTERS,
          'reference_selection':base['joint_sample_selection'],
          'pilot_draws_per_cluster':DRAWS,'seed':SEED,
          'pilot_selection':'iid categorical normalized nested/logwt; SeedSequence([seed, cluster list index]); retain entire joint row; repetitions retained; no replacement of failed rows',
          'pilot_interpretation':'Eight weighted draws per cluster are a sensitivity pilot, not converged credible intervals or posterior failure probabilities.',
          'reference_states':'Same 20 r/v/orientation states and 60-Msun probe as the original calculation',
          'radial_zero_speed_profile':CORE_RADII,
          'zero_speed_quadrature_orders':[64,128],
          'zero_speed_relative_refinement_limit':.005,
          'profile_interpolation':'Log cubic on monotone intervals; analytic derivative extrema select PCHIP fallback on intervals reversing the nodal trend; positive lowering-law outer taper.',
          'projection_refinement':'Starting with 128/256 nodes, double to at most 2048 when refinement fails; retain original difference and unchanged 0.2-percent threshold.',
          'normalization':'Per unit identical constant ln Lambda=1 for every species; not calibrated collision rates',
          'scope':'Host validation, local scattering and sparse-reservoir diagnostics; no remnant trajectories or residence',
          'models':'Adopted published standard spherical nonrotating multimass fits; NGC6397 uses the available standard-fit chain. Its distinct no-BH alternative has no chain in the inspected release and is not evaluated.',
          'not_claimed':['unique observed BH population','posterior-calibrated transport','diffusion failure','subsequent merger rates']}
    save(ANALYSIS/'protocol.json',plan)
    tasks=[]
    for ci,c in enumerate(CLUSTERS):
        src={}
        for suffix in ['sampler.hdf','observations.hdf']:
            name=f'gcfit_{c}_{suffix}';path=ROOT/'data/observed-clusters/raw'/name
            if sha(path)!=pinned[name]['sha256']:raise ValueError('Source hash mismatch '+name)
            src[str(path.relative_to(ROOT))]=pinned[name]['sha256']
        p={**base,'cluster':c,'sources':src}
        pp=ANALYSIS/'protocols'/f'{c}.json';save(pp,p)
        with h5py.File(ROOT/f'data/observed-clusters/raw/gcfit_{c}_sampler.hdf') as f:
            w,draws=selection(f['nested/logwt'][:],ci)
            rows=f['nested/samples'][:]
        records=[{'kind':'reference','label':'reference','row_index':int(w.argmax()),'normalized_source_weight':float(w.max())}]
        records += [{'kind':'pilot','label':f'pilot-{j:02d}','draw_index':j,'row_index':int(i),'normalized_source_weight':float(w[i])} for j,i in enumerate(draws)]
        for rec in records:
            rec['cluster']=c;rec['theta']=dict(zip(base['parameter_order'],rows[rec['row_index']].tolist()))
            rec['protocol_path']=str(pp.relative_to(ROOT));tasks.append(rec)
        save(OUT/c/'selection.json',{'cluster':c,'source_sha256':src,'weight_ESS':float(1/(w@w)), 'rows':records})
    return tasks

def zero_speed(host,r,order):
    cutoff=float(np.sqrt(2*host.psi(r)))
    if cutoff<=0:raise ValueError('Radius outside field distribution support')
    x,w=leggauss(order);v=(x+1)*cutoff/2
    vv,mu=np.meshgrid(v,x,indexing='ij')
    u=np.column_stack([(vv*np.sqrt(1-mu*mu)).ravel(),np.zeros(vv.size),(vv*mu).ravel()])
    # h_j(0)=int F_j(u)/|u| d^3u; angular azimuth integrates to 2pi.
    weight=(2*np.pi*(w*cutoff/2*v)[:,None]*w[None,:]).ravel()
    h=host.number_df(float(r),u)@weight
    return 8*np.pi*host.G**2*host.mj**2*h/.9777922216807892

def core_profile(host):
    rows=[];bh=(host.types=='BH')&host.physical
    total_bh_mass=float(host.Mj[bh].sum())
    count=float((host.Mj[bh]/host.mj[bh]).sum())
    for fraction in CORE_RADII:
        r=fraction*host.rh
        if r>=host.rt:
            rows.append({'r_over_rh':fraction,'r_pc':r,'status':'outside_DF_support','trace_D_per_lnLambda':None,'BH_trace_D_fraction':None,'BH_density_fraction':None,'refinement_relative_error':None})
            continue
        coarse,fine=[zero_speed(host,r,n) for n in [64,128]]
        denom=float(fine[host.physical].sum())
        err=float(np.sum(np.abs(fine-coarse)[host.physical])/denom)
        rows.append({'r_over_rh':fraction,'r_pc':r,'status':'passed' if err<.005 else 'quadrature_failed',
                     'trace_D_per_lnLambda':denom,'BH_trace_D_fraction':float(fine[bh].sum()/denom),
                     'BH_density_fraction':float(host.rho_spline(r)[bh].sum()/host.rho_spline(r)[host.physical].sum()),
                     'refinement_relative_error':err})
    return rows,{'BH_expected_number':count,'BH_mass_Msun':total_bh_mass,
                 'probe_to_BH_reservoir_mass_ratio':60/total_bh_mass if total_bh_mass>0 else None,
                 'BH_absent':total_bh_mass==0,
                 'BH_mean_individual_mass_Msun':total_bh_mass/count if count>0 else None}

def worker(task):
    started=time.monotonic();cluster,label=task['cluster'],task['label']
    export=OUT/cluster/label;export.mkdir(parents=True,exist_ok=True)
    # Large pilot numerical meshes are rebuildable scratch, not 48 duplicate archives.
    work=export if task['kind']=='reference' else ROOT/'tmp/cluster-comparison'/cluster/label
    work.mkdir(parents=True,exist_ok=True)
    log=ROOT/'tmp/cluster-comparison/logs'/f'{cluster}-{label}.log';log.parent.mkdir(parents=True,exist_ok=True)
    previous_result=json.loads((export/'result.json').read_text()) if (export/'result.json').exists() else None
    result={**task,'host_status':'not_run','transport_status':'not_run'}
    if previous_result and previous_result.get('previous_failed_attempt'):
        result['previous_failed_attempt']=previous_result['previous_failed_attempt']
    if previous_result and (previous_result.get('host_status')=='failed' or previous_result.get('transport_status')=='failed'):
        result['previous_failed_attempt']={k:previous_result[k] for k in ['host_status','host_summary','failure_log','transport_status','transport_summary','transport_failure_log'] if k in previous_result}
    cmd=[sys.executable,str(ROOT/'analysis/host_reconstruction/run.py'),'--protocol',str(ROOT/task['protocol_path']), '--output',str(work),'--row-index',str(task['row_index'])]
    cached=work/'summary.json'
    if cached.exists():
        prev=json.loads(cached.read_text())
        if prev.get('theta')==task['theta'] and prev.get('protocol_sha256')==sha(ROOT/task['protocol_path']):cmd.append('--reuse-solution')
    with log.open('w') as f:
        done=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    if cached.exists():
        hs=json.loads(cached.read_text());result['host_summary']=hs
        if task['kind']=='pilot':
            save(export/'summary.json',hs)
            for name in ['species.csv','df-checks.csv','projection-checks.csv']:
                if (work/name).exists():(export/name).write_bytes((work/name).read_bytes())
    if done.returncode:
        result.update(host_status='failed',failure_log=log.read_text()[-6000:].replace(str(ROOT),'<repo>'))
    else:
        result['host_status']='passed'
        try:
            host=Host(dict(np.load(work/'host.npz',allow_pickle=False)))
            profile,metrics=core_profile(host);csv_save(export/'core-profile.csv',profile)
            result['BH_metrics']=metrics;result['core_profile']=profile
            result['core_profile_status']='passed' if all(r['status']=='passed' for r in profile) else 'partial_or_failed'
            result['host_arrays_sha256']=sha(work/'host.npz')
            if task['kind']=='reference':
                with log.open('a') as f:
                    done=subprocess.run([sys.executable,str(ROOT/'analysis/transport/run.py'),'--host',str(work),'--output',str(export/'transport')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
                result['transport_status']='passed' if done.returncode==0 else 'failed'
                if (export/'transport/summary.json').exists():result['transport_summary']=json.loads((export/'transport/summary.json').read_text())
                if done.returncode:result['transport_failure_log']=log.read_text()[-4000:].replace(str(ROOT),'<repo>')
        except (ValueError,FloatingPointError) as exc:
            result['local_calculation_error']=str(exc)
    result['runtime_seconds']=round(time.monotonic()-started,3)
    save(export/'result.json',result)
    print(cluster,label,result['host_status'],result['transport_status'],result['runtime_seconds'],flush=True)
    return result

def collect():
    rows=[]
    for c in CLUSTERS:
        selection=json.loads((OUT/c/'selection.json').read_text())
        for rec in selection['rows']:
            path=OUT/c/rec['label']/'result.json'
            rows.append(json.loads(path.read_text()) if path.exists() else {**rec,'host_status':'pending','transport_status':'pending'})
    save(OUT/'summary.json',{'revision':'six-cluster-local-comparison-1','runs':rows,
         'attempted':sum(r['host_status']!='pending' for r in rows),'planned':len(rows),
         'remnant_trajectories_performed':False,'pilot_is_converged_posterior':False})

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['reference','pilot','all','collect','retry-failed'],default='all');parser.add_argument('--workers',type=int,default=2);parser.add_argument('--resume',action='store_true');parser.add_argument('--clusters',nargs='+',choices=CLUSTERS,default=CLUSTERS)
    args=parser.parse_args();tasks=prepare()
    if args.stage!='collect':
        selected=[t for t in tasks if t['cluster'] in args.clusters and (args.stage in ['all','retry-failed'] or t['kind']==args.stage)]
        if args.resume:
            selected=[t for t in selected if not (OUT/t['cluster']/t['label']/'result.json').exists()]
        if args.stage=='retry-failed':
            selected=[t for t in selected if (OUT/t['cluster']/t['label']/'result.json').exists() and json.loads((OUT/t['cluster']/t['label']/'result.json').read_text()).get('host_status')=='failed']
        with ProcessPoolExecutor(max_workers=args.workers) as pool:list(pool.map(worker,selected))
    collect()
if __name__=='__main__':main()
