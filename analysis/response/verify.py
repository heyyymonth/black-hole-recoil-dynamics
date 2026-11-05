"""Read-only response evidence audit; NumPy only, with no dynamics rerun.

Checks every retained path and reconstructs stage inventories, sample-count
selection, seed coordinates, and prediction-before-validation dependencies.
Recorded UTC ordering corroborates the digest chain; it is local provenance,
not an externally timestamped preregistration.
"""
from datetime import datetime
from pathlib import Path
import argparse
import json
import math
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analysis.matched.verify import check_hashes, compare_tree, digest, require
from analysis.matched.statistics import _normal_critical
from analysis.response.analyze import analyze

OUT = ROOT/'results/response'
PROTOCOL = ROOT/'analysis/response/protocol.json'
CODE = (
    'analysis/response/run.py','analysis/response/transport.py','analysis/response/cache.py',
    'analysis/matched/run.py','analysis/matched/transport.py','analysis/matched/kernel.py',
    'analysis/matched/statistics.py','analysis/host_reconstruction/host.py',
    'analysis/encounters/physics.py','analysis/encounters/run.py',
)
BRANCH_ARRAYS = ('states','residence','normalized_residence','exposure',
    'normalized_exposure','end_time','stop_cause','alive','elastic_events','evaluated_steps')
BRANCHES = {'weighted':('jump_zero','jump_positive','diffusion_plus','diffusion_minus'),
    'naive':('baseline','replacement'),'paired_mixture':('baseline','mixture')}


def read(path):
    return json.loads(Path(path).read_text())


def timestamp(value):
    parsed = datetime.fromisoformat(value)
    require(parsed.tzinfo is not None and parsed.utcoffset().total_seconds()==0,
            'timestamps must explicitly use UTC')
    return parsed


def summary(values):
    return {'count':len(values),'mean':float(values.mean()),
        'variance':float(values.var(ddof=1)),
        'standard_error':float(values.std(ddof=1)/np.sqrt(len(values))),
        'minimum':float(values.min()),'maximum':float(values.max())}


def expected_jobs(p,stage,counts=None):
    """Inventory rebuilt without importing the SciPy-backed execution driver."""
    for si,setting in enumerate(p['numerical_settings']):
        if stage=='pilot' and setting['name']!=p['pilot']['setting']:
            continue
        for li,launch in enumerate(p['launches']):
            if stage=='pilot':
                variants=[(method,None,p['pilot']['pairs_per_method_per_launch'])
                          for method in p['pilot']['methods']]
            elif stage=='prediction':
                n=counts[launch['name']] if setting['role']=='primary' else p['prediction']['control_pairs']
                variants=[('weighted',None,n)]
            else:
                values=p['validation']['primary_a_values'] if setting['role']=='primary' else p['validation']['control_a_values']
                variants=[('paired_mixture',a,p['validation']['paired_paths_per_setting']) for a in values]
            for vi,(method,a,n) in enumerate(variants):
                suffix=f'a{a:.3f}' if stage=='validation' else method
                name=f'{stage}-{setting["name"]}-{launch["name"]}-{suffix}'
                yield {'name':name,'stage':stage,'setting_index':si,'launch_index':li,
                    'variant_index':vi,'setting':setting,'launch':launch,'method':method,
                    'a':a,'paths':int(n),'seed':[p['seed'],{'pilot':10,'prediction':20,'validation':30}[stage],si,li,vi]}


def count_selection(variance,p):
    z=_normal_critical(1-(1-p['reporting']['confidence'])/p['reporting']['bonferroni_family_size'])
    target=p['prediction']['primary_halfwidth_target']
    needed=max(float(p['prediction']['primary_pairs_min']),z*z*variance/(target*target))
    uncapped=2**math.ceil(math.log2(needed))
    maximum=p['prediction']['primary_pairs_max']
    return {'pilot_variance':variance,'uncapped_variance_rule_count':uncapped,
        'selected_count':min(maximum,uncapped),'cap_hit':uncapped>maximum},z


def verify_storage_adapter(p,ids,lock):
    """Audit the separately recorded, storage-only continuation provenance."""
    path=OUT/'storage-adapter.json'
    if not path.exists():
        return None
    adapter=read(path)
    require(adapter['revision']=='lossless-disk-storage-1','unknown storage adapter revision')
    require(adapter['frozen_identity']==ids,'storage adapter changed frozen execution identity')
    require(set(adapter['source_sha256'])=={
        'analysis/response/execute_cached.py','analysis/response/disk_cache.py',
        'analysis/response/test_disk_cache.py'},'incomplete storage adapter source inventory')
    require(set(adapter['input_sha256'])=={
        'results/response/execution-design.json','results/response/sample-count-lock.json'},
        'incomplete storage adapter input inventory')
    check_hashes(adapter['source_sha256']);check_hashes(adapter['input_sha256'])
    evidence=adapter['test_evidence']
    require(evidence['returncode']==0,'storage equivalence tests did not succeed')
    check_hashes({evidence['path']:evidence['sha256']})
    require((ROOT/evidence['path']).stat().st_size>0,'empty storage test evidence')
    check_hashes(adapter['existing_record_sha256'])
    counts={name:row['selected_count'] for name,row in lock['selection'].items()}
    jobs={job['name']:job for stage in ['pilot','prediction','validation']
          for job in expected_jobs(p,stage,counts)}
    existing={Path(name).stem for name in adapter['existing_record_sha256']}
    adapted=set(adapter['adapted_jobs'])
    require(len(existing)==len(adapter['existing_record_sha256'])==14,'storage continuation must preserve all 14 completed records')
    require(len(adapted)==len(adapter['adapted_jobs'])==22,'storage continuation has missing or duplicate adapted jobs')
    require(existing.isdisjoint(adapted) and existing|adapted==set(jobs),'storage usage does not partition the registered 36 jobs')
    require(all(str((OUT/(name+'.json')).relative_to(ROOT)) in adapter['existing_record_sha256'] for name in existing),'wrong existing-record paths')
    require(sum(jobs[name]['stage']=='prediction' for name in adapted)==1 and
            sum(jobs[name]['stage']=='validation' for name in adapted)==21,
            'adapter changed the registered continuation point')
    created=timestamp(adapter['created_utc'])
    require(timestamp(lock['created_utc'])<=created,'storage manifest precedes the fixed sample counts')
    for name in adapter['existing_record_sha256']:
        require(timestamp(read(ROOT/name)['completed_utc'])<=created,'storage manifest predates a declared completed record')
    return {'manifest':adapter,'created':created,'adapted':adapted,
        'reference':{'path':str(path.relative_to(ROOT)),'sha256':digest(path)}}


def verify_adapter_record(adapter,name,record):
    if adapter is not None and name in adapter['adapted']:
        require(record.get('storage_adapter')==adapter['reference'],'missing or stale declared storage adapter')
        require(adapter['created']<=timestamp(record['started_utc']),'adapted run precedes its storage manifest')
        require(isinstance(record.get('disk_cache'),dict),'adapted run lacks disk-cache diagnostics')
    else:
        require('storage_adapter' not in record and 'disk_cache' not in record,'unregistered adapter usage')


def verify_sample_arrays(a,r):
    """Independent reward, stopping, insertion, weighting, and censor accounting."""
    method=r['job']['method'];n=r['job']['paths'];H=r['horizon'];steps=r['n_steps']
    branches=BRANCHES[method]
    require(r['branches']==list(branches),'wrong branch inventory')
    extra=['numerical_missing_bound']
    if method=='paired_mixture':
        extra+=['paired_difference']
    else:
        extra+=['response_sample','replacement_steps','active_at_replacement','replacement_states']
    if method=='weighted':
        extra+=['jump_zero_weight','jump_positive_weight','selected_poisson_mean','positive_branch_count']
    expected={f'{branch}__{key}' for branch in branches for key in BRANCH_ARRAYS}|set(extra)
    require(set(a)==expected,'unexpected or missing retained arrays')
    for key,value in a.items():
        shape=(n,6) if key.endswith('__states') or key=='replacement_states' else (n,)
        require(value.shape==shape,'unaligned path array '+key)
        if value.dtype.kind in 'fiu':
            require(np.isfinite(value).all(),'nonfinite path array '+key)
    compare_tree(r['array_inventory'],{key:{'shape':list(value.shape),'dtype':str(value.dtype)} for key,value in a.items()})
    bounds={};stops={}
    for branch in branches:
        b={key:a[f'{branch}__{key}'] for key in BRANCH_ARRAYS}
        end=b['end_time'];causes=b['stop_cause'];alive=b['alive']
        require(alive.dtype.kind=='b' and causes.dtype.kind in 'SU','invalid survival data types')
        require(np.isin(causes,['horizon','spatial_domain','velocity_domain','excluded_channel']).all(),'unregistered stop')
        require(np.array_equal(alive,causes=='horizon'),'inconsistent survival status')
        require(np.all((end>=0)&(end<=H+1e-12)),'invalid modeled duration')
        require(np.all(end[alive]==H),'surviving path has wrong horizon')
        require(np.all((b['residence']>=-1e-12)&(b['residence']<=end+1e-10)),'residence exceeds modeled duration')
        require(np.all((b['exposure']>=0)&(b['exposure']<=end+1e-7)),'normalized density exposure exceeds modeled duration')
        np.testing.assert_allclose(b['normalized_residence'],b['residence']/H,rtol=0,atol=1e-13)
        np.testing.assert_allclose(b['normalized_exposure'],b['exposure']/H,rtol=0,atol=1e-13)
        for key in ['elastic_events','evaluated_steps']:
            require(b[key].dtype.kind in 'iu' and np.all(b[key]>=0),'invalid integral count '+key)
        require(np.all(b['evaluated_steps']<=steps),'too many evaluated steps')
        if branch in ['baseline','diffusion_plus','diffusion_minus','jump_zero']:
            require(not b['elastic_events'].any(),'events in a zero-jump branch')
        spatial=causes=='spatial_domain';velocity=causes=='velocity_domain'
        radius=np.linalg.norm(b['states'][:,:3],axis=1);speed=np.linalg.norm(b['states'][:,3:],axis=1)
        np.testing.assert_allclose(radius[spatial],r['spatial_stop'],rtol=0,atol=1e-10)
        require(np.all(speed[velocity]>=r['velocity_stop']-1e-12),'velocity stop below threshold')
        modeled=~(spatial|velocity)
        require(np.all(radius[modeled]<r['spatial_stop']+1e-12),'modeled path outside spatial domain')
        require(np.all(speed[modeled]<r['velocity_stop']+1e-12),'modeled path outside velocity domain')
        bounds[branch]=np.where(spatial|velocity,np.maximum(0.,1-end/H),0.)
        labels,counts=np.unique(causes,return_counts=True)
        stops[branch]=dict(zip(labels.tolist(),counts.tolist()))
    compare_tree(r['stops'],stops)
    reward=lambda name:a[f'{name}__normalized_residence']
    if method=='weighted':
        p0,pp,mu=[a[key] for key in ['jump_zero_weight','jump_positive_weight','selected_poisson_mean']]
        require(np.all(mu>=0),'negative Poisson mean')
        require(np.all((p0>=0)&(p0<=1)&(pp>=0)&(pp<=1)),'invalid mixture weights')
        np.testing.assert_allclose(p0,np.exp(-mu),rtol=2e-13,atol=2e-15)
        np.testing.assert_allclose(pp,-np.expm1(-mu),rtol=2e-13,atol=2e-15)
        np.testing.assert_allclose(p0+pp,1.,rtol=0,atol=2e-15)
        count=a['positive_branch_count']
        require(count.dtype.kind in 'iu' and np.all(count[mu>0]>0) and np.all(count[mu==0]==0),'invalid conditional positive counts')
        applied=a['jump_positive__elastic_events']
        require(np.all((applied==0)|(applied==count)),'positive branch applies a partial or excess event batch')
        score=steps*(p0*reward('jump_zero')+pp*reward('jump_positive')-.5*(reward('diffusion_plus')+reward('diffusion_minus')))
        bound=steps*(p0*bounds['jump_zero']+pp*bounds['jump_positive']+.5*(bounds['diffusion_plus']+bounds['diffusion_minus']))
    elif method=='naive':
        score=steps*(reward('replacement')-reward('baseline'))
        bound=steps*(bounds['replacement']+bounds['baseline'])
    else:
        score=reward('mixture')-reward('baseline');bound=bounds['mixture']+bounds['baseline']
    score_key='paired_difference' if method=='paired_mixture' else 'response_sample'
    np.testing.assert_allclose(a[score_key],score,rtol=2e-11,atol=2e-13)
    np.testing.assert_allclose(a['numerical_missing_bound'],bound,rtol=2e-11,atol=2e-13)
    compare_tree(r['summary'],summary(score));compare_tree(r['numerical_missing_bound'],summary(bound))
    if method!='paired_mixture':
        chosen=a['replacement_steps'];active=a['active_at_replacement'];state=a['replacement_states']
        require(chosen.dtype.kind in 'iu' and np.all((chosen>=0)&(chosen<steps)),'invalid insertion step')
        require(active.dtype.kind=='b','invalid active insertion flags')
        require(np.all(state[~active]==0),'inactive insertion lacks zero placeholder')
        require(np.all(abs(score[~active])<2e-10),'dead-before-insertion response is nonzero')
        require(np.all(np.linalg.norm(state[active,:3],axis=1)<r['spatial_stop']),'insertion outside spatial domain')
        require(np.all(np.linalg.norm(state[active,3:],axis=1)<r['velocity_stop']),'insertion outside speed domain')
        for branch in branches:
            require(np.all(a[f'{branch}__evaluated_steps'][active]>=chosen[active]+1),'active insertion was never evaluated')
            require(np.all(a[f'{branch}__evaluated_steps'][~active]<=chosen[~active]),'inactive insertion was already passed')
            for key in BRANCH_ARRAYS:
                np.testing.assert_array_equal(a[f'{branch}__{key}'][~active],a[f'{branches[0]}__{key}'][~active])
        if method=='weighted':
            require(np.all(mu[~active]==0),'inactive insertion has nonzero event rate')


def verify():
    from analysis.response.coupling_diagnostic import verify as verify_coupling
    verify_coupling(read(OUT/'coupling-diagnostic.json'))
    p=read(PROTOCOL);matched=read(ROOT/p['matched_protocol'])
    design_path=OUT/'execution-design.json';lock_path=OUT/'sample-count-lock.json';register_path=OUT/'prediction-register.json'
    design=read(design_path);lock=read(lock_path);register=read(register_path)
    ids=design['identity']
    require(ids['protocol_sha256']==digest(PROTOCOL),'stale response protocol')
    required_inputs={str(PROTOCOL.relative_to(ROOT)),p['matched_protocol'],p['host_path'],
        matched['encounter_protocol'],matched['baseline_protocol'],'results/collisionless-residence/NGC5139.json'}
    for path in p['kernel_tables'].values():
        required_inputs.update([path,str(Path(path).with_suffix('.json'))])
    require(set(ids['input_sha256'])==required_inputs,'incomplete physical input provenance')
    require(set(ids['implementation_sha256'])==set(CODE),'incomplete implementation provenance')
    check_hashes(ids['input_sha256']);check_hashes(ids['implementation_sha256'])
    require(all(isinstance(ids[k],str) and ids[k] for k in ['numpy_version','scipy_version']),'missing execution versions')
    for item in [design,lock,register]:
        require(item['revision']==p['revision'] and item['identity']==ids,'inconsistent stage identity')
    require(design['seed_scheme']==p['seed_scheme'],'wrong seed scheme')
    require(design['stage_order']==['pilot','sample-count-lock','prediction','prediction-register','validation'],'wrong stage ordering')
    design_sha,count_sha,register_sha=map(digest,[design_path,lock_path,register_path])
    require(lock['execution_design_sha256']==register['execution_design_sha256']==design_sha,'stale execution design dependency')
    require(register['sample_count_lock_sha256']==count_sha,'stale sample-count dependency')
    times=[timestamp(item['created_utc']) for item in [design,lock,register]]
    require(times[0]<=times[1]<=times[2],'stage locks out of order')
    adapter=verify_storage_adapter(p,ids,lock)
    counts={};records={};all_seeds=[]
    table_meta={grid:read((ROOT/path).with_suffix('.json')) for grid,path in p['kernel_tables'].items()}
    # Reuse the earlier hash-bound launch normalization; this audit does not
    # rebuild a gravitational potential or silently approximate its spline.
    reference=read(ROOT/'results/matched-residence/fine-dt0.020-radial_outward-jump.json')
    require(reference['host_sha256']==ids['input_sha256'][p['host_path']],'launch reference host mismatch')
    require(reference['protocol_sha256']==ids['input_sha256'][p['matched_protocol']],'launch reference protocol mismatch')
    speed=float(np.linalg.norm(reference['initial_state'][3:]))
    for stage in ['pilot','prediction','validation']:
        inventory=[]
        dependencies={'execution_design_sha256':design_sha}
        if stage!='pilot':dependencies['sample_count_lock_sha256']=count_sha
        if stage=='validation':dependencies['prediction_register_sha256']=register_sha
        for job in expected_jobs(p,stage,counts):
            path=OUT/(job['name']+'.npz');meta=path.with_suffix('.json');r=read(meta)
            require(r['complete'] is True and r['revision']==p['revision'],'incomplete or wrong result revision')
            compare_tree(r['job'],job)
            require(r['identity']==ids and r['stage_dependencies']==dependencies,'stale run source/stage identity')
            require(r['samples_sha256']==digest(path),'changed retained arrays')
            start,end=timestamp(r['started_utc']),timestamp(r['completed_utc'])
            verify_adapter_record(adapter,job['name'],r)
            lower=times[['pilot','prediction','validation'].index(stage)]
            require(lower<=start<=end,'run precedes its stage lock or has negative duration')
            if stage!='validation':
                require(end<=times[['pilot','prediction'].index(stage)+1],'stage lock precedes a required completed run')
            require(np.isfinite(r['wall_seconds']) and r['wall_seconds']>0,'invalid run cost')
            require(r['horizon']==p['horizon'] and r['dt']==job['setting']['dt'],'wrong time controls')
            require(r['n_steps']==int(np.ceil(p['horizon']/job['setting']['dt'])),'wrong scheduled step count')
            require(r['region_radius']==matched['region_radius_over_density_radius'],'wrong residence region')
            scale=table_meta[job['setting']['grid']]['scales']
            compare_tree(r['scales'],{'density_radius_pc':scale['r_scale_pc'],'velocity_kms':scale['v_scale_kms'],'time_Myr':scale['time_scale_Myr']})
            inherited=read(ROOT/'results/collisionless-residence/NGC5139.json')
            force_domain=min(row['summary']['stop_radius_pc'] for row in inherited['records'])/scale['r_scale_pc']
            compare_tree(r['inherited_force_domain'],force_domain)
            compare_tree(r['spatial_stop'],min(matched['kernel']['coarse_axes']['radius'][-1],force_domain))
            compare_tree(r['velocity_stop'],float(matched['kernel']['coarse_axes']['speed'][-1]))
            angle=np.deg2rad(job['launch']['angle_from_radial_degrees'])
            initial=[matched['launch']['radius_over_density_radius'],0.,0.,speed*np.cos(angle),speed*np.sin(angle),0.]
            compare_tree(r['initial_state'],initial)
            with np.load(path,allow_pickle=False) as retained:a=dict(retained)
            verify_sample_arrays(a,r)
            records[job['name']]=r;all_seeds.append(tuple(job['seed']))
            inventory.append({'name':job['name'],'arrays':str(path.relative_to(ROOT)),
                'arrays_sha256':digest(path),'record':str(meta.relative_to(ROOT)),
                'record_sha256':digest(meta),'summary':r['summary']})
            if stage=='pilot' and job['method']=='weighted':
                selection,z=count_selection(float(a['response_sample'].var(ddof=1)),p)
                compare_tree(lock['selection'][job['launch']['name']],selection)
                compare_tree(lock['normal_critical'],z)
                counts[job['launch']['name']]=selection['selected_count']
        if stage=='pilot':compare_tree(lock['pilot_inventory'],inventory)
        if stage=='prediction':compare_tree(register['prediction_inventory'],inventory)
    require(lock['pilot_excluded_from_production'] is True,'pilot included in inference')
    compare_tree(lock['target_halfwidth'],p['prediction']['primary_halfwidth_target'])
    require(len(set(all_seeds))==36 and len(records)==36,'seed collision or incomplete experiment inventory')
    expected=set(records)
    for extension in ['npz','json']:
        actual={f.stem for stage in ['pilot','prediction','validation'] for f in OUT.glob(f'{stage}-*.{extension}')}
        if extension=='json':actual.discard('prediction-register')
        require(actual==expected,'extra or missing result inventory')
    result=analyze();compare_tree(read(OUT/'summary.json'),result)
    require(result['inferential_comparisons']==63<=result['bonferroni_family_size'],'inferential family not covered')
    print(f"Verified 6 pilot, 9 prediction and 21 independent-validation records; {result['counts']['branch_trajectories']} retained branch trajectories. No trajectories rerun.")
    return result


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args();verify()
