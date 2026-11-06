"""Read-only audit of the independently sampled symmetric-root response retest.

Reconstruct the fixed inventory and random-stream assignments independently
of the execution driver. Verify complete transition rewards, retained scores,
provenance, and prediction registration before validation. Parallel workers
may overlap; work assigned to the same worker may not overlap.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.response.verify import compare_tree, timestamp, verify_sample_arrays

OUT = ROOT/'results/retest'
PROTOCOL = ROOT/'analysis/retest/protocol.json'


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_record_arrays(arrays, record):
    """Recompute every reward/weight/stop and the scheduled insertion draws."""
    verify_sample_arrays(arrays, record)
    if record['job']['method']=='weighted':
        rng = np.random.default_rng(np.random.SeedSequence(record['job']['seed']).spawn(5)[4])
        expected = rng.integers(0, record['n_steps'], size=record['job']['paths'])
        np.testing.assert_array_equal(arrays['replacement_steps'], expected)


def verify_chronology(records, design, register):
    """Check stage gates while allowing independent workers to overlap."""
    designed = timestamp(design['created_utc'])
    registered = timestamp(register['created_utc'])
    require(designed<=registered, 'prediction register predates execution design')
    workers = {}
    for record in records:
        job = record['job']
        start, end = timestamp(record['started_utc']), timestamp(record['completed_utc'])
        require(designed<=start<=end, 'run predates design or has negative duration')
        if job['stage']=='prediction':
            require(end<=registered, 'prediction register precedes a completed prediction')
        elif job['stage']=='validation':
            require(registered<=start, 'validation precedes the prediction register')
        else:
            raise AssertionError('unknown retest stage')
        workers.setdefault(job['worker'], []).append((start, end, job['name']))
    for timeline in workers.values():
        ordered = sorted(timeline)
        for first, second in zip(ordered[:-1], ordered[1:]):
            require(first[1]<=second[0], 'overlapping runs assigned to one worker')


def verify_job(record, expected, identity, dependencies):
    compare_tree(record['job'], expected)
    require(record['identity']==identity, 'changed run identity')
    require(record['stage_dependencies']==dependencies, 'changed stage dependency')


def expected_jobs(p):
    jobs=[]
    for stage in ('prediction','validation'):
        index=0
        for si,s in enumerate(p['numerical_settings']):
            for li,l in enumerate(p['launches']):
                values=[None] if stage=='prediction' else (p['validation']['primary_a_values'] if s['role']=='primary' else p['validation']['control_a_values'])
                method='weighted' if stage=='prediction' else 'paired_mixture'
                n=(p['prediction']['primary_quartets'] if s['role']=='primary' else p['prediction']['control_quartets']) if stage=='prediction' else p['validation']['pairs_per_setting']
                for vi,a in enumerate(values):
                    suffix='weighted' if stage=='prediction' else f'a{a:.3f}'
                    jobs.append({'name':f"{stage}-{s['name']}-{l['name']}-{suffix}",
                        'stage':stage,'setting_index':si,'launch_index':li,'variant_index':vi,
                        'setting':s,'launch':l,'method':method,'a':a,'paths':n,
                        'seed':[p['seed'],50 if stage=='prediction' else 60,si,li,vi],
                        'worker':index%p['workers']})
                    index+=1
    return jobs


def verify():
    p=read(PROTOCOL);design=read(OUT/'execution-design.json');register=read(OUT/'prediction-register.json')
    ids=design['identity'];jobs=expected_jobs(p)
    require(len(jobs)==24 and sum(j['stage']=='prediction' for j in jobs)==9,'wrong fixed inventory')
    compare_tree(design['jobs'],jobs)
    require(design['revision']==p['revision'] and design['scope']==p['scope'],'wrong design revision')
    require(design['stage_order']==['prediction','prediction-register','validation'],'wrong stage order')
    required={str(PROTOCOL.relative_to(ROOT)),p['reference_design'],p['precision_design'],
        p['development_evidence'],'results/response/summary.json','results/retest/preflight.txt',
        'analysis/retest/run.py','analysis/retest/analyze.py','analysis/retest/statistics.py',
        'analysis/precision/coupling.py','analysis/response/disk_cache.py'}
    require(set(ids['source_sha256'])==required,'incomplete source inventory')
    for path,sha in ids['source_sha256'].items():require(digest(ROOT/path)==sha,'changed source '+path)
    reference=read(ROOT/p['reference_design'])['identity']
    require(ids['reference_identity']==reference,'changed numerical reference')
    for group in ('input_sha256','implementation_sha256'):
        for path,sha in reference[group].items():require(digest(ROOT/path)==sha,'changed reference source '+path)
    precision=read(ROOT/p['precision_design'])['identity']
    require(precision['reference_identity']==reference,'different precision reference')
    for path,sha in precision['source_sha256'].items():require(digest(ROOT/path)==sha,'changed precision source '+path)
    previous=read(ROOT/'results/response/summary.json')['input_sha256']
    require(register['revision']==p['revision'] and register['identity']==ids,'changed prediction registration')
    require(register['execution_design_sha256']==digest(OUT/'execution-design.json'),'changed design link')
    records=[];inventory=[];seeds=[]
    for job in jobs:
        path=OUT/(job['name']+'.npz');meta=path.with_suffix('.json');r=read(meta)
        deps={'execution_design_sha256':digest(OUT/'execution-design.json')}
        if job['stage']=='validation':deps['prediction_register_sha256']=digest(OUT/'prediction-register.json')
        verify_job(r,job,ids,deps);seeds.append(tuple(job['seed']))
        require(r['complete'] and r['revision']==p['revision'],'incomplete result')
        require(r['samples_sha256']==digest(path),'changed retained arrays')
        require(np.isfinite([r['cpu_seconds'],r['wall_seconds']]).all() and min(r['cpu_seconds'],r['wall_seconds'])>0,'invalid timing')
        require(r['horizon']==p['horizon'] and r['dt']==job['setting']['dt'],'changed time controls')
        require(r['n_steps']==int(np.ceil(r['horizon']/r['dt'])),'wrong step count')
        prior_path=f"results/response/prediction-{job['setting']['name']}-{job['launch']['name']}-weighted.json"
        require(digest(ROOT/prior_path)==previous[prior_path],'changed previous physical controls')
        prior=read(ROOT/prior_path)
        for key in ('initial_state','spatial_stop','velocity_stop','inherited_force_domain','region_radius','scales'):
            compare_tree(r[key],prior[key])
        require(r['coupling'].startswith('Symmetric PSD root; '),'wrong coupling')
        with np.load(path,allow_pickle=False) as loaded:arrays=dict(loaded)
        verify_record_arrays(arrays,r);records.append(r)
        if job['stage']=='prediction':
            inventory.append({'name':job['name'],'arrays':str(path.relative_to(ROOT)),
                'arrays_sha256':digest(path),'record':str(meta.relative_to(ROOT)),
                'record_sha256':digest(meta),'summary':r['summary']})
    require(len(seeds)==len(set(seeds))==24,'seed reuse')
    compare_tree(register['prediction_inventory'],inventory)
    verify_chronology(records,design,register)
    names={j['name'] for j in jobs}
    for suffix in ('json','npz'):
        actual={f.stem for stage in ('prediction','validation') for f in OUT.glob(stage+'-*.'+suffix) if f.stem!='prediction-register'}
        require(actual==names,'unexpected or missing outcome inventory')
    from analysis.retest.analyze import analyze
    result=analyze();compare_tree(read(OUT/'summary.json'),result)
    require(result['counts']=={'prediction_ensembles':9,'validation_ensembles':15,
        'prediction_quartets':245760,'validation_pairs':122880,'independent_groups':368640,
        'branch_trajectories':1228800,'normal_statements':51,'continuation_bounds':24},'wrong sample totals')
    print('Verified independent response retest: 9 registered predictions, 15 fresh validations, 368640 independent groups and 24 population continuation bounds; no trajectories rerun.')
    return result


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args();verify()
