"""Verify retained precision evidence with NumPy; no SciPy or dynamics required."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.precision.analyze import analyze, read, digest, OUT, PROTOCOL
from analysis.response.verify import compare_tree, verify_sample_arrays, timestamp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def verify_record_arrays(a, record, p):
    job = record['job']; labels = a['insertion_stratum']
    require(labels.dtype.kind in 'iu' and labels.shape==(job['paths'],), 'invalid stratum labels')
    if job['method']['selection']=='uniform':
        require(np.all(labels==-1), 'uniform rows have stratum labels')
        rng = np.random.default_rng(np.random.SeedSequence(job['seed']).spawn(5)[4])
        expected = rng.integers(0, record['n_steps'], size=job['paths'])
    else:
        per = job['paths']//4
        require(job['paths']%4==0, 'nonproportional allocation')
        np.testing.assert_array_equal(labels, np.repeat(np.arange(4), per))
        rng = np.random.default_rng(np.random.SeedSequence(job['selection_seed']))
        expected = np.concatenate([rng.integers(lo, hi, size=per)
            for lo,hi in zip(p['step_strata'][:-1],p['step_strata'][1:])])
    np.testing.assert_array_equal(a['replacement_steps'],expected)
    old_arrays = {k:v for k,v in a.items() if k!='insertion_stratum'}
    old_record = {**record, 'job': {**job, 'method':'weighted'},
        'array_inventory': {k:v for k,v in record['array_inventory'].items() if k!='insertion_stratum'}}
    verify_sample_arrays(old_arrays, old_record)
    compare_tree(record['array_inventory']['insertion_stratum'],
                 {'shape':list(labels.shape),'dtype':str(labels.dtype)})


def verify():
    p=read(PROTOCOL);design=read(OUT/'execution-design.json');ids=design['identity']
    require(design['revision']==p['revision'] and design['scope']==p['scope'],'wrong design')
    for path,sha in ids['source_sha256'].items():
        require(digest(ROOT/path)==sha,'changed precision source '+path)
    required={'analysis/precision/protocol.json','results/response/execution-design.json',
        'results/response/summary.json','analysis/precision/run.py','analysis/precision/coupling.py',
        'analysis/precision/statistics.py','analysis/precision/analyze.py','analysis/response/disk_cache.py',
        'results/precision/preflight.txt'}
    require(set(ids['source_sha256'])==required,'incomplete source inventory')
    ref=read(ROOT/p['reference_design'])['identity']
    require(ids['reference_identity']==ref,'different reference law')
    for group in ('input_sha256','implementation_sha256'):
        for path,sha in ref[group].items():require(digest(ROOT/path)==sha,'changed frozen source '+path)
    expected=[];seeds=[];start_time=timestamp(design['created_utc'])
    reference=read(ROOT/'results/response/prediction-fine-dt0.010-radial_outward-weighted.json')
    for block in range(p['blocks']):
        for offset in range(len(p['launches'])):
            li=(offset+block)%len(p['launches'])
            for mo in range(len(p['methods'])):
                mi=(mo+block+li)%len(p['methods']);method=p['methods'][mi];launch=p['launches'][li]
                job={'name':f"{launch}-{method['name']}-block{block}",'launch':launch,
                    'launch_index':li,'method':method,'method_index':mi,'block':block,
                    'paths':p['quartets_per_block'],'seed':[p['seed'],40,li,mi,block],
                    'selection_seed':[p['seed'],41,li,mi,block]}
                expected.append(job);seeds.extend([tuple(job['seed']),tuple(job['selection_seed'])])
                path=OUT/(job['name']+'.npz');r=read(path.with_suffix('.json'))
                require(r['complete'] and r['revision']==p['revision'],'incomplete result')
                compare_tree(r['job'],job)
                require(r['identity']==ids and r['design_sha256']==digest(OUT/'execution-design.json'),'stale identity')
                require(r['samples_sha256']==digest(path),'changed sample bytes')
                begin,end=timestamp(r['started_utc']),timestamp(r['completed_utc'])
                require(start_time<=begin<=end,'result predates design')
                start_time=end
                require(np.isfinite([r['cpu_seconds'],r['wall_seconds']]).all() and min(r['cpu_seconds'],r['wall_seconds'])>0,'invalid cost')
                for key in ('horizon','dt'):require(r[key]==p[key],'changed time controls')
                require(r['n_steps']==int(np.ceil(p['horizon']/p['dt'])),'wrong step count')
                require(r['n_steps']==p['step_strata'][-1] and p['step_strata'][0]==0,'incomplete strata')
                np.testing.assert_allclose(np.diff(p['step_strata'])/r['n_steps'],p['stratum_weights'],rtol=0,atol=0)
                prior=read(ROOT/f'results/response/prediction-fine-dt0.010-{launch}-weighted.json')
                for key in ('initial_state','spatial_stop','velocity_stop','region_radius'):
                    compare_tree(r[key],prior[key])
                with np.load(path,allow_pickle=False) as data:a=dict(data)
                verify_record_arrays(a,r,p)
    compare_tree(design['jobs'],expected)
    require(len(seeds)==len(set(seeds))==72,'seed reuse')
    names={job['name'] for job in expected}
    for suffix in ('npz','json'):
        require({f.stem for f in OUT.glob('*-block*.'+suffix)}==names,'unexpected/missing job inventory')
    result=analyze();compare_tree(read(OUT/'summary.json'),result)
    require(len(result['comparisons'])==p['reporting']['efficiency_family_size']==6,'wrong ratio family')
    print('Verified 36 precision ensembles, 36864 quartets, identical reference inputs and design-aware statistics; no trajectories rerun.')
    return result


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args();verify()
