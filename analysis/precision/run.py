"""Run the fixed covariance-coupling and insertion-stratification pilot.

The earlier experiment is immutable. This driver reuses its complete reward,
transition, stopping and exact sampler code through an isolated transport module.
Every method estimates the same finite-step response. No validation is rerun.
"""
from pathlib import Path
import argparse
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.response import run as reference
from analysis.response.disk_cache import DiskMatchedKernel
from analysis.precision.coupling import make_transport

OUT = ROOT/'results/precision'
PROTOCOL = ROOT/'analysis/precision/protocol.json'
CODE = ('analysis/precision/run.py', 'analysis/precision/coupling.py',
        'analysis/precision/statistics.py', 'analysis/precision/analyze.py',
        'analysis/response/disk_cache.py')
read = reference.read_json
digest = reference.digest


def jobs(p):
    for block in range(p['blocks']):
        for offset in range(len(p['launches'])):
            li = (offset+block) % len(p['launches'])
            for mo in range(len(p['methods'])):
                mi = (mo+block+li) % len(p['methods'])
                method = p['methods'][mi]
                yield {'name': f"{p['launches'][li]}-{method['name']}-block{block}",
                       'launch': p['launches'][li], 'launch_index': li,
                       'method': method, 'method_index': mi, 'block': block,
                       'paths': p['quartets_per_block'],
                       'seed': [p['seed'], 40, li, mi, block],
                       'selection_seed': [p['seed'], 41, li, mi, block]}


def selections(p, job):
    if job['method']['selection'] == 'uniform':
        return None, np.full(job['paths'], -1, dtype=np.int64)
    edges = np.asarray(p['step_strata'], dtype=np.int64)
    nstrata = len(edges)-1
    if job['paths'] % nstrata:
        raise ValueError('Equal allocation requires an integer count per stratum')
    rng = np.random.default_rng(np.random.SeedSequence(job['selection_seed']))
    n = job['paths']//nstrata
    selected = np.concatenate([rng.integers(a, b, size=n) for a, b in zip(edges[:-1], edges[1:])])
    return selected, np.repeat(np.arange(nstrata, dtype=np.int64), n)


def identity(p):
    rp = read(ROOT/p['response_protocol'])
    rid = reference.identity(rp)
    saved = read(ROOT/p['reference_design'])
    if saved['identity'] != rid:
        raise ValueError('The frozen response inputs or implementation changed')
    names = [str(PROTOCOL.relative_to(ROOT)), p['reference_design'],
             'results/response/summary.json', 'results/precision/preflight.txt', *CODE]
    return {'reference_identity': rid,
            'source_sha256': {name: digest(ROOT/name) for name in names}}


def ensure_design(p, ids):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT/'execution-design.json'
    core = {'revision': p['revision'], 'identity': ids, 'jobs': list(jobs(p)),
            'scope': p['scope']}
    if path.exists():
        old = read(path)
        if any(old.get(k) != v for k, v in core.items()):
            raise ValueError('Changed fixed design; do not reuse these records')
    else:
        if list(OUT.glob('*-block*.json')) or list(OUT.glob('*-block*.npz')):
            raise ValueError('Outcomes predate the execution design')
        reference.write_new_json(path, {**core, 'created_utc': reference.utc_now()})
    return digest(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    p = read(PROTOCOL); ids = identity(p); design = ensure_design(p, ids)
    rp = read(ROOT/p['response_protocol']); matched = read(ROOT/rp['matched_protocol'])
    table = ROOT/rp['kernel_tables']['fine']; metadata = read(table.with_suffix('.json'))
    if not (metadata['completed'] and metadata['all_vertex_quadrature_passed']
            and metadata['all_holdout_quadrature_passed']):
        raise ValueError('Accepted fine encounter table required')
    cache = p['cache']
    kernel = DiskMatchedKernel(ROOT/rp['host_path'], table, ROOT/rp['matched_protocol'],
        disk_cache_dir=ROOT/cache['directory'], namespace=digest(ROOT/cache['reuse_namespace']),
        grid='fine', cache_bytes=cache['retained_array_bytes'], cache_vertices=cache['vertex_cap'],
        disk_max_bytes=cache['disk_max_bytes'])
    force, density = reference.host_functions(kernel)
    inherited = read(ROOT/'results/collisionless-residence/NGC5139.json')
    force_domain = min(row['summary']['stop_radius_pc'] for row in inherited['records'])/kernel.r_scale_pc
    controls = {'spatial_stop': float(min(kernel.axes[0][-1], force_domain)),
                'velocity_stop': float(kernel.axes[1][-1]), 'density': density,
                'region_radius': matched['region_radius_over_density_radius']}
    modules = {factor: make_transport(factor) for factor in ('eigen', 'symmetric')}
    for job in jobs(p):
        path = OUT/(job['name']+'.npz'); meta = path.with_suffix('.json')
        if path.exists() or meta.exists():
            if not args.resume or not (path.exists() and meta.exists()):
                raise ValueError('Existing/incomplete outcome: '+job['name'])
            r = read(meta)
            if (not r['complete'] or r['identity'] != ids or r['job'] != job
                    or r['design_sha256'] != design or r['samples_sha256'] != digest(path)):
                raise ValueError('Changed retained outcome: '+job['name'])
            print('verified resume', job['name'], flush=True)
            continue
        if identity(p) != ids:
            raise ValueError('Sources changed before execution')
        launch = next(v for v in rp['launches'] if v['name'] == job['launch'])
        initial = reference.initial_state(kernel, matched, launch)
        selected, strata = selections(p, job)
        started_utc = reference.utc_now(); before = kernel.disk_statistics()
        wall = time.perf_counter(); cpu = time.process_time()
        print('starting', job['name'], flush=True)
        result = modules[job['method']['factor']].run_weighted_response(kernel, force,
            np.repeat(initial[None, :], job['paths'], axis=0), p['horizon'], p['dt'], job['seed'],
            replacement_steps=selected, **controls)
        cpu = time.process_time()-cpu; wall = time.perf_counter()-wall
        arrays, branches = reference.flatten_result(result, 'weighted', p['horizon'], p['dt'])
        arrays['insertion_stratum'] = strata
        if identity(p) != ids:
            raise ValueError('Sources changed during execution')
        after = kernel.disk_statistics()
        after['cache_path'] = str(Path(after['cache_path']).relative_to(ROOT))
        record = {'complete': True, 'revision': p['revision'], 'identity': ids,
            'design_sha256': design, 'job': job, 'started_utc': started_utc,
            'completed_utc': reference.utc_now(), 'cpu_seconds': cpu, 'wall_seconds': wall,
            'horizon': p['horizon'], 'dt': p['dt'], 'n_steps': result['n_steps'],
            'initial_state': initial.tolist(), 'spatial_stop': controls['spatial_stop'],
            'velocity_stop': controls['velocity_stop'], 'region_radius': controls['region_radius'],
            'branches': list(branches), 'summary': reference.sample_summary(arrays['response_sample']),
            'summary_scope': 'Descriptive raw-score statistics; design-aware uncertainty is computed separately',
            'numerical_missing_bound': reference.sample_summary(arrays['numerical_missing_bound']),
            'stops': {b: dict(zip(*[v.tolist() for v in np.unique(result[b]['stop_cause'], return_counts=True)])) for b in branches},
            'cache': after, 'cache_before': {k: v for k, v in before.items() if k != 'cache_path'},
            'array_inventory': {k: {'shape': list(v.shape), 'dtype': str(v.dtype)} for k, v in arrays.items()}}
        temporary = path.with_suffix('.npz.tmp')
        with temporary.open('wb') as f:
            np.savez_compressed(f, **arrays)
        temporary.replace(path); record['samples_sha256'] = digest(path)
        reference.write_new_json(meta, record)
        print('finished', job['name'], 'cpu', round(cpu, 2), 'wall', round(wall, 2),
              'mean', record['summary']['mean'], flush=True)
    print('completed precision pilot', flush=True)


if __name__ == '__main__':
    main()
