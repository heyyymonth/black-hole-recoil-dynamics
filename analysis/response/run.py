"""Execute the prospectively registered, staged residence-response test.

Pilot rows choose sample sizes from variance only and never enter production
estimates. All production predictions are locked before validation starts.
Resuming requires exact source/input hashes and intact retained arrays. This
driver does not reinterpret an incomplete precision test as a positive result.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import gc
import json
import math
import sys
import time

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.matched.kernel import digest
from analysis.matched.run import host_functions
from analysis.matched.statistics import _normal_critical
from analysis.response.cache import CompactMatchedKernel
from analysis.response.transport import (
    run_paired_mixture, run_paired_response, run_weighted_response,
)

PROTOCOL = ROOT/'analysis/response/protocol.json'
OUT = ROOT/'results/response'
CODE = (
    'analysis/response/run.py', 'analysis/response/transport.py',
    'analysis/response/cache.py', 'analysis/matched/run.py',
    'analysis/matched/transport.py', 'analysis/matched/kernel.py',
    'analysis/matched/statistics.py',
    'analysis/host_reconstruction/host.py', 'analysis/encounters/physics.py',
    'analysis/encounters/run.py',
)
ARRAYS = ('states', 'residence', 'normalized_residence', 'exposure',
          'normalized_exposure', 'end_time', 'stop_cause', 'alive',
          'elastic_events', 'evaluated_steps')
INSERTION_ARRAYS = ('replacement_steps', 'active_at_replacement', 'replacement_states')
WEIGHT_ARRAYS = ('jump_zero_weight', 'jump_positive_weight',
                'selected_poisson_mean', 'positive_branch_count')
STAGE_CODES = {'pilot': 10, 'prediction': 20, 'validation': 30}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_new_json(path, record):
    path = Path(path)
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite {path}')
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def identity(p):
    matched = read_json(ROOT/p['matched_protocol'])
    inputs = [str(PROTOCOL.relative_to(ROOT)), p['matched_protocol'], p['host_path'],
              matched['encounter_protocol'], matched['baseline_protocol'],
              'results/collisionless-residence/NGC5139.json']
    for table in p['kernel_tables'].values():
        inputs.extend([table, str(Path(table).with_suffix('.json'))])
    return {'protocol_sha256': digest(PROTOCOL),
            'input_sha256': {name: digest(ROOT/name) for name in inputs},
            'implementation_sha256': {name: digest(ROOT/name) for name in CODE},
            'numpy_version': np.__version__, 'scipy_version': scipy.__version__}


def job_name(stage, setting, launch, method=None, a=None):
    tail = f'a{a:.3f}' if stage == 'validation' else method
    return f'{stage}-{setting["name"]}-{launch["name"]}-{tail}'


def jobs(p, stage, counts=None):
    """Canonical inventory and stable seed coordinates, independent of run order."""
    for si, setting in enumerate(p['numerical_settings']):
        if stage == 'pilot' and setting['name'] != p['pilot']['setting']:
            continue
        for li, launch in enumerate(p['launches']):
            if stage == 'pilot':
                variants = [(method, None, p['pilot']['pairs_per_method_per_launch'])
                            for method in p['pilot']['methods']]
            elif stage == 'prediction':
                if counts is None:
                    raise ValueError('Prediction inventory requires locked sample counts')
                n = (counts[launch['name']] if setting['role'] == 'primary'
                     else p['prediction']['control_pairs'])
                variants = [('weighted', None, n)]
            else:
                values = (p['validation']['primary_a_values'] if setting['role'] == 'primary'
                          else p['validation']['control_a_values'])
                variants = [('paired_mixture', a, p['validation']['paired_paths_per_setting'])
                            for a in values]
            for vi, (method, a, n) in enumerate(variants):
                # Prediction has one registered estimator, hence method index 0.
                seed = [p['seed'], STAGE_CODES[stage], si, li, vi]
                yield {'name': job_name(stage, setting, launch, method, a),
                       'stage': stage, 'setting_index': si, 'launch_index': li,
                       'variant_index': vi, 'setting': setting, 'launch': launch,
                       'method': method, 'a': a, 'paths': int(n), 'seed': seed}


def initial_state(kernel, matched, launch):
    radius = matched['launch']['radius_over_density_radius']
    apocenter = matched['launch']['radial_apocenter_over_density_radius']
    speed = np.sqrt(2*(kernel.host.psi(radius*kernel.r_scale_pc)
                      - kernel.host.psi(apocenter*kernel.r_scale_pc)))/kernel.v_scale_kms
    angle = launch['angle_from_radial_degrees']
    if angle == 0:
        direction = np.array([1., 0., 0.])
    elif angle == 90:
        direction = np.array([0., 1., 0.])
    else:
        radians = np.deg2rad(angle)
        direction = np.array([np.cos(radians), np.sin(radians), 0.])
    return np.r_[radius, 0., 0., speed*direction]


def missing_residence(branch, horizon):
    numerical = np.isin(branch['stop_cause'], ['spatial_domain', 'velocity_domain'])
    return np.where(numerical, np.maximum(0., 1-branch['end_time']/horizon), 0.)


def flatten_result(result, method, horizon, dt):
    """Retain every branch; auxiliary contrasts remain independently recomputable."""
    if method == 'weighted':
        branches = ('jump_zero', 'jump_positive', 'diffusion_plus', 'diffusion_minus')
    elif method == 'naive':
        branches = ('baseline', 'replacement')
    elif method == 'paired_mixture':
        branches = ('baseline', 'mixture')
    else:
        raise ValueError(f'Unknown method {method}')
    arrays = {f'{branch}__{key}': result[branch][key]
              for branch in branches for key in ARRAYS}
    if method != 'paired_mixture':
        arrays.update({key: result[key] for key in INSERTION_ARRAYS})
    steps = int(np.ceil(horizon/dt))
    if method == 'weighted':
        arrays.update({key: result[key] for key in WEIGHT_ARRAYS})
        p0, pp = result['jump_zero_weight'], result['jump_positive_weight']
        reward = lambda key: result[key]['normalized_residence']
        arrays['response_sample'] = steps*(p0*reward('jump_zero')+pp*reward('jump_positive')
                - .5*(reward('diffusion_plus')+reward('diffusion_minus')))
        bound = {key: missing_residence(result[key], horizon) for key in branches}
        arrays['numerical_missing_bound'] = steps*(p0*bound['jump_zero']
                + pp*bound['jump_positive']+.5*(bound['diffusion_plus']+bound['diffusion_minus']))
    elif method == 'naive':
        arrays['response_sample'] = steps*(result['replacement']['normalized_residence']
                                          - result['baseline']['normalized_residence'])
        arrays['numerical_missing_bound'] = steps*sum(
            missing_residence(result[key], horizon) for key in branches)
    else:
        arrays['paired_difference'] = (result['mixture']['normalized_residence']
                                      - result['baseline']['normalized_residence'])
        arrays['numerical_missing_bound'] = sum(
            missing_residence(result[key], horizon) for key in branches)
    return arrays, branches


def sample_summary(values):
    values = np.asarray(values)
    return {'count': int(len(values)), 'mean': float(values.mean()),
            'variance': float(values.var(ddof=1)),
            'standard_error': float(values.std(ddof=1)/np.sqrt(len(values))),
            'minimum': float(values.min()), 'maximum': float(values.max())}


def checked_record(job, ids, dependencies):
    path = OUT/(job['name']+'.npz')
    record_path = path.with_suffix('.json')
    if not path.exists() or not record_path.exists():
        raise FileNotFoundError(f'Complete retained result required: {job["name"]}')
    record = read_json(record_path)
    if not record.get('complete') or record.get('job') != job:
        raise ValueError(f'Wrong or incomplete job record: {job["name"]}')
    if record.get('identity') != ids or record.get('stage_dependencies') != dependencies:
        raise ValueError(f'Stale source/input/stage identity: {job["name"]}')
    if record['samples_sha256'] != digest(path):
        raise ValueError(f'Changed samples: {job["name"]}')
    return record


def retained_manifest(job, record):
    path = OUT/(job['name']+'.npz')
    return {'name': job['name'], 'arrays': str(path.relative_to(ROOT)),
            'arrays_sha256': digest(path),
            'record': str(path.with_suffix('.json').relative_to(ROOT)),
            'record_sha256': digest(path.with_suffix('.json')),
            'summary': record['summary']}


def ensure_design(p, ids):
    path = OUT/'execution-design.json'
    core = {'revision': p['revision'], 'identity': ids,
            'seed_scheme': p['seed_scheme'],
            'stage_order': ['pilot', 'sample-count-lock', 'prediction', 'prediction-register', 'validation']}
    if path.exists():
        existing = read_json(path)
        if any(existing.get(key) != value for key, value in core.items()):
            raise ValueError('Execution design has stale sources or inputs; do not reuse outcomes')
    else:
        if list(OUT.glob('pilot-*.npz')) or list(OUT.glob('prediction-*.npz')) or list(OUT.glob('validation-*.npz')):
            raise ValueError('Research outcomes exist without an execution design')
        write_new_json(path, {**core, 'created_utc': utc_now()})
    return digest(path)


def lock_counts(p, ids, design_sha):
    """Only pilot weighted variance enters sample-size selection."""
    pilot = []
    selected = {}
    alpha = 1-p['reporting']['confidence']
    z = _normal_critical(1-alpha/p['reporting']['bonferroni_family_size'])
    minimum = p['prediction']['primary_pairs_min']
    maximum = p['prediction']['primary_pairs_max']
    target = p['prediction']['primary_halfwidth_target']
    for job in jobs(p, 'pilot'):
        record = checked_record(job, ids, {'execution_design_sha256': design_sha})
        pilot.append(retained_manifest(job, record))
        if job['method'] == 'weighted':
            with np.load(OUT/(job['name']+'.npz'), allow_pickle=False) as arrays:
                variance = float(arrays['response_sample'].var(ddof=1))
            needed = max(float(minimum), z*z*variance/(target*target))
            uncapped = 2**int(math.ceil(math.log2(needed)))
            selected[job['launch']['name']] = {
                'pilot_variance': variance, 'uncapped_variance_rule_count': uncapped,
                'selected_count': min(maximum, uncapped), 'cap_hit': uncapped > maximum}
    core = {'revision': p['revision'], 'identity': ids,
            'execution_design_sha256': design_sha,
            'pilot_inventory': pilot, 'selection': selected, 'normal_critical': z,
            'target_halfwidth': target, 'pilot_excluded_from_production': True,
            'selection_uses': 'Weighted pilot variance only; not its mean, sign, or validation outcomes'}
    path = OUT/'sample-count-lock.json'
    if path.exists():
        saved = read_json(path)
        if any(saved.get(key) != value for key, value in core.items()):
            raise ValueError('Sample-count lock no longer matches its registered pilot')
    else:
        if list(OUT.glob('prediction-*.npz')) or list(OUT.glob('validation-*.npz')):
            raise ValueError('Production outcomes exist before their sample-count lock')
        write_new_json(path, {**core, 'created_utc': utc_now()})
    return {key: value['selected_count'] for key, value in selected.items()}, digest(path)


def lock_predictions(p, ids, counts, design_sha, count_sha):
    dependencies = {'execution_design_sha256': design_sha, 'sample_count_lock_sha256': count_sha}
    inventory = [retained_manifest(job, checked_record(job, ids, dependencies))
                 for job in jobs(p, 'prediction', counts)]
    core = {'revision': p['revision'], 'identity': ids, **dependencies,
            'prediction_inventory': inventory,
            'claims': 'Frozen finite-step first-response predictions before independent validation',
            'known_endpoints': p['validation']['known_endpoints']}
    path = OUT/'prediction-register.json'
    if path.exists():
        saved = read_json(path)
        if any(saved.get(key) != value for key, value in core.items()):
            raise ValueError('Prediction register no longer matches retained predictions')
    else:
        if list(OUT.glob('validation-*.npz')):
            raise ValueError('Validation outcomes exist before a prediction register')
        write_new_json(path, {**core, 'created_utc': utc_now()})
    return digest(path)


class Runner:
    """Retain at most one kernel cache at once, reused across consecutive jobs."""
    def __init__(self, p, ids, resume):
        self.p, self.ids, self.resume = p, ids, resume
        self.matched = read_json(ROOT/p['matched_protocol'])
        self.kernel = None
        self.grid = None

    def get_kernel(self, grid):
        if self.grid != grid:
            self.kernel = None
            gc.collect()
            table = ROOT/self.p['kernel_tables'][grid]
            metadata = read_json(table.with_suffix('.json'))
            if not (metadata['completed'] and metadata['all_vertex_quadrature_passed']
                    and metadata['all_holdout_quadrature_passed']):
                raise ValueError('Accepted kernel quadrature and holdouts required')
            self.kernel = CompactMatchedKernel(ROOT/self.p['host_path'], table,
                ROOT/self.p['matched_protocol'],
                cache_bytes=self.p['cache']['retained_array_bytes'],
                cache_vertices=self.p['cache']['vertex_cap'])
            self.grid = grid
        return self.kernel

    def execute(self, job, dependencies):
        path = OUT/(job['name']+'.npz')
        if path.exists() or path.with_suffix('.json').exists():
            if not self.resume:
                raise FileExistsError(f'Result exists; use checked --resume: {job["name"]}')
            checked_record(job, self.ids, dependencies)
            print('verified resume', job['name'], flush=True)
            return
        if identity(self.p) != self.ids:
            raise ValueError('Sources or inputs changed before run')
        kernel = self.get_kernel(job['setting']['grid'])
        force, density = host_functions(kernel)
        inherited = read_json(ROOT/'results/collisionless-residence/NGC5139.json')
        force_domain = min(row['summary']['stop_radius_pc'] for row in inherited['records'])/kernel.r_scale_pc
        spatial = float(min(kernel.axes[0][-1], force_domain))
        velocity = float(kernel.axes[1][-1])
        initial = initial_state(kernel, self.matched, job['launch'])
        horizon, dt = self.p['horizon'], job['setting']['dt']
        controls = {'spatial_stop': spatial, 'velocity_stop': velocity, 'density': density,
                    'region_radius': self.matched['region_radius_over_density_radius']}
        args = (kernel, force, np.repeat(initial[None, :], job['paths'], axis=0),
                horizon, dt, job['seed'])
        started_utc = utc_now()
        started = time.monotonic()
        print('starting', job['name'], 'paths', job['paths'], flush=True)
        if job['method'] == 'weighted':
            result = run_weighted_response(*args, **controls)
        elif job['method'] == 'naive':
            result = run_paired_response(*args, **controls)
        else:
            result = run_paired_mixture(*args, mixture_probability=job['a'], **controls)
        elapsed = time.monotonic()-started
        arrays, branches = flatten_result(result, job['method'], horizon, dt)
        if identity(self.p) != self.ids:
            raise ValueError('Sources or inputs changed during run; outputs not registered')
        score = 'paired_difference' if job['stage'] == 'validation' else 'response_sample'
        summary = sample_summary(arrays[score])
        record = {'complete': True, 'revision': self.p['revision'], 'identity': self.ids,
                  'job': job, 'stage_dependencies': dependencies,
                  'started_utc': started_utc, 'completed_utc': utc_now(), 'wall_seconds': elapsed,
                  'horizon': horizon, 'dt': dt, 'n_steps': int(np.ceil(horizon/dt)),
                  'initial_state': initial.tolist(), 'spatial_stop': spatial,
                  'velocity_stop': velocity, 'inherited_force_domain': float(force_domain),
                  'region_radius': controls['region_radius'],
                  'scales': {'density_radius_pc': kernel.r_scale_pc,
                             'velocity_kms': kernel.v_scale_kms, 'time_Myr': kernel.time_scale_Myr},
                  'units': 'r_rho, v_dyn, t_dyn; normalized rewards divided by full horizon',
                  'coupling': result['coupling'], 'branches': list(branches), 'summary': summary,
                  'numerical_missing_bound': sample_summary(arrays['numerical_missing_bound']),
                  'stops': {branch: dict(zip(*[values.tolist() for values in
                      np.unique(result[branch]['stop_cause'], return_counts=True)])) for branch in branches},
                  'cache': kernel.cache_statistics(),
                  'array_inventory': {key: {'shape': list(value.shape), 'dtype': str(value.dtype)}
                                      for key, value in arrays.items()}}
        temporary = path.with_suffix('.npz.tmp')
        with temporary.open('wb') as handle:
            np.savez_compressed(handle, **arrays)
        temporary.replace(path)
        record['samples_sha256'] = digest(path)
        write_new_json(path.with_suffix('.json'), record)
        print('finished', job['name'], 'seconds', round(elapsed, 2),
              'estimate', summary['mean'], 'se', summary['standard_error'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['pilot', 'prediction', 'validation', 'all'], required=True)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    p = read_json(PROTOCOL)
    OUT.mkdir(parents=True, exist_ok=True)
    ids = identity(p)
    design_sha = ensure_design(p, ids)
    runner = Runner(p, ids, args.resume)
    if args.stage in ('pilot', 'all'):
        for job in jobs(p, 'pilot'):
            runner.execute(job, {'execution_design_sha256': design_sha})
    counts, count_sha = lock_counts(p, ids, design_sha)
    if args.stage in ('prediction', 'all'):
        for job in jobs(p, 'prediction', counts):
            runner.execute(job, {'execution_design_sha256': design_sha,
                                 'sample_count_lock_sha256': count_sha})
    if args.stage in ('prediction', 'validation', 'all'):
        register_sha = lock_predictions(p, ids, counts, design_sha, count_sha)
    if args.stage in ('validation', 'all'):
        for job in jobs(p, 'validation'):
            runner.execute(job, {'execution_design_sha256': design_sha,
                                 'sample_count_lock_sha256': count_sha,
                                 'prediction_register_sha256': register_sha})
    print('completed stage', args.stage, flush=True)


if __name__ == '__main__':
    main()
