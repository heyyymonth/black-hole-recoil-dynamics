"""Analyze the fixed precision pilot from retained samples, without dynamics."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.precision.statistics import estimate, efficiency_bootstrap

OUT = ROOT/'results/precision'
PROTOCOL = ROOT/'analysis/precision/protocol.json'
BRANCHES = ('jump_zero', 'jump_positive', 'diffusion_plus', 'diffusion_minus')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def analyze():
    p = read(PROTOCOL); design = read(OUT/'execution-design.json')
    inputs = {str(path.relative_to(ROOT)): digest(path) for path in
              [PROTOCOL, OUT/'execution-design.json', Path(__file__), ROOT/'analysis/precision/statistics.py']}
    groups = {}; rows = []
    for launch in p['launches']:
        for method in p['methods']:
            key = (launch, method['name']); arrays = []; records = []; blocks = []
            for block in range(p['blocks']):
                name = f"{launch}-{method['name']}-block{block}"
                path = OUT/(name+'.npz'); meta = path.with_suffix('.json'); record = read(meta)
                if not record['complete'] or record['samples_sha256'] != digest(path):
                    raise ValueError('Incomplete or changed precision samples: '+name)
                inputs[str(path.relative_to(ROOT))] = digest(path)
                inputs[str(meta.relative_to(ROOT))] = digest(meta)
                with np.load(path, allow_pickle=False) as data: a = dict(data)
                strata = a['insertion_stratum'] if method['selection']=='stratified' else None
                blocks.append({'name': name, 'estimate': estimate(a['response_sample'], strata),
                               'cpu_seconds': record['cpu_seconds'], 'wall_seconds': record['wall_seconds']})
                arrays.append(a); records.append(record)
            a = {k: np.concatenate([v[k] for v in arrays]) for k in arrays[0]}
            strata = a['insertion_stratum'] if method['selection']=='stratified' else None
            stats = estimate(a['response_sample'], strata)
            cpu = sum(r['cpu_seconds'] for r in records); wall = sum(r['wall_seconds'] for r in records)
            covariance = np.cov(np.stack([a[b+'__normalized_residence'] for b in BRANCHES]), ddof=1)
            # This pooled covariance is descriptive, not a variance formula for
            # the weighted score (whose Poisson weights vary by path).
            weighted_jump = a['jump_zero_weight']*a['jump_zero__normalized_residence'] + a['jump_positive_weight']*a['jump_positive__normalized_residence']
            diffusion = .5*(a['diffusion_plus__normalized_residence']+a['diffusion_minus__normalized_residence'])
            stops = {}
            for b in BRANCHES:
                labels, counts = np.unique(a[b+'__stop_cause'], return_counts=True)
                stops[b] = dict(zip(labels.tolist(), counts.tolist()))
            time_contributions = []
            for lo, hi in zip(p['step_strata'][:-1], p['step_strata'][1:]):
                mask = (a['replacement_steps']>=lo)&(a['replacement_steps']<hi)
                time_contributions.append({'steps': [lo, hi], 'count': int(mask.sum()),
                    'signed_contribution': float(a['response_sample'][mask].sum()/len(a['response_sample']))})
            row = {'launch': launch, 'method': method['name'], 'estimate': stats, 'blocks': blocks,
                'cpu_seconds': cpu, 'wall_seconds': wall,
                'variance_times_cpu': float(stats['variance_of_mean']*cpu),
                'branch_covariance': covariance.tolist(), 'branch_order': list(BRANCHES),
                'weighted_jump_diffusion_covariance': float(np.cov(weighted_jump, diffusion, ddof=1)[0,1]),
                'weighted_jump_diffusion_correlation': (float(np.corrcoef(weighted_jump, diffusion)[0,1]) if np.std(weighted_jump)>0 and np.std(diffusion)>0 else None),
                'numerical_continuation_sample_mean_bound': float(a['numerical_missing_bound'].mean()),
                'numerical_continuation_nonzero_rows': int(np.count_nonzero(a['numerical_missing_bound'])),
                'stops': stops, 'time_contributions': time_contributions,
                'covariance_scope': 'Descriptive pooled branch statistics; no additional confidence claims'}
            rows.append(row); groups[key] = (a['response_sample'], strata, row)
    comparisons = []
    for li, launch in enumerate(p['launches']):
        for pi, (first, second) in enumerate(p['reporting']['efficiency_pairs']):
            a, sa, ra = groups[(launch, first)]; b, sb, rb = groups[(launch, second)]
            ratio = efficiency_bootstrap(a, sa, b, sb,
                reps=p['reporting']['bootstrap_replicates'],
                seed=[p['reporting']['bootstrap_seed'], li, pi],
                family_size=p['reporting']['efficiency_family_size'],
                cost_a=ra['cpu_seconds'], cost_b=rb['cpu_seconds'])
            mean_difference = ra['estimate']['mean']-rb['estimate']['mean']
            se = np.sqrt(ra['estimate']['variance_of_mean']+rb['estimate']['variance_of_mean'])
            comparisons.append({'launch': launch, 'first': first, 'second': second,
                'efficiency': ratio, 'mean_difference': float(mean_difference),
                'mean_difference_standard_error': float(se),
                'mean_comparison_scope': 'Descriptive independent-arm difference, not an equivalence test'})
    n = sum(r['estimate']['n'] for r in rows)
    return {'revision': p['revision'], 'scope': p['scope'], 'complete': True,
        'counts': {'ensembles': len(design['jobs']), 'quartets': n, 'branch_trajectories': 4*n},
        'reporting': p['reporting'], 'rows': rows, 'comparisons': comparisons,
        'input_sha256': inputs,
        'interpretation': 'Precision for the same finite-step response only. Earlier endpoint validation is not reused to certify prediction accuracy. Numerical continuation, refinement and physical calibration remain separate.'}


if __name__ == '__main__':
    result = analyze()
    (OUT/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Analyzed', result['counts'], 'without rerunning trajectories.')
