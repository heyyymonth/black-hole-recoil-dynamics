"""Summarize retained response and independent validation samples; NumPy only.

No dynamics run here. Intervals concern the fixed, numerically stopped process.
The continuation bounds describe the retained finite sample, not an upper
confidence bound on missing residence in the population.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.response.statistics import summarize_response, prediction_error, _interval

OUT = ROOT/'results/response'
PROTOCOL = ROOT/'analysis/response/protocol.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def result_name(stage, setting, launch, method='weighted', a=None):
    suffix = f'a{a:.3f}' if stage == 'validation' else method
    return f'{stage}-{setting}-{launch}-{suffix}'


def classification(interval, tolerance):
    low, high = interval
    if low >= -tolerance and high <= tolerance:
        return 'within_tolerance'
    if low > tolerance or high < -tolerance:
        return 'outside_tolerance'
    return 'inconclusive'


def contribution_bins(a, score, edges, coordinate):
    active = a['active_at_replacement']
    rows = []
    covered = np.zeros(len(score), dtype=bool)
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = active & (coordinate >= lo) & (coordinate < hi)
        if index == len(edges)-2:
            mask = active & (coordinate >= lo) & (coordinate <= hi)
        covered |= mask
        rows.append({'low': lo, 'high': hi, 'count': int(mask.sum()),
                     'contribution': float(score[mask].sum()/len(score))})
    return {'bins': rows, 'denominator': len(score),
            'outside_count': int((active & ~covered).sum()),
            'outside_contribution': float(score[active & ~covered].sum()/len(score)),
            'inactive_count': int((~active).sum()),
            'inactive_contribution': float(score[~active].sum()/len(score)),
            'scope': 'Descriptive signed decomposition, not separate confidence statements'}


def analyze():
    p = json.loads(PROTOCOL.read_text())
    confidence = p['reporting']['confidence']
    family = p['reporting']['bonferroni_family_size']
    tol = p['reporting']['prediction_absolute_error_target']
    inputs = {str(PROTOCOL.relative_to(ROOT)): digest(PROTOCOL)}
    for name in ['execution-design.json', 'sample-count-lock.json', 'prediction-register.json']:
        path = OUT/name
        inputs[str(path.relative_to(ROOT))] = digest(path)
    adapter_path = OUT/'storage-adapter.json'
    if adapter_path.exists():
        adapter = json.loads(adapter_path.read_text())
        inputs[str(adapter_path.relative_to(ROOT))] = digest(adapter_path)
        for name in adapter['source_sha256']:
            inputs[name] = digest(ROOT/name)
        evidence = adapter['test_evidence']
        inputs[evidence['path']] = digest(ROOT/evidence['path'])
    count_lock = json.loads((OUT/'sample-count-lock.json').read_text())
    primary = next(s for s in p['numerical_settings'] if s['role'] == 'primary')
    samples, records = {}, {}

    def load(stage, setting, launch, method='weighted', a=None):
        name = result_name(stage, setting['name'], launch['name'], method, a)
        path = OUT/(name+'.npz'); meta = path.with_suffix('.json')
        record = json.loads(meta.read_text())
        if digest(path) != record['samples_sha256']:
            raise ValueError('Changed retained samples: '+name)
        for f in (path, meta): inputs[str(f.relative_to(ROOT))] = digest(f)
        arrays = dict(np.load(path, allow_pickle=False))
        records[name], samples[name] = record, arrays
        return name, record, arrays

    pilot = []
    for launch in p['launches']:
        methods = {}
        for method in p['pilot']['methods']:
            name, r, a = load('pilot', primary, launch, method)
            v = a['response_sample']
            methods[method] = {'record': name, 'n': len(v), 'mean': float(v.mean()),
                'variance': float(v.var(ddof=1)), 'wall_seconds': r['wall_seconds'],
                'variance_times_seconds_per_pair': float(v.var(ddof=1)*r['wall_seconds']/len(v))}
        pilot.append({'launch': launch['name'], 'methods': methods,
            'variance_ratio_naive_over_weighted': methods['naive']['variance']/methods['weighted']['variance'] if methods['weighted']['variance'] else None,
            'cost_adjusted_ratio_naive_over_weighted': methods['naive']['variance_times_seconds_per_pair']/methods['weighted']['variance_times_seconds_per_pair'] if methods['weighted']['variance_times_seconds_per_pair'] else None,
            'count_selection': count_lock['selection'][launch['name']],
            'scope': 'Independent exploratory variance estimates; pilot excluded from inference; cache warmup affects measured costs'})

    predictions = []
    for setting in p['numerical_settings']:
        for launch in p['launches']:
            name, r, a = load('prediction', setting, launch)
            score = a['response_sample']; bound = a['numerical_missing_bound']
            stats = summarize_response(score, confidence, family)
            absolute = np.sort(abs(score))[::-1]
            top = max(1, int(np.ceil(.01*len(score))))
            predictions.append({'name': name, 'setting': setting['name'], 'launch': launch['name'],
                'grid': setting['grid'], 'dt': setting['dt'], 'role': setting['role'],
                'response': stats,
                'precision_target_met_by_observed_interval': bool((stats['interval'][1]-stats['interval'][0])/2 <= p['prediction']['primary_halfwidth_target']),
                'top_one_percent_absolute_contribution_fraction': float(absolute[:top].sum()/absolute.sum()) if absolute.sum() else 0.,
                'numerical_continuation': {'sample_mean_bound': float(bound.mean()),
                    'sample_maximum_bound': float(bound.max()), 'nonzero_rows': int(np.count_nonzero(bound)),
                    'budget': p['reporting']['numerical_missing_residence_budget'],
                    'scope': 'Finite-sample pathwise bound; not a population confidence bound'},
                'time_contributions': contribution_bins(a, score, p['prediction']['time_bin_edges'], a['replacement_steps']*setting['dt']),
                'radius_contributions': contribution_bins(a, score, p['prediction']['radius_bin_edges'], np.linalg.norm(a['replacement_states'][:, :3], axis=1)),
                'stops': r['stops']})

    validations = []
    for setting in p['numerical_settings']:
        values = p['validation']['primary_a_values'] if setting['role'] == 'primary' else p['validation']['control_a_values']
        for launch in p['launches']:
            pred_name = result_name('prediction', setting['name'], launch['name'])
            slope = samples[pred_name]['response_sample']
            for value in values:
                name, r, a = load('validation', setting, launch, a=value)
                error = prediction_error(a['mixture__normalized_residence'], a['baseline__normalized_residence'], slope, value, confidence, family, paired_validation=True)
                shift = summarize_response(a['paired_difference'], confidence, family)
                shift['interpretation'] = 'Paired mixture-minus-diffusion mean at declared a'
                shift['normalization'] = '(C_mixture-C_diffusion)/T'
                validations.append({'name': name, 'setting': setting['name'], 'launch': launch['name'], 'a': value,
                    'role': setting['role'], 'observed_shift': shift, 'prediction_error': error,
                    'scaled_remainder': {'estimate': error['estimate']/value,
                        'interval': [x/value for x in error['interval']],
                        'scope': 'Same residual event rescaled by a; finite-a remainder is not the exact derivative'},
                    'endpoint_classification': classification(error['interval'], tol) if value == 1. else None,
                    'numerical_continuation_sample_mean_bound': float(a['numerical_missing_bound'].mean()),
                    'stops': r['stops']})

    refinements = []
    coarse, fine, primary = p['numerical_settings']
    for launch in p['launches']:
        for kind, first, second in [('grid', fine, coarse), ('time', primary, fine)]:
            for observable in ['response', 'endpoint_shift']:
                stage = 'prediction' if observable == 'response' else 'validation'
                variant = None if observable == 'response' else 1.
                names = [result_name(stage, s['name'], launch['name'], a=variant) for s in [first, second]]
                key = 'response_sample' if observable == 'response' else 'paired_difference'
                x, y = [samples[n][key] for n in names]
                estimate = float(x.mean()-y.mean())
                variance = float(x.var(ddof=1)/len(x)+y.var(ddof=1)/len(y))
                stats = _interval(estimate, variance, confidence, family)
                refinements.append({'kind': kind, 'launch': launch['name'], 'observable': observable,
                    'first': names[0], 'second': names[1], 'difference': stats,
                    'classification': classification(stats['interval'], p['reporting']['numerical_sensitivity_budget']),
                    'scope': 'Independent finite-resolution sensitivity contrast; not a truncation bound'})
    for f in ['analysis/response/analyze.py', 'analysis/response/statistics.py', 'analysis/matched/statistics.py']:
        inputs[f] = digest(ROOT/f)
    diagnostic_path = OUT/'coupling-diagnostic.json'
    diagnostic = json.loads(diagnostic_path.read_text())
    inputs[str(diagnostic_path.relative_to(ROOT))] = digest(diagnostic_path)
    inputs.update(diagnostic['input_and_source_sha256'])
    return {'revision': p['revision'], 'cluster': p['host'],
        'cluster_residence_response_evaluated': True, 'discrete_response_validation_performed': True,
        'continuum_accuracy_classification': 'unresolved',
        'scope': p['scope'], 'reporting': p['reporting'],
        'inferential_comparisons': 63, 'bonferroni_family_size': family,
        'counts': {'pilot_ensembles': len(pilot)*2, 'prediction_ensembles': len(predictions),
            'validation_ensembles': len(validations),
            'independent_groups': sum(r['job']['paths'] for r in records.values()),
            'branch_trajectories': sum(r['job']['paths']*len(r['branches']) for r in records.values())},
        'pilot': pilot, 'predictions': predictions, 'validations': validations, 'refinements': refinements,
        'coupling_diagnostic': {'revision': diagnostic['revision'], 'scope': diagnostic['scope'],
            'interpretation': diagnostic['interpretation'],
            'counts': {'synthetic': len(diagnostic['synthetic']),
                       'local_pairs': sum(len(row['probes']) for row in diagnostic['local'])}},
        'input_sha256': inputs,
        'interpretation': 'Fixed numerical-process test of a mean-response extrapolation. Normal intervals are approximate; rare influential scores, sampling, grid/time sensitivity and numerical continuation remain separate. No physical-cluster or merger-rate calibration.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true'); args = parser.parse_args()
    result = analyze(); path = OUT/'summary.json'
    if args.verify:
        from analysis.matched.verify import compare_tree
        compare_tree(json.loads(path.read_text()), result)
        print('Verified response summaries from retained samples; no trajectories rerun.')
    else:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        print('Summarized', result['counts'], 'conditional response experiment.')


if __name__ == '__main__':
    main()
