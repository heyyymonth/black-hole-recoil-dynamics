"""Audit BBH abundance in the available late-time CMC summaries."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fetch_e005 import fetch
from recoil_dynamics.catalog_counts import load_param, summarize_models

SOURCES = ['experiments/e005.json', 'experiments/e005-inputs.json', 'docs/e005-protocol.md',
           'experiments/fetch_e005.py', 'experiments/run_e005.py', 'experiments/verify_e005.py',
           'src/recoil_dynamics/catalog_counts.py', 'tests/test_catalog_counts.py',
           'tests/test_e005_fetch.py', 'results/e004/summary.json', 'results/e004/population.csv',
           'requirements-lock.txt']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distribution(values):
    return {'minimum': min(values), 'median': statistics.median(values), 'maximum': max(values)} if values else None


def calculate(config):
    paths = fetch(config)
    rows = [load_param(local, model=name.split('/')[1], source_path=name)
            for name, local in sorted(paths.items())]
    models = summarize_models(rows, reference_age_gyr=config['comparison_age_gyr'],
                              tolerance_gyr=config['age_tolerance_gyr'],
                              minimum_bbh_count=config['minimum_bbh_count_screen'])
    by_path = {row['path']: row for row in rows}
    selected = [dict(by_path[model['selected_path']],
                     age_offset_gyr=by_path[model['selected_path']]['age_gyr']-config['comparison_age_gyr'])
                for model in models if model['selected_path'] is not None]
    missing = [model['model'] for model in models if model['selected_path'] is None]
    ref = by_path[config['e004_reference_path']]
    old = json.loads((ROOT/'results/e004/summary.json').read_text())
    with (ROOT/'results/e004/population.csv').open() as stream:
        old_counts = {row['population']: int(row['count']) for row in csv.DictReader(stream)}
    comparisons = {
        'age_gyr': (ref['age_gyr'], old['scales']['age_gyr'], 1e-12),
        'host_mass_msun': (ref['host_mass_msun'], old['scales']['total_mass_msun'], 1e-7),
        'system_count': (ref['system_count'], old_counts['all_systems'], 0),
        'binary_count': (ref['binary_count'], old_counts['all_binaries'], 0),
        'single_bh_count': (ref['single_bh_count'], old_counts['single_black_holes'], 0),
        'bh_bearing_binary_count': (ref['bh_nonbbh_binary_count']+ref['bbh_count'], old_counts['bh_bearing_binaries'], 0),
        'bbh_count': (ref['bbh_count'], old_counts['binary_black_holes'], 0),
    }
    checks = [{'quantity': name, 'catalog': a, 'full_snapshot': b, 'relative_tolerance': tol,
               'passed': math.isclose(a, b, rel_tol=tol, abs_tol=0)} for name, (a,b,tol) in comparisons.items()]
    if not all(c['passed'] for c in checks):
        raise ValueError('E004 cross-check failed')
    counts = [r['bbh_count'] for r in selected]
    screen = config['minimum_bbh_count_screen']
    summary = {
        'experiment': 'e005', 'snapshot_records': len(rows), 'model_histories': len(models),
        'available_age_gyr': distribution([r['age_gyr'] for r in rows]),
        'records_at_or_below_10_gyr': sum(r['age_gyr'] <= 10 for r in rows),
        'count_accounting_valid_records': len(rows),
        'comparison_age_gyr': config['comparison_age_gyr'], 'age_tolerance_gyr': config['age_tolerance_gyr'],
        'selected_models': len(selected), 'models_without_age_coverage': missing,
        'selected_age_gyr': distribution([r['age_gyr'] for r in selected]),
        'selected_bbh_counts': distribution(counts),
        'selected_total_bh_counts': distribution([r['bh_count'] for r in selected]),
        'selected_all_binary_counts': distribution([r['binary_count'] for r in selected]),
        'selected_zero_bbh_models': sum(c == 0 for c in counts),
        'selected_at_most_two_bbh_models': sum(c <= 2 for c in counts),
        'selected_below_count_screen_models': sum(c < screen for c in counts),
        'selected_meeting_count_screen_models': sum(c >= screen for c in counts),
        'minimum_bbh_count_screen': screen,
        'saved_state_bbh_count_range': [min(r['bbh_count'] for r in rows), max(r['bbh_count'] for r in rows)],
        'histories_with_any_saved_state_meeting_screen': len({r['model'] for r in rows if r['bbh_count'] >= screen}),
        'e004_cross_checks': checks,
        'e004_model_history': next(m for m in models if m['model'] == ref['model']),
        'joint_companion_distribution_calibrated': False,
        'live_recoil_trajectories_validated': False,
        'merger_conditioned_population_measured': False,
        'interpretation': 'Equal-weight available model states near 12 Gyr; saved-state history ranges. '
                          'No cosmological weights, time averages, distinct-object census, encounter rates, or live recoil validation.'}
    return {'snapshots': rows, 'models': models, 'selected': selected}, summary, paths


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def plot(tables, summary, out):
    plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 10,
                         'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5), layout='constrained')
    values = [r['bbh_count'] for r in tables['selected']]
    maximum = max(values, default=0)
    bins = range(0, max(32, maximum+2))
    axes[0].hist(values, bins=[b-.5 for b in bins], color='#14677d', edgecolor='white')
    axes[0].axvline(summary['minimum_bbh_count_screen']-.5, color='#be6d22', ls='--', label='Count-30 screen')
    axes[0].set(xlabel='BH–BH binaries in selected state', ylabel='Number of model histories',
                title=f"(a) One state near 12 Gyr per model (N={len(values)})")
    axes[0].legend(frameon=False, fontsize=8)
    for row in tables['selected']:
        axes[1].scatter(row['bh_count'], row['bbh_count'], s=16, color='#14677d', alpha=.65)
    axes[1].axhline(summary['minimum_bbh_count_screen'], color='#be6d22', ls='--')
    axes[1].set(xlabel='All black-hole objects in selected state', ylabel='BH–BH binary systems',
                title='(b) All BHs are not available BH–BH targets', ylim=(-1, max(32, maximum+2)))
    for ax in axes:
        ax.grid(axis='y', color='#dddddd', linewidth=.4)
        ax.set_axisbelow(True)
    fig.savefig(out/'target_counts.png', dpi=200)
    fig.savefig(out/'target_counts.pdf', metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)


def main():
    dirty = subprocess.check_output(['git', 'status', '--porcelain', '--', *SOURCES], cwd=ROOT, text=True)
    if dirty.strip():
        raise SystemExit('Commit E005 source inputs before generating recorded outcomes.')
    config = json.loads((ROOT/'experiments/e005.json').read_text())
    tables, summary, paths = calculate(config)
    out = ROOT/'results/e005'
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        write_csv(out/(name+'.csv'), rows)
    (out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    plot(tables, summary, out)
    manifest = {'experiment': 'e005', 'researcher': 'Venkata Hemanth Athota',
                'run_utc': datetime.now(timezone.utc).isoformat(), 'python': platform.python_version(),
                'versions': {'matplotlib': importlib.metadata.version('matplotlib')},
                'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'protocol_revision': '72916d1', 'external_source': config['source_repository'],
                'external_revision': config['source_revision'],
                'external_sha256': {name: sha(path) for name, path in sorted(paths.items())},
                'source_sha256': {name: sha(ROOT/name) for name in SOURCES},
                'output_sha256': {str(path.relative_to(ROOT)): sha(path) for path in sorted(out.iterdir())
                                  if path.suffix in ['.csv', '.json', '.png', '.pdf']},
                'rows': {name: len(rows) for name, rows in tables.items()}}
    (ROOT/'results/run_e005.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
