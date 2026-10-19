"""Run the predeclared E003 close-passage benchmark; retain baseline releases."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from recoil_dynamics.encounters import experiment_rows

SOURCE_FILES = [
    'experiments/e003.json', 'docs/e003-protocol.md',
    'experiments/run_e003.py', 'experiments/verify_e003.py',
    'src/recoil_dynamics/encounters.py', 'src/recoil_dynamics/plummer.py',
    'src/recoil_dynamics/__init__.py', 'tests/test_encounters.py',
    'pyproject.toml', 'requirements.txt', 'requirements-lock.txt',
]


def sha(relative):
    return hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def save_plot(rows, config, path):
    plt.rcParams.update({
        'font.family': 'STIXGeneral', 'font.size': 10, 'mathtext.fontset': 'stix',
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'svg.hashsalt': 'recoil-e003', 'axes.spines.top': False,
        'axes.spines.right': False, 'axes.linewidth': .7,
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), layout='constrained')
    positions = np.arange(len(rows))
    labels = [f"{row['eta']:g}" for row in rows]
    axes[0].bar(positions-.17, [r['density_ratio'] for r in rows], .32,
                facecolor='white', edgecolor='#444444', hatch='///', linewidth=.8,
                label='Density proxy')
    axes[0].bar(positions+.17, [r['encounter_ratio'] for r in rows], .32,
                color='#14677d', edgecolor='#14677d', linewidth=.8,
                label='Velocity-weighted passages')
    axes[0].set(xticks=positions, xticklabels=labels, ylim=(0, 1.09),
                xlabel=r'Kick fraction $\eta$ (selected cases)',
                ylabel='Fraction of stationary control', title='(a) Common time window')
    axes[0].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    axes[0].legend(loc='upper right', frameon=False, fontsize=8)
    errors = [r['proxy_relative_error'] for r in rows[1:]]
    positions = np.arange(len(errors))
    axes[1].bar(positions, errors, .55, color='#be6d22', edgecolor='#70441b', linewidth=.7)
    for pos, error in zip(positions, errors):
        small = abs(error) < config['proxy_relative_error_tolerance']
        axes[1].annotate(f'{error:+.1%}', (pos+.28 if small else pos, error),
                         xytext=(5, 0) if small else (0, 5 if error >= 0 else -12),
                         textcoords='offset points', ha='left' if small else 'center',
                         va='center' if small else 'baseline', fontsize=9)
    tolerance = config['proxy_relative_error_tolerance']
    axes[1].axhline(0, color='#444444', lw=.7)
    axes[1].axhline(tolerance, color='#444444', ls='--', lw=.9, label='10% screening tolerance')
    if min(errors) < 0:
        axes[1].axhline(-tolerance, color='#444444', ls='--', lw=.9)
    axes[1].set(xticks=positions, xticklabels=labels[1:],
                xlabel=r'Kick fraction $\eta$ (selected cases)',
                ylabel='(Density proxy / passage ratio) − 1', title='(b) Proxy error')
    axes[1].set_ylim(min(0, min(errors)*1.25), max(tolerance, max(errors))*1.28)
    axes[1].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    axes[1].legend(loc='upper left', frameon=False, fontsize=8)
    for ax in axes:
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#dddddd', linewidth=.45)
    fig.savefig(path.with_suffix('.png'), dpi=220)
    fig.savefig(path.with_suffix('.pdf'), metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(path.with_suffix('.svg'), metadata={'Date': None})
    svg = path.with_suffix('.svg')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)


def main():
    # Bind outputs to committed code, including any independent validation edits.
    if git('status', '--porcelain', '--', *SOURCE_FILES):
        raise SystemExit('Commit E003 source inputs before generating recorded results.')
    revision = git('rev-parse', 'HEAD')
    protocol_revision = git('log', '--reverse', '--format=%H', '--',
                            'docs/e003-protocol.md').splitlines()[0]
    config = json.loads((ROOT/'experiments/e003.json').read_text())
    scales, rows = experiment_rows(config)
    out = ROOT/'results/e003'
    out.mkdir(parents=True, exist_ok=True)
    with (out/'encounters.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        'scales': scales,
        'proxy_relative_error_tolerance': config['proxy_relative_error_tolerance'],
        'nonzero_kicks_exceeding_tolerance': [r['eta'] for r in rows if r['eta'] and r['exceeds_tolerance']],
        'maximum_absolute_proxy_relative_error': max(abs(r['proxy_relative_error']) for r in rows),
        'interpretation': 'Conditional close-passage optical depths in an imposed tracer bath; no acquisition or merger outcomes.',
    }
    (out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    save_plot(rows, config, out/'proxy_comparison')
    outputs = ['results/e003/'+name for name in (
        'encounters.csv', 'summary.json', 'proxy_comparison.png',
        'proxy_comparison.pdf', 'proxy_comparison.svg')]
    manifest = {
        'experiment': 'e003', 'researcher': 'Venkata Hemanth Athota',
        'run_utc': datetime.now(timezone.utc).isoformat(),
        'source_revision': revision, 'protocol_revision': protocol_revision,
        'python': platform.python_version(),
        'versions': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'matplotlib']},
        'source_sha256': {name: sha(name) for name in SOURCE_FILES},
        'output_sha256': {name: sha(name) for name in outputs},
        'rows': len(rows), 'stochastic_sampling': False,
    }
    (ROOT/'results/run_e003.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
