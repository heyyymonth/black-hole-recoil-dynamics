"""Audit a pinned CMC snapshot; compare frozen host monopoles at fixed kicks."""
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fetch_e004 import fetch
from recoil_dynamics.cluster_snapshot import (
    load_snapshot, population_audit, compress_shells, check_upstream_summary)
from recoil_dynamics.snapshot_potential import ShellPotential
from recoil_dynamics.plummer import PlummerOrbit

SOURCES = ['experiments/e004.json', 'docs/e004-protocol.md',
           'experiments/fetch_e004.py', 'experiments/run_e004.py', 'experiments/verify_e004.py',
           'src/recoil_dynamics/cluster_snapshot.py', 'src/recoil_dynamics/snapshot_potential.py',
           'src/recoil_dynamics/plummer.py', 'src/recoil_dynamics/__init__.py',
           'tests/test_cluster_snapshot.py', 'tests/test_snapshot_potential.py',
           'requirements-lock.txt']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calculate(config):
    paths = fetch(config)
    snapshot = load_snapshot(paths['snapshot'], paths['units'])
    upstream_checks = check_upstream_summary(snapshot, paths['reference'])
    scales, population = population_audit(snapshot, config['G_pc_kms2_per_msun'])
    b, mass = scales['matched_plummer_scale_pc'], scales['total_mass_msun']
    time = b*config['pc_in_km']/scales['velocity_unit_kms']/config['seconds_per_myr']
    scales['dynamical_time_myr'] = time
    trajectories = []
    for count in config['shell_resolutions']:
        r, m = compress_shells(snapshot.radius/b, snapshot.mass/mass, count)
        potential = ShellPotential(r, m)
        for eta in config['kick_fractions_of_matched_plummer_escape']:
            kick = math.sqrt(2)*eta
            bound = kick < potential.escape_speed
            apocenter = potential.apocenter(kick) if bound else None
            period = potential.first_return(kick) if bound else None
            plummer = PlummerOrbit(eta)
            trajectories.append({
                'shell_count': len(r), 'eta_plummer': eta,
                'kick_kms': kick*scales['velocity_unit_kms'],
                'eta_snapshot': kick/potential.escape_speed,
                'bound_in_isolated_snapshot': int(bound),
                'snapshot_apocenter_pc': apocenter*b if bound else None,
                'plummer_apocenter_pc': plummer.apocenter*b,
                'snapshot_return_myr': period*time if bound else None,
                'plummer_return_myr': plummer.period*time,
                'plummer_return_relative_error': plummer.period/period-1 if bound else None,
                'central_potential_over_GM_b': potential.central_potential,
            })
    finest = trajectories[-len(config['kick_fractions_of_matched_plummer_escape']):]
    previous = trajectories[-2*len(finest):-len(finest)]
    convergence = []
    for coarse, fine in zip(previous, finest):
        errors = {}
        for field in ['snapshot_apocenter_pc', 'snapshot_return_myr']:
            errors[field+'_relative_change'] = (abs(coarse[field]/fine[field]-1)
                                               if fine[field] is not None else None)
        converged = all(value is not None and value <= config['shell_convergence_relative_tolerance']
                        for value in errors.values())
        convergence.append({'eta_plummer': fine['eta_plummer'], **errors, 'converged': converged})

    # Adjacent compression levels need not bound error. Compare the final
    # approximation directly with every system shell in this same snapshot.
    full_potential = ShellPotential(snapshot.radius/b, snapshot.mass/mass)
    full_checks = []
    for fine in finest:
        kick = math.sqrt(2)*fine['eta_plummer']
        full_radius = full_potential.apocenter(kick)*b
        full_return = full_potential.first_return(kick)*time
        radius_error = fine['snapshot_apocenter_pc']/full_radius-1
        return_error = fine['snapshot_return_myr']/full_return-1
        full_checks.append({
            'eta_plummer': fine['eta_plummer'],
            'uncompressed_apocenter_pc': full_radius,
            'uncompressed_return_myr': full_return,
            'compressed_apocenter_relative_error': radius_error,
            'compressed_return_relative_error': return_error,
            'passed': max(abs(radius_error), abs(return_error)) <= config['shell_convergence_relative_tolerance']})

    # A figure grid describes the snapshot, not additional observations.
    x = np.geomspace(.001, 20, 160)  # radius / snapshot half-mass radius
    order = np.argsort(snapshot.radius)
    cmass = np.r_[0., np.cumsum(snapshot.mass[order])/mass]
    host_cdf = cmass[np.searchsorted(snapshot.radius[order], x*scales['half_mass_radius_pc'], side='right')]
    bbmask = snapshot.cohorts()['binary_black_holes']
    bbr = np.sort(snapshot.radius[bbmask]/scales['half_mass_radius_pc'])
    bb_cdf = np.searchsorted(bbr, x, side='right')/len(bbr) if len(bbr) else np.full(len(x), np.nan)
    y = x*scales['half_mass_radius_pc']/b
    profiles = [{'radius_over_rh': float(xx), 'host_mass_cdf': float(cc),
                 'plummer_mass_cdf': float(yy**3/(1+yy*yy)**1.5),
                 'bbh_number_cdf': float(bb) if np.isfinite(bb) else None,
                 'snapshot_potential_over_GM_b': potential.potential(float(yy)),
                 'plummer_potential_over_GM_b': -1/math.sqrt(1+yy*yy)}
                for xx, yy, cc, bb in zip(x, y, host_cdf, bb_cdf)]
    bbh = next(row for row in population if row['population'] == 'binary_black_holes')
    summary = {'scales': scales, 'upstream_summary_checks': upstream_checks,
               'binary_black_hole_count': bbh['count'],
               'minimum_bbh_count_for_distribution_calibration': config['minimum_bbh_count_for_distribution_calibration'],
               'companion_distribution_calibrated': False,
               'sample_count_screen_passed': bbh['count'] >= config['minimum_bbh_count_for_distribution_calibration'],
               'live_recoil_trajectories_validated': False,
               'snapshot_central_escape_kms': potential.escape_speed*scales['velocity_unit_kms'],
               'shell_convergence': convergence,
               'uncompressed_shell_count': len(full_potential.radii),
               'full_shell_reference': full_checks,
               'return_time_screen': config['return_time_relative_error_screen'],
               'kicks_exceeding_return_time_screen': [r['eta_plummer'] for r in finest
                   if r['plummer_return_relative_error'] is not None and
                   abs(r['plummer_return_relative_error']) > config['return_time_relative_error_screen']],
               'interpretation': 'One simulation snapshot. Population audit and frozen ballistic-potential comparison; no live kick or friction calibration.'}
    return {'population': population, 'trajectories': trajectories, 'profiles': profiles}, summary


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def plot(tables, summary, out):
    plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix',
                         'font.size': 10, 'pdf.fonttype': 42, 'svg.hashsalt': 'e004',
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1), layout='constrained')
    p = tables['profiles']
    x = [r['radius_over_rh'] for r in p]
    axes[0].semilogx(x, [r['host_mass_cdf'] for r in p], color='#333333', label='Snapshot host mass')
    axes[0].semilogx(x, [r['plummer_mass_cdf'] for r in p], color='#14677d', ls='--', label='Matched Plummer mass')
    axes[0].step(x, [r['bbh_number_cdf'] for r in p], where='post', color='#be6d22',
                 label=f"BH binaries (N={summary['binary_black_hole_count']})")
    axes[0].set(xlabel=r'Radius / host $r_h$', ylabel='Cumulative fraction',
                title='(a) Population concentration', ylim=(0, 1.03))
    axes[0].legend(fontsize=7, frameon=False, loc='lower right')
    axes[1].semilogx(x, [-r['snapshot_potential_over_GM_b'] for r in p], color='#333333', label='Frozen snapshot')
    axes[1].semilogx(x, [-r['plummer_potential_over_GM_b'] for r in p], color='#14677d', ls='--', label='Matched Plummer')
    axes[1].set(xlabel=r'Radius / host $r_h$', ylabel=r'Potential depth / $(GM/b)$', title='(b) Host potential', ylim=(0, None))
    axes[1].legend(fontsize=8, frameon=False)
    fine = tables['trajectories'][-len(summary['shell_convergence']):]
    pos = np.arange(len(fine))
    axes[2].bar(pos-.17, [r['snapshot_return_myr'] or 0 for r in fine], .32, color='#333333', label='Frozen snapshot')
    axes[2].bar(pos+.17, [r['plummer_return_myr'] for r in fine], .32, facecolor='white', edgecolor='#14677d', hatch='///', label='Matched Plummer')
    axes[2].set(xticks=pos, xticklabels=[str(r['eta_plummer']) for r in fine],
                xlabel=r'Kick / Plummer escape speed', ylabel='First ballistic return (Myr)', title='(c) Identical physical kicks')
    axes[2].legend(fontsize=8, frameon=False)
    for ax in axes:
        ax.grid(axis='y', color='#dddddd', linewidth=.4)
        ax.set_axisbelow(True)
    fig.savefig(out/'snapshot_comparison.png', dpi=200)
    fig.savefig(out/'snapshot_comparison.pdf', metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)


def main():
    dirty = subprocess.check_output(['git', 'status', '--porcelain', '--', *SOURCES], cwd=ROOT, text=True)
    if dirty.strip():
        raise SystemExit('Commit E004 source inputs before generating recorded results.')
    config = json.loads((ROOT/'experiments/e004.json').read_text())
    tables, summary = calculate(config)
    out = ROOT/'results/e004'
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in tables.items():
        write_csv(out/(name+'.csv'), rows)
    (out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    plot(tables, summary, out)
    manifest = {'experiment': 'e004', 'researcher': 'Venkata Hemanth Athota',
                'run_utc': datetime.now(timezone.utc).isoformat(), 'python': platform.python_version(),
                'versions': {n: importlib.metadata.version(n) for n in ['numpy', 'scipy', 'matplotlib']},
                'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'protocol_revision': '50cf824', 'external_source': config['source_repository'],
                'external_revision': config['source_revision'],
                'external_sha256': {key: value['sha256'] for key, value in config['files'].items()},
                'source_sha256': {name: sha(ROOT/name) for name in SOURCES},
                'output_sha256': {str(path.relative_to(ROOT)): sha(path) for path in sorted(out.iterdir())
                                  if path.suffix in ['.csv', '.json', '.png', '.pdf']},
                'rows': {name: len(rows) for name, rows in tables.items()}}
    (ROOT/'results/run_e004.json').write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
