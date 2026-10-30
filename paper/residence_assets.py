"""Present collisionless reference paths and encounter integrals, without reruns."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT/'paper'
RESULTS = ROOT/'results/collisionless-residence'
CLUSTERS = ('NGC0104', 'NGC1851', 'NGC3201', 'NGC5139', 'NGC6121', 'NGC6397')
ORIENTATIONS = ('radial_outward', 'tangential')
NAMES = {'NGC0104': '47 Tuc', 'NGC1851': 'NGC 1851', 'NGC3201': 'NGC 3201',
         'NGC5139': r'$\omega$ Cen', 'NGC6121': 'M4', 'NGC6397': 'NGC 6397 (standard)'}
INPUTS = ['paper/residence_assets.py', 'analysis/residence/protocol.json',
          'analysis/encounters/protocol.json',
          'results/collisionless-residence/summary.json',
          'results/collisionless-residence/comparison.csv',
          *(f'results/collisionless-residence/{c}.json' for c in CLUSTERS)]
OUTPUTS = ['paper/tables/collisionless-residence.tex',
           'paper/figures/collisionless-residence.pdf',
           'results/collisionless-residence/diagnostics.pdf',
           'results/collisionless-residence/diagnostics.png']


def hashes(paths):
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths}


def read_rows():
    with (RESULTS/'comparison.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((RESULTS/'summary.json').read_text())
    if (len(rows) != 12 or summary['trajectories'] != 12
            or summary['accepted_trajectories'] != 12
            or not summary['all_numerical_checks_passed']):
        raise ValueError('All twelve declared reference trajectories must be accepted')
    if (not summary['collisionless_trajectories_performed']
            or summary['collisional_trajectories_performed']
            or summary['residence_response_evaluated']):
        raise ValueError('Presentation scope requires collisionless controls only')
    selected = {}
    for cluster in CLUSTERS:
        data = json.loads((RESULTS/f'{cluster}.json').read_text())
        recorded = {r['summary']['orientation']: r['summary'] for r in data['records']}
        if set(recorded) != set(ORIENTATIONS):
            raise ValueError('Unexpected launch inventory')
        for orientation in ORIENTATIONS:
            match = [r for r in rows if r['cluster']==cluster and r['orientation']==orientation]
            if len(match) != 1 or match[0]['accepted'] != 'True':
                raise ValueError('A unique accepted trajectory is required')
            row = match[0]
            for key, value in recorded[orientation].items():
                if row[key] != str(value):
                    raise ValueError(f'CSV and JSON differ: {cluster}, {orientation}, {key}')
            if not 0 <= float(row['residence_fraction']) <= 1:
                raise ValueError('Residence fraction lies outside the physical interval')
            if float(row['allowed_rate_integral']) <= 0:
                raise ValueError('Positive encounter integral is required for this plot')
            selected[cluster, orientation] = row
    return selected


def table(rows):
    text = r'''\begin{table*}[t]
\caption{Collisionless reference residence and BH encounter integrals. Each launch begins at $0.5r_\rho$ with the speed giving a radial apocenter of $2r_\rho$; the tangential launch has the same energy. Residence counts every visit to $r<r_\rho$ during $T=10t_{\rm dyn}$. The gravitational field includes all physical host species. The allowed BH elastic rate uses $b_{\max}=0.1L$ and the declared close-passage exclusion, evaluated along each prescribed path. Its time integral is a diagnostic that does not alter the orbit. These are conditional collisionless controls, not collisional residence predictions or merger probabilities.}\label{tab:collisionless-residence}
\centering\begin{tabular}{lrrrrr}
\toprule
 & & \multicolumn{2}{c}{$C^{(0)}_T/T$} & \multicolumn{2}{c}{$\int_0^{T_*}\Gamma_{\rm el}\,dt$} \\
Cluster & $T$ (Myr) & Radial & Tangential & Radial & Tangential \\
\midrule
'''
    for cluster in CLUSTERS:
        radial, tangential = (rows[cluster, orientation] for orientation in ORIENTATIONS)
        text += (f"{NAMES[cluster]} & {float(radial['horizon_Myr']):.3g} & {float(radial['residence_fraction']):.4f} & "
                 f"{float(tangential['residence_fraction']):.4f} & "
                 f"{float(radial['allowed_rate_integral']):.3g} & "
                 f"{float(tangential['allowed_rate_integral']):.3g}" + r' \\' + '\n')
    return text + r'''\bottomrule
\end{tabular}
\end{table*}
'''


def figures(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                                'pdf.fonttype': 42, 'ps.fonttype': 42,
                                'axes.spines.top': False, 'axes.spines.right': False,
                                'axes.linewidth': .65, 'lines.linewidth': 1.,
                                'savefig.facecolor': 'white'})
    fig, axes = plt.subplots(2, 1, figsize=(7., 4.5), sharex=True,
                             gridspec_kw={'hspace': .13})
    x = np.arange(len(CLUSTERS))
    for orientation, color, marker, offset, label in [
            ('radial_outward', '#0072B2', 'o', -.09, 'Radial outward'),
            ('tangential', '#D55E00', 's', .09, 'Tangential')]:
        residence = [float(rows[c, orientation]['residence_fraction']) for c in CLUSTERS]
        integral = [float(rows[c, orientation]['allowed_rate_integral']) for c in CLUSTERS]
        axes[0].plot(x+offset, residence, linestyle='none', marker=marker, markersize=5,
                     color=color, label=label)
        axes[1].plot(x+offset, integral, linestyle='none', marker=marker, markersize=5,
                     color=color)
    axes[0].set_ylabel(r'Residence fraction $C^{(0)}_T/T$')
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc='upper left', ncol=2, fontsize=8, frameon=False,
                   bbox_to_anchor=(0., 1.23), borderaxespad=0.)
    axes[1].set_yscale('log')
    axes[1].set_ylabel(r'Allowed BH integral $\int\Gamma_{\rm el}\,dt$')
    axes[1].set_xticks(x, [NAMES[c].replace(' (standard)', '\n(standard)') for c in CLUSTERS])
    axes[0].text(.99, .95, '(a) Collisionless reference paths',
                 transform=axes[0].transAxes, ha='right', va='top', fontsize=8)
    axes[1].text(.015, .95, r'(b) Path diagnostic, $b_{\max}=0.1L$',
                 transform=axes[1].transAxes, ha='left', va='top', fontsize=8)
    for ax in axes:
        ax.set_xlim(-.55, len(CLUSTERS)-.45)
        ax.grid(axis='y', color='#dddddd', linewidth=.5)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.12, right=.98, top=.89, bottom=.13)
    metadata = {'Creator': 'collisionless reference presentation',
                'CreationDate': None, 'ModDate': None}
    fig.savefig(RESULTS/'diagnostics.pdf', metadata=metadata)
    fig.savefig(RESULTS/'diagnostics.png', dpi=180)
    plt.close(fig)
    shutil.copyfile(RESULTS/'diagnostics.pdf', PAPER/'figures/collisionless-residence.pdf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = PAPER/'residence-assets.json'
    rows = read_rows()
    if args.verify:
        stored = json.loads(path.read_text())
        assert stored['input_sha256'] == hashes(INPUTS)
        assert stored['output_sha256'] == hashes(OUTPUTS)
        assert (PAPER/'tables/collisionless-residence.tex').read_text() == table(rows)
        assert (PAPER/'figures/collisionless-residence.pdf').read_bytes() == (RESULTS/'diagnostics.pdf').read_bytes()
        print('Verified collisionless residence presentation; no calculations rerun.')
        return
    (PAPER/'tables/collisionless-residence.tex').write_text(table(rows))
    figures(rows)
    path.write_text(json.dumps({'scope': 'Collisionless reference paths and BH encounter integrals; no collisional dynamics',
                               'input_sha256': hashes(INPUTS), 'output_sha256': hashes(OUTPUTS)},
                              indent=2)+'\n')
    print('Generated collisionless residence presentation from recorded calculations.')


if __name__ == '__main__':
    main()
