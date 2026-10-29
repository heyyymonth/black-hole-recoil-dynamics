"""Present recorded finite-angle encounter calculations; never evolve a path.

The read-only verifier uses the standard library. Regeneration needs the
repository's Matplotlib environment and creates vector PDF presentation assets.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'paper'
RESULTS = ROOT / 'results/finite-angle-encounters'
CLUSTERS = ('NGC0104', 'NGC1851', 'NGC3201', 'NGC5139', 'NGC6121', 'NGC6397')
NAMES = {'NGC0104': '47 Tuc', 'NGC1851': 'NGC 1851', 'NGC3201': 'NGC 3201',
         'NGC5139': r'$\omega$ Cen', 'NGC6121': 'M4', 'NGC6397': 'NGC 6397 (standard)'}
INPUTS = ['paper/encounter_assets.py', 'analysis/encounters/protocol.json',
          'results/finite-angle-encounters/summary.json',
          'results/finite-angle-encounters/states.csv',
          *(f'results/finite-angle-encounters/{c}.json' for c in CLUSTERS)]
OUTPUTS = ['paper/tables/finite-angle-encounters.tex',
           'paper/figures/finite-angle-encounters.pdf',
           'results/finite-angle-encounters/diagnostics.pdf',
           'results/finite-angle-encounters/diagnostics.png']


def hashes(paths):
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths}


def read_rows():
    with (RESULTS/'states.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    summary = json.loads((RESULTS/'summary.json').read_text())
    if len(rows) != summary['states'] or summary['clusters'] != len(CLUSTERS):
        raise ValueError('Recorded summary and state inventory differ')
    if summary['transport_simulations_performed'] or summary['residence_response_evaluated']:
        raise ValueError('These presentation assets concern local encounters only')
    for cluster in CLUSTERS:
        recorded = json.loads((RESULTS/f'{cluster}.json').read_text())['states']
        selected = [row for row in rows if row['cluster'] == cluster]
        if len(selected) != len(recorded):
            raise ValueError('CSV and host state inventory differ')
        for csv_row, json_row in zip(selected, recorded, strict=True):
            for key in json_row:
                if str(json_row[key]) != csv_row[key]:
                    raise ValueError(f'CSV and JSON differ: {cluster}, {key}')
    return rows


def selected_rows(rows, eta):
    selected = {}
    for cluster in CLUSTERS:
        match = [row for row in rows if row['cluster'] == cluster
                 and float(row['r_over_density_radius']) == 1.
                 and float(row['speed_over_sigma']) == 1.
                 and row['orientation'] == 'radial'
                 and float(row['bmax_over_local_scale']) == eta]
        if len(match) != 1 or match[0]['accepted'] != 'True':
            raise ValueError(f'A unique accepted presentation state is required: {cluster}, {eta}')
        selected[cluster] = match[0]
    return selected


def table(rows):
    selected = selected_rows(rows, .1)
    text = r'''\begin{table*}[t]
\caption{Finite-angle BH encounters in the adopted reference hosts. At radius $r_\rho$, BH mass density is half its central value; this is not a half-mass radius. Each row uses $r=r_\rho$, radial remnant speed $v=\sigma_{\rm BH}$, remnant mass $60\,M_\odot$ and $b_{\max}=0.1L$. The one-dimensional scale $\sigma_{\rm BH}$ is number weighted, and $L$ is the smaller of $r_\rho$ and the local BH density scale. $\Gamma_{\rm el}$ and $\langle|\Delta\mathbf v|^2\rangle_{\rm el}^{1/2}$ refer to the allowed elastic domain. $f_{90}$ is the fraction of its diffusion trace arising from relative velocities with $b_{90}>L$. The latter is a scale diagnostic, not an error bound or a residence prediction.}\label{tab:finite-angle-encounters}
\centering\begin{tabular}{lrrrr}
\toprule
Cluster & $r_\rho$ (pc) & $\Gamma_{\rm el}$ (Myr$^{-1}$) & $\langle|\Delta\mathbf v|^2\rangle_{\rm el}^{1/2}/\sigma_{\rm BH}$ & $f_{90}$ (\%) \\
\midrule
'''
    for cluster, row in selected.items():
        text += (f"{NAMES[cluster]} & {float(row['density_radius_pc']):.4g} & "
                 f"{float(row['allowed_rate_per_Myr']):.4g} & "
                 f"{float(row['rms_increment_over_sigma']):.4g} & "
                 f"{100*float(row['b90_above_L_trace_fraction']):.2f}" + r' \\' + '\n')
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
                             gridspec_kw={'hspace': .12})
    x = np.arange(len(CLUSTERS))
    variants = [(.05, '#0072B2', 'o', -.17), (.1, '#D55E00', 's', 0.),
                (.2, '#009E73', '^', .17)]
    for eta, color, marker, offset in variants:
        selected = selected_rows(rows, eta)
        rates = [float(selected[c]['allowed_rate_per_Myr']) for c in CLUSTERS]
        fractions = [100*float(selected[c]['b90_above_L_trace_fraction']) for c in CLUSTERS]
        axes[0].plot(x+offset, rates, linestyle='none', marker=marker, markersize=5,
                     color=color, label=rf'$b_{{\max}}/L={eta:g}$')
        axes[1].plot(x+offset, fractions, linestyle='none', marker=marker, markersize=5,
                     color=color)
    axes[0].set_yscale('log')
    axes[0].set_ylabel(r'Allowed rate $\Gamma_{\rm el}$ (Myr$^{-1}$)')
    axes[0].legend(loc='upper left', ncol=3, fontsize=8, frameon=False,
                   bbox_to_anchor=(0., 1.23), borderaxespad=0.)
    axes[1].set_ylabel(r'Diffusion trace from $b_{90}>L$ (%)')
    axes[1].set_ylim(-3, 103)
    axes[1].set_yticks([0, 25, 50, 75, 100])
    axes[1].set_xticks(x, [NAMES[c].replace(' (standard)', '\n(standard)') for c in CLUSTERS])
    axes[0].text(.99, .95, r'(a) $r=r_\rho$, radial $v=\sigma_{\rm BH}$',
                 transform=axes[0].transAxes, ha='right', va='top', fontsize=8)
    axes[1].text(.015, .95, '(b) Conditional scale diagnostic',
                 transform=axes[1].transAxes, ha='left', va='top', fontsize=8)
    for ax in axes:
        ax.set_xlim(-.55, len(CLUSTERS)-.45)
        ax.grid(axis='y', color='#dddddd', linewidth=.5)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.12, right=.98, top=.89, bottom=.13)
    metadata = {'Creator': 'finite-angle encounter presentation',
                'CreationDate': None, 'ModDate': None}
    fig.savefig(RESULTS/'diagnostics.pdf', metadata=metadata)
    fig.savefig(RESULTS/'diagnostics.png', dpi=180)
    plt.close(fig)
    shutil.copyfile(RESULTS/'diagnostics.pdf', PAPER/'figures/finite-angle-encounters.pdf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    path = PAPER/'encounter-assets.json'
    rows = read_rows()
    if args.verify:
        stored = json.loads(path.read_text())
        assert stored['input_sha256'] == hashes(INPUTS)
        assert stored['output_sha256'] == hashes(OUTPUTS)
        assert (PAPER/'tables/finite-angle-encounters.tex').read_text() == table(rows)
        assert (PAPER/'figures/finite-angle-encounters.pdf').read_bytes() == (RESULTS/'diagnostics.pdf').read_bytes()
        print('Verified finite-angle encounter table and vector figure; no trajectories computed.')
        return
    (PAPER/'tables/finite-angle-encounters.tex').write_text(table(rows))
    figures(rows)
    path.write_text(json.dumps({'scope': 'Conditional local finite-angle BH encounter rates and scale diagnostics',
                               'input_sha256': hashes(INPUTS), 'output_sha256': hashes(OUTPUTS)},
                              indent=2)+'\n')
    print('Generated finite-angle encounter presentation from recorded calculations.')


if __name__ == '__main__':
    main()
