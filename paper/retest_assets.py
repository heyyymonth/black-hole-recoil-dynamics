"""Present retained response predictions and independent validation; no dynamics."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/retest'
PAPER = ROOT/'paper'
OUTPUTS = ['paper/tables/independent-response.tex', 'paper/tables/response-continuation.tex',
           'paper/figures/independent-response.pdf',
           'results/retest/comparison.pdf', 'results/retest/comparison.png']
LAUNCHES = [('radial_outward', 'Radial'), ('tangential', 'Tangential'),
            ('oblique', 'Oblique')]
VERDICTS = {'within_tolerance': 'Within', 'outside_tolerance': 'Outside',
            'inconclusive': 'Inconclusive'}


def hashes(paths):
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def inputs(summary):
    return sorted({'paper/retest_assets.py', 'results/retest/summary.json',
                   *summary['input_sha256']})


def primary_rows(summary):
    """Require a complete recorded primary comparison; never invent a row."""
    rows = []
    for name, label in LAUNCHES:
        predictions = [r for r in summary['predictions']
                       if r['role'] == 'primary' and r['launch'] == name]
        validation = sorted((r for r in summary['validations']
                             if r['role'] == 'primary' and r['launch'] == name),
                            key=lambda r: r['a'])
        if len(predictions) != 1 or [r['a'] for r in validation] != [.05, .1, 1.]:
            raise ValueError('Incomplete primary response evidence: '+name)
        rows.append((label, predictions[0], validation))
    return rows


def cell(stat):
    low, high = stat['interval']
    return (r"\shortstack[r]{"+f"${stat['estimate']:+.4f}$"+r"\\"
            +f"$[{low:+.4f}, {high:+.4f}]$"+"}")


def table(summary):
    text = r'''\begin{table*}[t]
\caption{Symmetric-root first-response prediction and independent endpoint validation in the conditional $\omega$ Cen reference, using the refined encounter table and $\Delta t=0.01t_{\rm dyn}$. All quantities concern normalized residence $C_T/T$. Each response uses 65,536 independent weighted samples; each endpoint shift uses 8192 mixture/diffusion pairs independent of the prediction. Intervals allocate error probability $0.025$ over a conservative Bonferroni family of 64 normal approximations. A separate $0.025$ budget covers 24 population continuation bounds. The endpoint error is the measured shift minus $S_h$; its classification uses the declared $0.001$ absolute tolerance. These sampling intervals exclude numerical refinement and continuation uncertainty, so a finite-resolution classification is not a continuum or physical-cluster verdict.}\label{tab:independent-response}
\centering\small
\begin{tabular}{lrrrl}
\toprule
Launch & $S_h$ [interval] & Endpoint shift [interval] & Prediction error [interval] & Test \\
\midrule
'''
    for label, prediction, validation in primary_rows(summary):
        endpoint = validation[-1]
        text += (f"{label} & "
                 f"{cell(prediction['response'])} & {cell(endpoint['observed_shift'])} & "
                 f"{cell(endpoint['prediction_error'])} & "
                 f"{VERDICTS[endpoint['finite_model_classification']]}"+r' \\[3pt]'+'\n')
    return text+r'''\bottomrule
\end{tabular}
\end{table*}
'''



def continuation_table(summary):
    text = r'''\begin{table*}[t]
\caption{Numerical-continuation allowances at the primary resolution. $\overline B$ is the finite-sample mean allowance; $U_B$ is the population-mean upper bound obtained from a one-sided binomial occurrence limit with $\delta=0.025/24$. The deterministic support maxima are 2000 for response quartets and 2 for validation pairs. Endpoint validation uses $a=1$. These conservative bounds do not supply dynamics beyond the numerical domain. The continuation budget is $0.0002$; a zero observed allowance is not zero population uncertainty.}\label{tab:response-continuation}
\centering\small
\begin{tabular}{lrrrr}
\toprule
Launch & Response $\overline B$ & Response $U_B$ & Endpoint $\overline B$ & Endpoint $U_B$ \\
\midrule
'''
    for label,prediction,validation in primary_rows(summary):
        a=prediction['continuation'];b=validation[-1]['continuation']
        def f(x):
            if x==0:return '$0$'
            if x<.001:
                m,e=f'{x:.2e}'.split('e')
                return '$'+m+r'\times10^{'+str(int(e))+'}$'
            return f'${x:.5f}$'
        text+=label+' & '+' & '.join(f(v) for v in [a['sample_mean_bound'],a['upper_bound'],b['sample_mean_bound'],b['upper_bound']])+r' \\'+ '\n'
    return text+r'''\bottomrule
\end{tabular}
\end{table*}
'''


def figures(summary):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
        'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .65})
    fig, axes = plt.subplots(2, 3, figsize=(7., 5.2), sharex=True, sharey='row')
    tolerance = summary['reporting']['prediction_absolute_error_target']
    for col, (label, prediction, validation) in enumerate(primary_rows(summary)):
        upper, lower = axes[:, col]
        response = prediction['response']
        a = np.array([r['a'] for r in validation])
        fine_a = np.linspace(0, 1, 101)
        upper.plot(fine_a, fine_a*response['estimate'], color='#0072B2',
                   linewidth=1.1, label='Symmetric-root response')
        upper.fill_between(fine_a, fine_a*response['interval'][0],
                           fine_a*response['interval'][1], color='#0072B2',
                           alpha=.16, linewidth=0)
        shift = np.array([r['observed_shift']['estimate'] for r in validation])
        shift_interval = np.array([r['observed_shift']['interval'] for r in validation])
        upper.errorbar(a, shift, yerr=np.stack([shift-shift_interval[:, 0],
            shift_interval[:, 1]-shift]), linestyle='none', marker='o', markersize=3.5,
            capsize=2, linewidth=.8, color='#D55E00', label='Independent mixture shift')
        residual = np.array([r['prediction_error']['estimate'] for r in validation])
        interval = np.array([r['prediction_error']['interval'] for r in validation])
        lower.axhspan(-tolerance, tolerance, color='#777777', alpha=.13,
                      linewidth=0, label='Endpoint tolerance')
        lower.errorbar(a, residual, yerr=np.stack([residual-interval[:, 0],
            interval[:, 1]-residual]), linestyle='none', marker='s', markersize=3.,
            capsize=2, linewidth=.8, color='#222222')
        upper.set_title(f'({chr(97+col)}) {label}', loc='left', fontsize=9)
        lower.set_title(f'({chr(100+col)}) Residual', loc='left', fontsize=9)
        for ax in (upper, lower):
            ax.axhline(0, color='#999999', linewidth=.6, linestyle=':')
            ax.grid(axis='y', color='#dddddd', linewidth=.5)
            ax.set_xlim(-.035, 1.045)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        lower.set_xticks([0, .25, .5, .75, 1])
        lower.set_xlabel('Mixture fraction $a$')
    axes[0, 0].set_ylabel('Mean residence shift')
    axes[1, 0].set_ylabel('Prediction residual')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    extra_handles, extra_labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles+extra_handles, labels+extra_labels, loc='upper center',
               ncol=2, fontsize=7.5, frameon=False, bbox_to_anchor=(.52, 1.01))
    fig.subplots_adjust(left=.105, right=.985, bottom=.09, top=.85,
                        hspace=.39, wspace=.20)
    fig.savefig(OUT/'comparison.pdf', metadata={'Creator': 'independent residence response presentation',
                                               'CreationDate': None, 'ModDate': None})
    fig.savefig(OUT/'comparison.png', dpi=180)
    plt.close(fig)
    shutil.copyfile(OUT/'comparison.pdf', PAPER/'figures/independent-response.pdf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    summary = json.loads((OUT/'summary.json').read_text())
    if not summary.get('complete'):
        raise ValueError('Completed recorded response validation is required')
    path = PAPER/'retest-assets.json'
    if args.verify:
        recorded = json.loads(path.read_text())
        if recorded['input_sha256'] != hashes(inputs(summary)):
            raise ValueError('Response-presentation input hashes changed')
        if recorded['output_sha256'] != hashes(OUTPUTS):
            raise ValueError('Response-presentation output hashes changed')
        if (PAPER/'tables/independent-response.tex').read_text() != table(summary):
            raise ValueError('Response table differs from retained summary')
        if (PAPER/'tables/response-continuation.tex').read_text()!=continuation_table(summary):
            raise ValueError('Continuation table differs from retained summary')
        if (PAPER/'figures/independent-response.pdf').read_bytes() != (OUT/'comparison.pdf').read_bytes():
            raise ValueError('Response figure copy differs')
        print('Verified response table and vector figure from retained samples.')
        return
    (PAPER/'tables/independent-response.tex').write_text(table(summary))
    (PAPER/'tables/response-continuation.tex').write_text(continuation_table(summary))
    figures(summary)
    path.write_text(json.dumps({'scope': summary['scope'],
        'input_sha256': hashes(inputs(summary)), 'output_sha256': hashes(OUTPUTS)},
        indent=2, allow_nan=False)+'\n')
    print('Generated response presentation without rerunning trajectories.')


if __name__ == '__main__':
    main()
