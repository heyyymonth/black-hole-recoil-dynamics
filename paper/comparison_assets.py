"""Typeset recorded six-cluster findings; does not execute physical calculations."""
from pathlib import Path
import argparse,csv,hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1];PAPER=ROOT/'paper'
INPUTS=['paper/comparison_assets.py','results/cluster-comparison/reference-comparison.csv','results/cluster-comparison/posterior-inventory.json','results/cluster-comparison/summary.json','results/cluster-comparison/provenance.json','results/cluster-comparison/diffusion-comparison.pdf','results/cluster-comparison-environment/comparison.csv']
OUTPUTS=['tables/cluster-transport-comparison.tex','figures/cluster-transport-comparison.pdf']
NAMES={'NGC0104':'47 Tuc','NGC1851':'NGC 1851','NGC3201':'NGC 3201','NGC5139':r'$\omega$ Cen','NGC6121':'M4','NGC6397':'NGC 6397 (standard)'}
def hashes(paths,base):return {p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in paths}
def count_label(value):
    rounded=float(f"{value:.3g}")
    return f"{rounded:,.0f}" if abs(rounded)>=1000 else f"{rounded:.3g}"

def table():
    with (ROOT/INPUTS[1]).open() as f:rows=list(csv.DictReader(f))
    published={r['cluster']:r for r in json.loads((ROOT/INPUTS[2]).read_text())['clusters']}
    body=r'''\begin{table*}[t]
\caption{Conditional standard-fit comparison. $N_{\rm BH}^{\rm ref}$ is the expected count in the selected maximum-quadrature-weight joint row. Published intervals are empirical 16th--84th percentiles of 4,000 released derived model realizations per cluster. $f_D$ is the BH contribution to the local diffusion trace at zero remnant velocity, with a common constant Coulomb logarithm. The probe mass is $60\,M_\odot$. These are model quantities, not directly observed counts or residence fractions. The NGC 6397 standard fit is a sensitivity case; its alternative no-BH chain is unavailable.}\label{tab:cluster-transport-comparison}
\centering\begin{tabular}{lrrrrr}
\toprule
Cluster & $N_{\rm BH}^{\rm ref}$ & Published $N_{\rm BH}$ & $m_{\rm r}/M_{\rm BH}^{\rm ref}$ & $f_D(0)$ (\%) & $f_D(0.1r_h)$ (\%) \\
\midrule
'''
    for row in rows:
        c=row['cluster'];q=published[c]['BH_expected_number'];far=100*float(row['BH_D_fraction_at_0p1rh_v0'])
        farstr=f'{far:.4g}' if far>=.0001 else r'$<10^{-4}$'
        body+=f"{NAMES[c]} & {count_label(float(row['BH_expected_number']))} & {count_label(q['median'])} [{count_label(q['p16'])}, {count_label(q['p84'])}] & {float(row['probe_to_BH_reservoir_mass_ratio']):.4g} & {100*float(row['BH_D_fraction_center_v0']):.3f} & {farstr}"+r' \\'+'\n'
    return body+r'''\bottomrule
\end{tabular}
\end{table*}
'''
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();path=PAPER/'comparison-assets.json'
    if args.verify:
        stored=json.loads(path.read_text());assert stored['input_sha256']==hashes(INPUTS,ROOT);assert stored['output_sha256']==hashes(OUTPUTS,PAPER)
        assert (PAPER/OUTPUTS[0]).read_text()==table();print('Verified six-cluster table and vector figure.');return
    (PAPER/OUTPUTS[0]).write_text(table());shutil.copyfile(ROOT/INPUTS[5],PAPER/OUTPUTS[1])
    path.write_text(json.dumps({'scope':'conditional six-cluster host and local scattering comparison','input_sha256':hashes(INPUTS,ROOT),'output_sha256':hashes(OUTPUTS,PAPER)},indent=2)+'\n')
if __name__=='__main__':main()
