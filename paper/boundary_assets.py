"""Present recorded boundary results; --verify uses only the standard library."""
from pathlib import Path
import argparse,csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1];PAPER=ROOT/'paper'
INPUTS=['paper/boundary_assets.py','results/outer-boundary/provenance.json','results/outer-boundary/summary.json','results/outer-boundary/scenarios.csv','results/outer-boundary/radial-profile.csv']
OUTPUTS=['tables/outer-boundary.tex','figures/outer-boundary.pdf']

def rows(name):
    with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))

def table():
    body=r'''\begin{table}[t]
\caption{Conditional boundary calculations in a spherical logarithmic Galactic field, keeping the NGC 3201 host fixed. $r_J$ uses enclosed mass. $r_{1\%}$ bounds the maximum neglected Galactic acceleration by 1\% of the cluster acceleration. Cases are sensitivity settings, not uncertainty intervals; the 13 kpc case is a circular surrogate.}\label{tab:outer-boundary}
\centering\begin{tabular}{rrrrr}
\toprule
$R$ (kpc) & $v_c$ (km/s) & $r_J$ (pc) & $r_{1\%}$ (pc) & $\eta(3r_h)$ \\
\midrule
'''
    for row in rows('results/outer-boundary/scenarios.csv'):
        r={k:float(v) for k,v in row.items() if k!='scenario'}
        body+=f"{r['R_pc']/1000:.2f} & {r['vc_kms']:.0f} & {r['rJ_extended_pc']:.2f} & {r['r_valid_0.01_pc']:.2f} & {100*r['eta_at_3rh']:.2f}\\% "+r' \\'+'\n'
    return body+'\\bottomrule\n\\end{tabular}\n\\end{table}\n'

def figure():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    data=rows('results/outer-boundary/radial-profile.csv');s=json.loads((ROOT/'results/outer-boundary/summary.json').read_text());b=s['baseline']
    r=np.array([float(d['r_pc']) for d in data]);eta=np.array([float(d['max_external_to_internal_force']) for d in data]);mass=np.array([float(d['enclosed_mass_fraction']) for d in data])
    fig,axes=plt.subplots(1,2,figsize=(7,2.7))
    axes[0].loglog(r,eta,color='#245b82',label='Exact maximum force ratio')
    axes[0].axhline(.01,color='#a95b17',ls='--',lw=1,label='Chosen 1% tolerance')
    axes[0].axvline(b['r_valid_0.01_pc'],color='#a95b17',ls=':',lw=1)
    axes[0].plot([s['host_rh_pc'],3*s['host_rh_pc']],[b['eta_at_rh'],b['eta_at_3rh']],'o',color='#222222',ms=4,label=r'Local checks: $r_h$, $3r_h$')
    axes[0].set_xlim(.1,200);axes[0].set_ylim(1e-5,20);axes[0].set_ylabel(r'$\eta=\max |a_{\rm ext}|/|a_{\rm cl}|$');axes[0].legend(fontsize=6,frameon=False,loc='upper left')
    axes[1].semilogx(r,mass,color='#245b82',label='Unmodified host')
    axes[1].axvline(b['rJ_extended_pc'],color='#8d4566',ls='--',lw=1,label=r'Circular $r_J=84.13$ pc')
    axes[1].axvline(s['host_formal_rt_pc'],color='#72767a',ls=':',lw=1,label=r'Formal $r_t=981$ pc')
    axes[1].scatter([b['rJ_extended_pc']],[1-b['mass_outside_rJ_fraction']],s=12,color='#8d4566')
    axes[1].set_xlim(.1,1200);axes[1].set_ylim(0,1.04);axes[1].set_ylabel('Enclosed fraction of host mass');axes[1].legend(fontsize=6,frameon=False,loc='upper left')
    for ax in axes:ax.set_xlabel('Cluster-relative radius (pc)');ax.grid(alpha=.15)
    fig.tight_layout(pad=.7,w_pad=1.3)
    fig.savefig(PAPER/OUTPUTS[1],metadata={'CreationDate':None,'ModDate':None});plt.close(fig)

def hashes(paths,base):return {p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in paths}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args();path=PAPER/'boundary-assets.json'
    if args.verify:
        m=json.loads(path.read_text());assert m['input_sha256']==hashes(INPUTS,ROOT);assert m['output_sha256']==hashes(OUTPUTS,PAPER);assert (PAPER/OUTPUTS[0]).read_text()==table();print('Verified boundary table and vector figure.');return
    (PAPER/OUTPUTS[0]).write_text(table());figure();path.write_text(json.dumps({'scope':'recorded static boundary diagnostics','input_sha256':hashes(INPUTS,ROOT),'output_sha256':hashes(OUTPUTS,PAPER)},indent=2)+'\n')
if __name__=='__main__':main()
