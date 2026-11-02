"""Present the retained matched pilot; never run trajectories during a paper build."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from analysis.matched.analyze import name
OUT=ROOT/'results/matched-residence';PAPER=ROOT/'paper'


def hashes(paths):return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
def inputs(summary):
    return sorted({'paper/matched_assets.py','results/matched-residence/summary.json',
                   'results/collisionless-residence/NGC5139.json',*summary['input_sha256']})
OUTPUTS=['paper/tables/matched-residence.tex','paper/figures/matched-residence.pdf',
         'results/matched-residence/comparison.pdf','results/matched-residence/comparison.png']


def table(s):
    text=r'''\begin{table*}[t]
\caption{Conditional matched encounter pilot in the $\omega$ Cen reference. Each entry compares independent ensembles of 4096 paths per arm over $T=10t_{\rm dyn}$. The jump and diffusion collision increments have the same local mean and covariance. Only BH collisions are active; all physical species contribute smooth gravity. Means include stopped paths and use the full horizon. The signed contrast $\delta_C=(\langle C_T\rangle_J-\langle C_T\rangle_D)/T$ has a descriptive 95\% normal interval; the CDF discrepancy has a 95\% union-DKW interval. These are conditional numerical results, not calibrated cluster predictions.}\label{tab:matched-residence}
\centering\begin{tabular}{llrrrrr}
\toprule
Launch & Grid & $\Delta t/t_{\rm dyn}$ & $\langle C_T/T\rangle_J$ & $\langle C_T/T\rangle_D$ & $\delta_C$ [95\% interval] & $\Delta_F$ [95\% interval] \\
\midrule
'''
    for o in ['radial_outward','tangential']:
        for c in [c for c in s['comparisons'] if c['orientation']==o]:
            g,dt=c['grid'],c['dt'];j=s['settings'][name(g,dt,o,'jump')]['statistics'];d=s['settings'][name(g,dt,o,'diffusion')]['statistics']
            m=c['residence']['mean_difference'];k=c['residence']['empirical_cdf_distance']
            text+=(f"{'Radial' if o=='radial_outward' else 'Tangential'} & {g.capitalize()} & {dt:.2f} & "
                   f"{j['residence']['mean']:.4f} & {d['residence']['mean']:.4f} & "
                   f"${m['estimate']:+.4f}$ [$ {m['interval'][0]:+.4f}, {m['interval'][1]:+.4f}$] & "
                   f"{k['estimate']:.3f} [{k['interval'][0]:.3f}, {k['interval'][1]:.3f}]"+r' \\'+'\n')
    return text+r'''\bottomrule
\end{tabular}
\end{table*}
'''


def figures(s):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,
                               'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65})
    fig,axes=plt.subplots(2,2,figsize=(7.,5.0),gridspec_kw={'width_ratios':[1,1.15]})
    baseline=json.loads((ROOT/'results/collisionless-residence/NGC5139.json').read_text())
    base={r['summary']['orientation']:r['summary']['residence_fraction'] for r in baseline['records']}
    settings=[('coarse',.04),('fine',.04),('coarse',.02),('fine',.02)]
    for row,o in enumerate(['radial_outward','tangential']):
        left,right=axes[row]
        for method,color,marker,offset in [('jump','#0072B2','o',-.09),('diffusion','#D55E00','s',.09)]:
            means=[];errors=[]
            for g,dt in settings:
                m=s['settings'][name(g,dt,o,method)]['statistics']['residence']
                means.append(m['mean']);errors.append(m['interval'][1]-m['mean'])
            left.errorbar(np.arange(4)+offset,means,yerr=errors,color=color,marker=marker,
                          linestyle='none',capsize=2,markersize=3.5,label=method.capitalize())
            data=np.load(OUT/f'{name("fine",.02,o,method)}.npz',allow_pickle=False)
            x=np.sort(data['normalized_residence']);y=np.arange(1,len(x)+1)/len(x)
            right.step(x,y,where='post',color=color,linewidth=1.,label=method.capitalize())
        left.axhline(base[o],color='#777777',linestyle='--',linewidth=.8,label='Collisionless')
        right.axvline(base[o],color='#777777',linestyle='--',linewidth=.8)
        left.set_xticks(range(4),['C / .04','F / .04','C / .02','F / .02'],fontsize=7)
        left.set_ylabel('Mean residence fraction')
        left.set_title(('(a) Radial outward' if row==0 else '(c) Tangential'),loc='left',fontsize=9)
        right.set_title('(b) Refined-grid CDF' if row==0 else '(d) Refined-grid CDF',loc='left',fontsize=9)
        right.set_ylabel('Cumulative probability');right.set_xlim(0,1);right.set_ylim(0,1.02)
        left.grid(axis='y',color='#dddddd',linewidth=.5);right.grid(color='#dddddd',linewidth=.5)
        if row==1:
            left.set_xlabel('Grid / step in dynamical-time units')
            right.set_xlabel('Residence fraction')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=3,fontsize=8,frameon=False,
               bbox_to_anchor=(.52,1.),borderaxespad=.2)
    fig.subplots_adjust(left=.09,right=.985,bottom=.09,top=.89,wspace=.35,hspace=.34)
    fig.savefig(OUT/'comparison.pdf',metadata={'Creator':'matched residence presentation','CreationDate':None,'ModDate':None})
    fig.savefig(OUT/'comparison.png',dpi=180);plt.close(fig)
    shutil.copyfile(OUT/'comparison.pdf',PAPER/'figures/matched-residence.pdf')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    s=json.loads((OUT/'summary.json').read_text());path=PAPER/'matched-assets.json'
    if args.verify:
        stored=json.loads(path.read_text());assert stored['input_sha256']==hashes(inputs(s));assert stored['output_sha256']==hashes(OUTPUTS)
        assert (PAPER/'tables/matched-residence.tex').read_text()==table(s)
        assert (PAPER/'figures/matched-residence.pdf').read_bytes()==(OUT/'comparison.pdf').read_bytes()
        print('Verified matched-pilot table and vector figure from recorded ensembles.');return
    (PAPER/'tables/matched-residence.tex').write_text(table(s));figures(s)
    path.write_text(json.dumps({'scope':s['scope'],'input_sha256':hashes(inputs(s)),'output_sha256':hashes(OUTPUTS)},indent=2)+'\n')
    print('Generated matched-pilot presentation without rerunning trajectories.')
if __name__=='__main__':main()
