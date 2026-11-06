"""Present retained precision-pilot outcomes; never execute trajectories."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/precision'
PAPER = ROOT/'paper'
LAUNCHES = [('radial_outward','Radial'),('tangential','Tangential'),('oblique','Oblique')]
OUTPUTS = ['paper/tables/precision-pilot.tex','paper/figures/precision-pilot.pdf',
           'results/precision/comparison.pdf','results/precision/comparison.png']
INPUTS = ['results/precision/summary.json','analysis/precision/protocol.json','paper/precision_assets.py']


def hashes(names):
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def comparisons(s):
    rows=[]
    for launch,label in LAUNCHES:
        subset=[r for r in s['comparisons'] if r['launch']==launch]
        if len(subset)!=2:raise ValueError('Incomplete precision comparison')
        a=next(r for r in subset if r['first']=='eigen_uniform')
        b=next(r for r in subset if r['first']=='symmetric_uniform')
        rows.append((launch,label,a,b))
    return rows


def ratio_cell(r):
    if r['estimate'] is None:return 'Unresolved'
    interval=r['interval']
    def number(x):
        if x!=0 and (abs(x)>=1000 or abs(x)<.01):
            mantissa,exponent=f'{x:.2e}'.split('e')
            return mantissa+r'\times10^{'+str(int(exponent))+'}'
        return f'{x:.2f}'
    if interval is None:return f"${number(r['estimate'])}$ [unresolved]"
    return (r'\shortstack[r]{'+f"${number(r['estimate'])}$"+r'\\'
            +f"$[{number(interval[0])}, {number(interval[1])}]$"+'}')


def table(s):
    t=r'''\begin{table*}[t]
\caption{Precision comparison at the same fine-table setting and $\Delta t=0.01t_{\rm dyn}$. Each method uses 4096 independent response quartets per launch, divided into four execution batches. E denotes eigenbasis coupling with uniform insertion, S symmetric-root coupling with uniform insertion, and T symmetric-root coupling with proportional time stratification. Standard errors concern the estimated response slope. Variance ratios above one favor the denominator method. Bracketed intervals are exploratory percentile bootstrap intervals adjusted over six ratios; they do not cover unseen rare scores, numerical bias or physical-model uncertainty.}\label{tab:precision-pilot}
\centering\small
\begin{tabular}{lrrrrr}
\toprule
Launch & SE(E) & SE(S) & SE(T) & Variance E/S [interval] & Variance S/T [interval] \\
\midrule
'''
    for launch,label,a,b in comparisons(s):
        rows={r['method']:r for r in s['rows'] if r['launch']==launch}
        errors=[rows[m]['estimate']['standard_error'] for m in ['eigen_uniform','symmetric_uniform','symmetric_stratified']]
        t+=(f'{label} & '+ ' & '.join(f'${v:.5f}$' for v in errors)+' & '
            +ratio_cell(a['efficiency']['variance_ratio'])+' & '+ratio_cell(b['efficiency']['variance_ratio'])+r' \\[3pt]'+'\n')
    return t+r'''\bottomrule
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
    fig,axes=plt.subplots(1,2,figsize=(7.,3.7),sharey=True)
    specs=[]
    for _,label,a,b in comparisons(s):
        specs.extend([(label+' E/S',a),(label+' S/T',b)])
    for ax,key,title in zip(axes,['variance_ratio','cost_adjusted'],
        ['(a) Variance efficiency','(b) Including recorded CPU cost']):
        for i,(label,row) in enumerate(specs):
            r=row['efficiency'][key]
            if r is None or r['estimate'] is None:
                ax.text(1,i,'Unresolved',ha='left',va='center');continue
            color='#0072B2' if row['first']=='eigen_uniform' else '#D55E00'
            marker='o' if row['first']=='eigen_uniform' else 's'
            value=r['estimate'];interval=r['interval']
            if value<=0:
                ax.text(1,i,'Zero ratio',ha='left',va='center',fontsize=7)
                continue
            ax.plot(value,i,marker,color=color,markersize=4)
            if interval is not None and interval[0]>0:
                # Percentile intervals need not contain the original estimate.
                ax.plot(interval,[i,i],color=color,linewidth=1)
                ax.plot(interval,[i,i],'|',color=color,markersize=6)
            elif interval is not None:
                ax.annotate('Interval reaches zero',(value,i),xytext=(3,5),
                            textcoords='offset points',fontsize=6,color=color)
            else:
                ax.annotate('Interval unresolved',(value,i),xytext=(3,5),
                            textcoords='offset points',fontsize=6,color=color)
        ax.set_xscale('log');ax.axvline(1,color='#777777',ls='--',lw=.8)
        ax.set_title(title,loc='left',fontsize=9);ax.set_xlabel('Efficiency ratio (above 1 favors denominator)')
        ax.grid(axis='x',color='#dddddd',lw=.5)
    axes[0].set_yticks(range(len(specs)),[label for label,_ in specs]);axes[0].invert_yaxis()
    fig.subplots_adjust(left=.18,right=.985,top=.90,bottom=.19,wspace=.18)
    fig.savefig(OUT/'comparison.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(OUT/'comparison.png',dpi=180);plt.close(fig)
    shutil.copyfile(OUT/'comparison.pdf',PAPER/'figures/precision-pilot.pdf')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    s=json.loads((OUT/'summary.json').read_text());path=PAPER/'precision-assets.json'
    if not s['complete'] or s['counts']['ensembles']!=36:raise ValueError('Complete pilot required')
    if args.verify:
        m=json.loads(path.read_text())
        if m['input_sha256']!=hashes(INPUTS) or m['output_sha256']!=hashes(OUTPUTS):raise ValueError('Changed precision presentation')
        if (PAPER/'tables/precision-pilot.tex').read_text()!=table(s):raise ValueError('Table differs from recorded summary')
        if (OUT/'comparison.pdf').read_bytes()!=(PAPER/'figures/precision-pilot.pdf').read_bytes():raise ValueError('Different vector figure copy')
        print('Verified precision-pilot table and figure from retained evidence.');return
    (PAPER/'tables/precision-pilot.tex').write_text(table(s));figures(s)
    path.write_text(json.dumps({'input_sha256':hashes(INPUTS),'output_sha256':hashes(OUTPUTS),
        'scope':s['scope']},indent=2)+'\n')
    print('Generated precision-pilot presentation without dynamics.')


if __name__=='__main__':main()
