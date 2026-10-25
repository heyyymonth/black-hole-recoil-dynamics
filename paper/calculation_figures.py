"""Vector figures of stored reconstruction and kernel records."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def rows(root,path):
    with (root/path).open(newline='') as f:return list(csv.DictReader(f))

def render(root,out):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    data=rows(root,'results/host-reconstruction/published-profile-comparison.csv')
    fig,axes=plt.subplots(1,3,figsize=(7,2.5),sharex=True)
    for ax,key,title in zip(axes,['LOS','pm_R','pm_T'],['Line of sight','Radial proper motion','Tangential proper motion']):
        selected=[r for r in data if r['profile']==key]
        x=np.array([float(r['R_pc']) for r in selected])
        low,med,high,pred=[np.array([float(r[k]) for r in selected]) for k in ['published_lower_95_45','published_median','published_upper_95_45','reconstructed']]
        ax.fill_between(x,low,high,color='#aec5d5',alpha=.55,label='Published 95.45% envelope')
        ax.plot(x,med,color='#245b82',ls='--',lw=1,label='Published median')
        ax.plot(x,pred,color='#a95b17',lw=1.3,label='Fixed joint-row reconstruction')
        ax.set_xscale('log');ax.set_xlabel('Projected radius (pc)');ax.set_title(title,fontsize=8.5);ax.grid(alpha=.15)
    axes[0].set_ylabel(r'Dispersion (km s$^{-1}$)');axes[0].legend(fontsize=6,frameon=False,loc='lower left')
    fig.tight_layout(pad=.7,w_pad=.8);fig.savefig(out/'host-projections.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    host=dict(np.load(root/'results/host-reconstruction/host.npz',allow_pickle=False))
    species=rows(root,'results/local-transport-kernels/species.csv');states=rows(root,'results/local-transport-kernels/states.csv')
    colors=['#245b82','#c28c28','#72767a','#8d4566'];types=['MS','WD','NS','BH'];labels=['Main sequence','White dwarfs','Neutron stars','Black holes']
    fig,axes=plt.subplots(1,2,figsize=(7,2.75))
    selection=(host['r']>0)&(host['r']/float(host['rh'])<=3)
    x=host['r'][selection]/float(host['rh']);rho=host['rhoj'][:,selection]
    fractions=np.array([rho[host['types']==t].sum(axis=0)/rho.sum(axis=0) for t in types])
    axes[0].stackplot(x,fractions,labels=labels,colors=colors);axes[0].set_xscale('log');axes[0].set_xlim(1e-4,3);axes[0].set_ylim(0,1);axes[0].set_xlabel(r'Radius $r/r_h$');axes[0].set_ylabel('Local mass-density fraction');fig.legend(*axes[0].get_legend_handles_labels(),frameon=False,fontsize=7,loc='upper center',ncol=4,bbox_to_anchor=(.5,1.01))
    rest=[s for s in states if float(s['v_over_cutoff'])==0];bottom=np.zeros(len(rest))
    for typ,col in zip(types,colors):
        vals=[sum(float(r['trace_D_fraction']) for r in species if r['state']==s['state'] and r['type']==typ) for s in rest]
        axes[1].bar(np.arange(len(rest)),vals,bottom=bottom,color=col,width=.65);bottom+=vals
    axes[1].set_xticks(np.arange(len(rest)),[f"{float(s['r_over_rh']):g}" for s in rest]);axes[1].set_ylim(0,1);axes[1].set_xlabel(r'Radius $r/r_h$ at zero test velocity');axes[1].set_ylabel('Fraction of diffusion trace')
    fig.tight_layout(pad=.7,w_pad=1.2,rect=(0,0,1,.92));fig.savefig(out/'scattering-composition.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
