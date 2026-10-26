"""Typeset and plot the recorded Galactic reference-orbit audit."""
from pathlib import Path
import argparse,csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1];PAPER=ROOT/'paper'
INPUTS=['paper/orbit_assets.py','results/orbit/provenance.json','results/orbit/summary.json','results/orbit/tracks.csv','results/orbit/events.csv']
OUTPUTS=['tables/galactic-orbit.tex','figures/galactic-orbit.pdf']
def rows(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))
def table():
    s=json.loads((ROOT/'results/orbit/summary.json').read_text())
    labels=['Catalog distance','Selected-host distance','Catalog distance minus quoted error','Catalog distance plus quoted error']
    body=r'''\begin{table*}[t]
\caption{Galactic reference-orbit diagnostics over the defined 1 Gyr horizon in Irrgang Model I. The host mass profile is fixed. $D$ is heliocentric distance; mean pericenter and apocenter are spherical Galactocentric radii from detected turning events. The intervals show orbit-phase variation of the sufficient 1\% isolated-force radius $r_\epsilon$ and frozen-frame diagnostic $r_f$, not posterior uncertainties.}\label{tab:galactic-orbit}
\centering\begin{tabular}{lrrrrr}
\toprule
Case & $D$ (kpc) & $\langle R_p\rangle$ (kpc) & $\langle R_a\rangle$ (kpc) & $r_\epsilon$ (pc) & $r_f$ (pc) \\
\midrule
'''
    for label,c in zip(labels,s['cases']):
        body+=f"{label} & {c['distance_kpc']:.3f} & {c['pericenter_mean_kpc']:.2f} & {c['apocenter_mean_kpc']:.2f} & {c['r_sufficient_1pct_pc_min']:.1f}--{c['r_sufficient_1pct_pc_max']:.1f} & {c['frozen_frame_scale_pc_min']:.1f}--{c['frozen_frame_scale_pc_max']:.1f} "+r' \\'+'\n'
    return body+'\\bottomrule\n\\end{tabular}\n\\end{table*}\n'
def figure():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    data=[r for r in rows('results/orbit/tracks.csv') if r['case']=='host_distance']
    v=lambda k:np.array([float(r[k]) for r in data])
    t=v('t_Myr');distance=np.sqrt(v('X_kpc')**2+v('Y_kpc')**2+v('Z_kpc')**2)
    fig,axs=plt.subplots(2,2,figsize=(7,4.7),sharex=True);a,b,c,d=axs.ravel()
    blue,orange,purple='#245b82','#a95b17','#8d4566'
    a.plot(t,distance,color=blue,lw=1,label='Galactocentric distance');a.plot(t,np.abs(v('Z_kpc')),color=orange,lw=.9,ls='--',label='Height above/below disk');a.set_ylabel('Distance (kpc)');a.set_title('(a) Galactic reference orbit',fontsize=8.5)
    b.semilogy(t,v('frozen_frame_scale_pc'),color=purple,lw=.9,label='Frozen-frame scale');b.semilogy(t,v('r_sufficient_1pct_pc'),color=blue,lw=1,label='Sufficient 1% force radius');b.axhline(8,color='black',ls=':',lw=.8,label='Fixed 8 pc control');b.set_ylabel('Cluster-relative radius (pc)');b.set_title('(b) Distinct boundary diagnostics',fontsize=8.5)
    c.semilogy(t,v('eta_bound_3rh'),color=orange,lw=1,label=r'At $3r_h$');c.semilogy(t,v('eta_bound_rh'),color=blue,lw=1,label=r'At $r_h$');c.axhline(.01,color='black',ls=':',lw=.8,label='1% tolerance');c.set_ylabel('Bound on omitted-force ratio');c.set_title('(c) Isolated-force accuracy',fontsize=8.5)
    d.semilogy(t,v('timescale_ratio_frozen_scale'),color=purple,lw=.8,label=r'At $r_f$');d.semilogy(t,v('timescale_ratio_rh'),color=blue,lw=.9,label=r'At $r_h$');d.axhline(1,color='black',ls=':',lw=.8);d.set_ylim(.001,100);d.set_ylabel(r'$t_{\rm dyn}\,\|\dot{T}\|_2/\|T\|_2$');d.set_title('(d) Forcing versus host timescale',fontsize=8.5)
    for ax in axs.ravel():
        ax.set_xlim(0,1000);ax.grid(alpha=.15);ax.legend(loc='upper right',fontsize=6,frameon=True,facecolor='white',framealpha=.9,edgecolor='none')
    for ax in axs[-1]:ax.set_xlabel('Elapsed reference-orbit time (Myr)')
    fig.tight_layout(pad=.7,h_pad=1.2,w_pad=1.6)
    fig.savefig(PAPER/OUTPUTS[1],metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
def hashes(paths,base):return {p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in paths}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');a=parser.parse_args();m=PAPER/'orbit-assets.json'
    if a.verify:
        rec=json.loads(m.read_text());assert rec['input_sha256']==hashes(INPUTS,ROOT);assert rec['output_sha256']==hashes(OUTPUTS,PAPER);assert (PAPER/OUTPUTS[0]).read_text()==table();print('Verified Galactic-orbit table and vector figure.');return
    (PAPER/OUTPUTS[0]).write_text(table());figure();m.write_text(json.dumps({'scope':'recorded Galactic reference orbits and host-support diagnostics','input_sha256':hashes(INPUTS,ROOT),'output_sha256':hashes(OUTPUTS,PAPER)},indent=2)+'\n')
if __name__=='__main__':main()
