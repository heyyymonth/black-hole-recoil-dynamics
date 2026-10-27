"""Build comparison tables, a figure and research notes from recorded results only."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/cluster-comparison'
NAMES={'NGC0104':'47 Tuc','NGC1851':'NGC 1851','NGC3201':'NGC 3201','NGC5139':'Omega Centauri','NGC6121':'M4','NGC6397':'NGC 6397 (standard fit)'}

def load(path):return json.loads(path.read_text())
def build():
    allruns=load(OUT/'summary.json')['runs'];published=load(OUT/'posterior-inventory.json')['clusters']
    refs=[r for r in allruns if r['kind']=='reference'];flat=[]
    for r in refs:
        if r['host_status']!='passed':
            flat.append({'cluster':r['cluster'],'host_passed':False});continue
        h=r['host_summary']['summary'];m=r['BH_metrics'];p=r['core_profile']
        flat.append({'cluster':r['cluster'],'host_passed':True,'transport_passed':r['transport_status']=='passed',
                     'M_Msun':h['M_Msun'],'rh_pc':h['rh_pc'],'rt_pc':h['rt_pc'],**m,
                     'BH_D_fraction_center_v0':p[0]['BH_trace_D_fraction'],
                     'BH_D_fraction_at_0p1rh_v0':next(x['BH_trace_D_fraction'] for x in p if x['r_over_rh']==.1),
                     'central_trace_D_per_lnLambda':p[0]['trace_D_per_lnLambda']})
    keys=list(dict.fromkeys(k for r in flat for k in r))
    with (OUT/'reference-comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');w.writeheader();w.writerows(flat)
    fig,axes=plt.subplots(2,3,figsize=(7.1,4.5),sharex=True,sharey=True,layout='constrained')
    for ax,ref in zip(axes.flat,refs):
        c=ref['cluster'];pilots=[r for r in allruns if r['cluster']==c and r['kind']=='pilot']
        passed=0
        for r in pilots:
            if r['host_status']!='passed' or r.get('core_profile_status')!='passed':continue
            passed+=1
            p=r['core_profile'];ax.plot([x['r_over_rh'] for x in p],[100*x['BH_trace_D_fraction'] for x in p],color='#777777',alpha=.38,lw=.8)
        if ref['host_status']=='passed':
            p=ref['core_profile'];ax.plot([x['r_over_rh'] for x in p],[100*x['BH_trace_D_fraction'] for x in p],color='#1b4f72',lw=1.4,marker='o',ms=2.7,label='Reference')
        ax.set_xscale('symlog',linthresh=.001);ax.set_xlim(0,3);ax.set_ylim(-2,102)
        ax.set_xticks([0,.001,.01,.1,1,3]);ax.set_xticklabels(['0','.001','.01','.1','1','3'])
        ax.set_title(NAMES[c],fontsize=9)
        absent=sum(r.get('BH_metrics',{}).get('BH_absent',False) for r in pilots)
        note=f'{passed}/8 pilot rows'+(f'; {absent} BH-free' if absent else '')
        ax.text(.03,.08,note,ha='left',va='bottom',transform=ax.transAxes,fontsize=7)
        ax.grid(axis='y',color='#dddddd',lw=.5);ax.tick_params(labelsize=7)
        ax.spines[['top','right']].set_visible(False)
    for ax in axes[:,0]:ax.set_ylabel('BH fraction of diffusion trace (%)',fontsize=8)
    for ax in axes[1,:]:ax.set_xlabel(r'Radius $r/r_h$; remnant velocity $v=0$',fontsize=8)
    fig.savefig(OUT/'diffusion-comparison.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(OUT/'diffusion-comparison.png',dpi=180);plt.close(fig)
    lines=['# Six-cluster host and local-transport comparison','',
           'The same conditional host and local scattering experiments are evaluated for all six shortlisted clusters. These calculations test numerical consistency and expose differences between candidate environments; they do not compare diffusion with resolved encounters or calculate remnant residence.','',
           '## Calculation design','',
           '- One maximum-posterior-quadrature-weight joint row per cluster repeats the original reference selection. It is not a posterior mode or a typical host.',
           '- Each accepted reference repeats the same 20 local states: four radii, three speed fractions and radial/tangential directions, with a 60-solar-mass diagnostic probe.',
           '- Eight fixed-seed weighted joint draws per cluster form a sensitivity pilot. Failed rows are retained and not replaced. These eight draws do not establish converged credible intervals.',
           '- Nine zero-speed radii resolve the inner region more finely. Full anisotropic number distributions determine the integrals; all coefficients use the same unit constant Coulomb logarithm.',
           '- The released 4,000 derived realizations per cluster are summarized separately; they are not aligned by row number to the nested chains or treated as observed BH counts.','',
           '## Reference results','',
           '| Cluster | Expected BH count | BH mass (solar masses) | Probe / BH mass | Central BH diffusion share | Share at 0.1 rh |',
           '|---|---:|---:|---:|---:|---:|']
    for r in flat:
        if not r['host_passed']:
            lines.append(f"| {NAMES[r['cluster']]} | Numerically unresolved | — | — | — | — |");continue
        ratio=r['probe_to_BH_reservoir_mass_ratio']
        lines.append(f"| {NAMES[r['cluster']]} | {r['BH_expected_number']:.3g} | {r['BH_mass_Msun']:.4g} | {ratio:.4g} | {100*r['BH_D_fraction_center_v0']:.3f}% | {100*r['BH_D_fraction_at_0p1rh_v0']:.4g}% |")
    lines += ['', 'Shares refer to the diffusion trace at zero remnant velocity, with equal species Coulomb-log normalization. They are local coefficients, not fractions of the time a remnant spends encountering BHs. Expected counts need not be integers; they do not specify a discrete count distribution.','',
              '![Reference and joint-row sensitivity profiles](../results/cluster-comparison/diffusion-comparison.png)','',
              'Blue: selected reference. Gray: accepted pilot rows. The connecting lines locate sampled values; they are not an adopted trajectory coefficient interpolator or a confidence band.','',
              '## Published model variation','',
              '| Cluster | Expected BH count: median [16th, 84th percentile] | Zero-BH derived rows |',
              '|---|---:|---:|']
    for r in published:
        q=r['BH_expected_number'];lines.append(f"| {NAMES[r['cluster']]} | {q['median']:.3g} [{q['p16']:.3g}, {q['p84']:.3g}] | {r['zero_BH_rows']}/{r['derived_rows']} |")
    lines += ['', 'These are empirical summaries of the released derived model draws, conditional on each standard model. In particular, the selected NGC 3201 reference is not representative of the zero-BH part of that distribution. Subsystem radii and probe/reservoir ratios remain undefined in absent-BH rows.','',
              '## Reconstructed joint-row sensitivity','',
              '| Cluster | Numerically accepted hosts | Accepted zero-speed profiles | Central BH share: pilot min–max | Zero-BH pilot hosts |','|---|---:|---:|---:|---:|']
    for c in NAMES:
        pilots=[r for r in allruns if r['cluster']==c and r['kind']=='pilot'];good=[r for r in pilots if r['host_status']=='passed'];local=[r for r in good if r.get('core_profile_status')=='passed']
        vals=[100*r['core_profile'][0]['BH_trace_D_fraction'] for r in local]
        limits=f'{min(vals):.3f}–{max(vals):.3f}%' if vals else 'Not resolved'
        zero=sum(r['BH_metrics']['BH_absent'] for r in local)
        lines.append(f'| {NAMES[c]} | {len(good)}/{len(pilots)} | {len(local)}/{len(pilots)} | {limits} | {zero} |')
    lines += ['', 'The min–max columns describe these eight attempts only. They are not posterior credible bounds, and a small pilot cannot establish a rare-event probability. Numerical checks are applied before a row contributes to local coefficient summaries.','',
              'One M4 pilot row (nested row 6903) failed to converge to a finite model extent under the compatible reconstruction. Its original weight, full parameters and solver error are retained. It contributes no local coefficient estimate, is not replaced, and is not counted as a zero-BH host. This unresolved numerical reconstruction is distinct from a conclusion that M4 observations are unusable.','',
              '## Interpretation and suitability','',
              'The comparison distinguishes reservoir size, spatial concentration and local diffusion importance. A large central BH share does not by itself demonstrate that continuous diffusion fails; nor does a small global BH count prove failure. Their consequences require a matched encounter and residence calculation. The 60-solar-mass probe has different ratios to each cluster and BH reservoir, so these controls do not hold those ratios fixed.','',
              'The comparison adopts six conditional spherical, nonrotating standard-fit backgrounds. NGC 3201 supplies a sparse-population reference under that host prescription. NGC 6397 uses its available standard-fit chain; its distinct no-BH alternative remains unevaluated because no corresponding chain is available in the inspected release. Rotation, flattening and population complexity lie outside the adopted Omega Centauri background. The joint-row pilot explores variation within each adopted family, while the [data-quality assessment](six-cluster-data-quality.md) retains source-specific applicability information.','',
              'The [environmental comparison](six-cluster-environment.md) repeats static boundary controls and observed-phase-point Galactic reference orbits separately. Its radius bounds concern neglected instantaneous external acceleration, not escape or cumulative trajectory error.','',
              '## Reproducibility','',
              '[Protocol](../analysis/cluster_comparison/protocol.json) · [Code and commands](../analysis/cluster_comparison/README.md) · [Reference table](../results/cluster-comparison/reference-comparison.csv) · [All attempts](../results/cluster-comparison/summary.json) · [Published-draw inventory](../results/cluster-comparison/posterior-inventory.json) · [Executed notebook](../analysis/cluster_comparison/cluster_comparison.ipynb).',
              '', 'Primary source products and their rights remain those in the [source manifest](../data/observed-clusters/sources.json). The comparison does not create new observational fits.']
    (ROOT/'docs/six-cluster-comparison.md').write_text('\n'.join(lines)+'\n')
    return flat
if __name__=='__main__':build()
