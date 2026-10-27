"""Create and execute the compact recorded-results companion notebook."""
from pathlib import Path
import sys
import nbformat as nbf
from nbclient import NotebookClient
from jupyter_client import KernelManager
ROOT=Path(__file__).resolve().parents[2]

def main():
    md=nbf.v4.new_markdown_cell;code=nbf.v4.new_code_cell
    cells=[md('''# Six-cluster comparison: recorded calculations

**Scope.** Conditional host reconstruction, local scattering kernels and Galactic environmental controls for the six shortlisted clusters. No remnant trajectories, residence distributions or diffusion-versus-encounter test have been performed.

**Readout.** Compare reservoir size, central and radial scattering importance, posterior sensitivity and numerical/data suitability. Selected reference rows are not typical hosts. Eight joint draws per cluster are a pilot, not a credible interval.

Sources: the hash-pinned Dickson/GCfit standard fits and Baumgardt orbital table in this repository. The six standard-fit backgrounds define the conditional comparison. NGC 6397 uses its available standard-fit chain; a distinct no-BH alternative has no corresponding chain in the inspected release and is not evaluated. Source-specific applicability information is retained in the data-quality assessment.

The cells read retained products and execute top-to-bottom. Rebuilding the expensive hosts is a separate command in the adjacent README.'''),
    code('''from pathlib import Path
import csv, json, sys
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, Markdown
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / 'results/cluster-comparison').is_dir())
DATA = ROOT / 'results/cluster-comparison'
load = lambda p: json.loads(p.read_text())
summary = load(DATA/'summary.json')
runs = summary['runs']
assert summary['attempted'] == summary['planned'] == 54
names = {'NGC0104':'47 Tuc','NGC1851':'NGC 1851','NGC3201':'NGC 3201','NGC5139':'Omega Centauri','NGC6121':'M4','NGC6397':'NGC 6397 (standard fit)'}
print(f"{sum(r['kind']=='reference' for r in runs)} reference hosts; {sum(r['kind']=='pilot' for r in runs)} pilot attempts")'''),
    md('## Reference models and local scattering\nAll six references use the same selection rule, 60-solar-mass diagnostic probe and unit constant Coulomb logarithm. Diffusion shares refer to zero remnant velocity, not orbital residence.'),
    code('''with (DATA/'reference-comparison.csv').open() as f:
    rows = list(csv.DictReader(f))
text = '| Cluster | Expected BH count | Central BH diffusion share | Share at 0.1 rh |\\n|---|---:|---:|---:|'
for r in rows:
    text += f"\\n| {names[r['cluster']]} | {float(r['BH_expected_number']):.3g} | {100*float(r['BH_D_fraction_center_v0']):.3f}% | {100*float(r['BH_D_fraction_at_0p1rh_v0']):.4g}% |"
display(Markdown(text))'''),
    md('## Radial variation across reconstructed joint rows\nBlue curves are reference rows; gray curves are accepted pilot realizations. Curves join sampled points and are not a confidence band or a validated trajectory interpolator.'),
    code('''fig, axes = plt.subplots(2,3,figsize=(11,6),sharex=True,sharey=True,layout='constrained')
for ax, (cluster, name) in zip(axes.flat,names.items()):
    selected = [r for r in runs if r['cluster']==cluster and r['host_status']=='passed' and r.get('core_profile_status')=='passed']
    for row in sorted(selected,key=lambda r:r['kind']=='reference'):
        profile = row['core_profile']; reference = row['kind']=='reference'
        ax.plot([p['r_over_rh'] for p in profile],[100*p['BH_trace_D_fraction'] for p in profile],color='#1b4f72' if reference else '#777777',alpha=1 if reference else .35,lw=1.8 if reference else .9)
    ax.set_title(name); ax.set_xscale('symlog',linthresh=.001); ax.set_ylim(-2,102); ax.set_xlim(0,3)
    ax.grid(axis='y',alpha=.2); ax.set_xlabel('r / half-mass radius')
for ax in axes[:,0]: ax.set_ylabel('BH diffusion trace share (%)')
plt.show()'''),
    md('## Released derived-model distribution\nThese 4,000 realizations per cluster are a separate published product. Their original generating seed and raw-row mapping are unavailable. No new observational fit or direct BH census is implied.'),
    code('''inventory = load(DATA/'posterior-inventory.json')['clusters']
text = '| Cluster | BH count median [16th,84th] | Zero-BH rows |\\n|---|---:|---:|'
for row in inventory:
    q=row['BH_expected_number']
    text += f"\\n| {names[row['cluster']]} | {q['median']:.3g} [{q['p16']:.3g}, {q['p84']:.3g}] | {row['zero_BH_rows']}/{row['derived_rows']} |"
display(Markdown(text))'''),
    md('## Numerical acceptance and missing information\nAn absent BH population is retained; undefined subsystem radii and reservoir ratios are not filled with zero. Failed numerical attempts remain inspectable. Source/model suitability is separate from numerical acceptance.'),
    code('''for cluster,name in names.items():
    pilot=[r for r in runs if r['cluster']==cluster and r['kind']=='pilot']
    accepted=sum(r['host_status']=='passed' for r in pilot)
    absent=sum(r.get('BH_metrics',{}).get('BH_absent',False) for r in pilot)
    print(f'{name}: {accepted}/{len(pilot)} pilot hosts accepted; {absent} absent-BH cases')
    for row in pilot:
        if row['host_status']!='passed': print('  unresolved:',row['label'],row.get('failure_log','')[-300:])'''),
    md('## Environmental controls\nFour distance cases per cluster and three static circular controls are conditional calculations. The force-domain bound is instantaneous and sufficient; it is not a remnant escape boundary.'),
    code('''with (ROOT/'results/cluster-comparison-environment/comparison.csv').open() as f: environment=list(csv.DictReader(f))
text='| Cluster | Static Jacobi scale (pc) | Refined 1% force domain (pc) | Numerical checks | Catalog extrema within 5% |\\n|---|---:|---:|---|---|'
for row in environment:
    text+=f"\\n| {names[row['cluster']]} | {float(row['static_rJ_pc']):.2f} | {float(row['refined_minimum_1pct_radius_pc']):.3f} | {row['all_numerical_checks_passed']} | {row['catalog_orbit_within_5pct']} |"
display(Markdown(text))'''),
    md('''## Interpretation

Population size and local diffusion importance are distinct. The radial profiles identify where a later encounter comparison should be resolved; they do not establish whether the diffusion approximation fails or changes residence. The pilot variation cannot be treated as converged posterior uncertainty. M4 has one unresolved pilot reconstruction; its catalog-distance orbit passes numerical checks but differs by 20.6% from the published mean pericenter. Neither discrepancy is concealed by selecting another case.

The adopted backgrounds are spherical, nonrotating standard fits. NGC 3201 provides a conditional sparse-population reference; NGC 6397 uses the available standard-fit chain, with its distinct no-BH alternative left unevaluated. Rotation, binary populations, discrete counts, backreaction and physical scattering cutoffs require separate specifications for a matched encounter experiment. The data-quality assessment records source-specific applicability information alongside these experimental assumptions.

See `docs/six-cluster-comparison.md`, `docs/six-cluster-data-quality.md`, and `docs/six-cluster-environment.md`. Source URLs and exact hashes are retained in `data/observed-clusters/sources.json`, the per-cluster selection records, and the comparison provenance manifest.''')]
    nb=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'name':'python3','display_name':'Python 3','language':'python'}})
    manager=KernelManager(kernel_name='python3');manager.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
    NotebookClient(nb,km=manager,timeout=180,resources={'metadata':{'path':str(ROOT)}},record_timing=False).execute()
    nbf.validate(nb);path=Path(__file__).with_name('cluster_comparison.ipynb');nbf.write(nb,path)
    print('Executed',path.name)
if __name__=='__main__':main()
