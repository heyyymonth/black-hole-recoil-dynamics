"""Describe all recorded environmental diagnostics, including unsupported cases."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/cluster-comparison-environment'
NAMES={'NGC0104':'47 Tuc','NGC1851':'NGC 1851','NGC3201':'NGC 3201','NGC5139':'Omega Centauri','NGC6121':'M4','NGC6397':'NGC 6397 (standard fit)'}
def main():
    rows=[];lines=['# Environmental controls for six cluster models','',
    'Each cluster uses its own phase-space row from the same pinned Baumgardt orbital table. The existing Irrgang Model I force law, reference-frame transformation and integration method are retained. These are Galactic orbits of cluster centers, not remnant trajectories.','',
    '## Common experiment','',
    'Four distance cases are evaluated where a numerically accepted host exists: catalog distance, reconstructed-host distance and the catalog distance plus or minus its quoted uncertainty. Distance changes both position and transverse velocity. The one-gigayear integration samples Galactic phases; it does not assume that the cluster structure survives unchanged for that time. No missing phase-space covariance is invented.','',
    'Static circular logarithmic controls use each cluster\'s own catalog Galactocentric radius and circular speeds of 220, 198 and 242 km/s. These are deliberately common sensitivity assumptions, distinct from the noncircular orbit in the adopted Galactic potential. The spherical host field continues as a point-mass field outside its distribution-function support, allowing roots there when appropriate.','',
    'The source-specific distant circular scenario formerly considered for NGC 3201 is not assigned to the other clusters.','',
    '## Conditional results','',
    '| Cluster | Numerical status | Static Jacobi scale (pc) | Host-distance peri–apo (kpc) | Refined minimum 1% force domain (pc) |',
    '|---|---|---:|---:|---:|']
    for c in NAMES:
        path=OUT/c/'summary.json'
        if not path.exists():lines.append(f'| {NAMES[c]} | Not executed | — | — | — |');continue
        d=json.loads(path.read_text());case=next((x for x in d['cases'] if x['case']=='host_distance'),None)
        static=next((x for x in d['static_controls'] if x.get('vc_kms')==220),{})
        r={'cluster':c,'all_numerical_checks_passed':d['all_numerical_checks_passed'],
           'static_rJ_pc':static.get('rJ_extended_pc'),
           'host_distance_pericenter_kpc':case.get('pericenter_mean_kpc') if case else None,
           'host_distance_apocenter_kpc':case.get('apocenter_mean_kpc') if case else None,
           'refined_minimum_1pct_radius_pc':case.get('refined_minimum_sufficient_radius_pc') if case else None,
           'catalog_orbit_within_5pct':d['catalog_orbit_comparison']['within_5_percent'],
           'number_numerically_accepted_orbits':sum(x['all_orbit_numerical_checks_passed'] for x in d['cases'])}
        rows.append(r)
        if r['static_rJ_pc'] is not None and r['refined_minimum_1pct_radius_pc'] is not None:
            lines.append(f"| {NAMES[c]} | {'Passed' if d['all_numerical_checks_passed'] else 'Unresolved checks'} | {r['static_rJ_pc']:.2f} | {r['host_distance_pericenter_kpc']:.3f}–{r['host_distance_apocenter_kpc']:.3f} | {r['refined_minimum_1pct_radius_pc']:.3f} |")
        else:lines.append(f"| {NAMES[c]} | {d['host_status']} | — | — | — |")
    lines += ['', 'The force domain is a sufficient instantaneous bound on the neglected external acceleration. It is not a physical cluster edge, an escape radius, a residence result or an accumulated trajectory-error bound. The static Jacobi scale and the eccentric-orbit frozen-frame diagnostic have different meanings.','',
              '## Resolution and suitability','',
              'Integration checks cover energy, axial angular momentum, independent step/tolerance refinement and time reversal. If a check fails, steps are halved and tolerances tightened, at most three times, while the original acceptance limits remain fixed. Preceding failed attempts remain in the case record. The minimum force-domain radius is separately checked under successively halved output cadence; all refinement levels are retained.','',
              'Agreement with catalog orbital extrema within 5% is a descriptive cross-check, not a numerical acceptance condition and not a new fit. Numerical convergence does not validate the chosen Galactic potential. The finite output mesh may still miss features narrower than the tested cadence. Host re-equilibration, stripping and the response of the BH population remain uncomputed.','',
              'M4 is the only cluster outside that 5% catalog comparison: the catalog-distance integration gives mean pericenter 0.772 kpc versus the published 0.64 kpc (20.6% difference), and apocenter 6.316 kpc versus 6.48 kpc (2.54%). Its orbit passes the numerical integration checks. The mismatch remains unresolved; the finite integration and rounded phase point do not reproduce the catalog averaging procedure. The host-distance case in the table is a different sensitivity case and must not be used to conceal this discrepancy.','',
              'The same source-quality restrictions apply as in the [cluster suitability assessment](six-cluster-data-quality.md). Missing host calculations permit catalog-only orbital results, not invented internal structure.','',
              '## Evidence','',
              '[Input rows and protocol](../data/cluster-comparison-environment/) · [Per-cluster orbit and static-control records](../results/cluster-comparison-environment/) · [Comparison CSV](../results/cluster-comparison-environment/comparison.csv) · [Reproduction](../analysis/cluster_comparison/README.md).','',
              'Primary sources: [Baumgardt orbital table](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/orbits.html), [Irrgang et al., Galactic potential models](https://arxiv.org/abs/1211.4353).']
    (ROOT/'docs/six-cluster-environment.md').write_text('\n'.join(lines)+'\n')
    with (OUT/'comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n');w.writeheader();w.writerows(rows)
    from environment import record_provenance
    record_provenance()
if __name__=='__main__':main()
