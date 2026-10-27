# Environmental controls for six cluster models

Each cluster uses its own phase-space row from the same pinned Baumgardt orbital table. The existing Irrgang Model I force law, reference-frame transformation and integration method are retained. These are Galactic orbits of cluster centers, not remnant trajectories.

## Common experiment

Four distance cases are evaluated where a numerically accepted host exists: catalog distance, reconstructed-host distance and the catalog distance plus or minus its quoted uncertainty. Distance changes both position and transverse velocity. The one-gigayear integration samples Galactic phases; it does not assume that the cluster structure survives unchanged for that time. No missing phase-space covariance is invented.

Static circular logarithmic controls use each cluster's own catalog Galactocentric radius and circular speeds of 220, 198 and 242 km/s. These are deliberately common sensitivity assumptions, distinct from the noncircular orbit in the adopted Galactic potential. The spherical host field continues as a point-mass field outside its distribution-function support, allowing roots there when appropriate.

The source-specific distant circular scenario formerly considered for NGC 3201 is not assigned to the other clusters.

## Conditional results

| Cluster | Numerical status | Static Jacobi scale (pc) | Host-distance peri–apo (kpc) | Refined minimum 1% force domain (pc) |
|---|---|---:|---:|---:|
| 47 Tuc | Passed | 131.38 | 5.461–7.487 | 12.569 |
| NGC 1851 | Passed | 159.65 | 0.525–20.392 | 1.017 |
| NGC 3201 | Passed | 84.13 | 8.215–23.810 | 8.173 |
| Omega Centauri | Passed | 179.48 | 1.360–7.060 | 8.748 |
| M4 | Passed | 55.39 | 0.650–6.415 | 1.031 |
| NGC 6397 (standard fit) | Passed | 55.29 | 2.644–6.331 | 3.477 |

The force domain is a sufficient instantaneous bound on the neglected external acceleration. It is not a physical cluster edge, an escape radius, a residence result or an accumulated trajectory-error bound. The static Jacobi scale and the eccentric-orbit frozen-frame diagnostic have different meanings.

## Resolution and suitability

Integration checks cover energy, axial angular momentum, independent step/tolerance refinement and time reversal. If a check fails, steps are halved and tolerances tightened, at most three times, while the original acceptance limits remain fixed. Preceding failed attempts remain in the case record. The minimum force-domain radius is separately checked under successively halved output cadence; all refinement levels are retained.

Agreement with catalog orbital extrema within 5% is a descriptive cross-check, not a numerical acceptance condition and not a new fit. Numerical convergence does not validate the chosen Galactic potential. The finite output mesh may still miss features narrower than the tested cadence. Host re-equilibration, stripping and the response of the BH population remain uncomputed.

M4 is the only cluster outside that 5% catalog comparison: the catalog-distance integration gives mean pericenter 0.772 kpc versus the published 0.64 kpc (20.6% difference), and apocenter 6.316 kpc versus 6.48 kpc (2.54%). Its orbit passes the numerical integration checks. The mismatch remains unresolved; the finite integration and rounded phase point do not reproduce the catalog averaging procedure. The host-distance case in the table is a different sensitivity case and must not be used to conceal this discrepancy.

The same source-quality restrictions apply as in the [cluster suitability assessment](six-cluster-data-quality.md). Missing host calculations permit catalog-only orbital results, not invented internal structure.

## Evidence

[Input rows and protocol](../data/cluster-comparison-environment/) · [Per-cluster orbit and static-control records](../results/cluster-comparison-environment/) · [Comparison CSV](../results/cluster-comparison-environment/comparison.csv) · [Reproduction](../analysis/cluster_comparison/README.md).

Primary sources: [Baumgardt orbital table](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/orbits.html), [Irrgang et al., Galactic potential models](https://arxiv.org/abs/1211.4353).
