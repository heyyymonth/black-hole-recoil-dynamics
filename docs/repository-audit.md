# Repository and calculation audit

Review date: 26 October 2026. Scope: the active data, host, collision, boundary and Galactic-orbit calculations, their recorded outputs, and the two-step manuscript. Earlier studies are preserved in [immutable reference records](computational-records.md).

## Findings and repairs

| Finding | Correction and verification |
|---|---|
| The default test command selected only historical experiments. | Default discovery and CI now run every active scientific test and the regression checks under `tests/`. |
| Main dependencies described the older Plummer package. | The main requirements now match the active numerical-check environment; full reconstruction retains its separate scientific lock. |
| Negative host radii were clamped in the potential but extrapolated in density. Nonfinite radii could propagate silently. | Host interfaces consistently reject invalid radii. Projected moments reject the zero-surface-density boundary where their ratio is undefined. |
| Configuration fields could be ignored or yield misleading labels. | Galactic potential normalization is read from its declared input. The named 1% diagnostic rejects other tolerances, and per-unit collision outputs require unit Coulomb logarithm. Unsupported runner methods, unit conversions and orientations fail explicitly. |
| Host and collision runners could report failed validation while returning success. | Failed checks now terminate with an error. Cached host reuse also requires the same protocol, not only the same selected parameter row. |
| Missing refinement results could be treated as zero error by the orbit verifier. | Required refinement fields must exist for the two refined cases; the proposed control radius is checked against its recorded construction. |
| Old experiments and unused figures/tables obscured the active work; the paper archive duplicated the 5.8 MB host array. | Historical implementations remain in immutable revisions. The active branch retains the current stages, and shared paper evidence is packaged once. |

## Evidence reviewed

The six-cluster data audit retains its original hashes and revision. All 51 locally cached original products were hash-checked. The chosen NGC 3201 row, all thirteen parameters, and normalized quadrature weight match the original posterior exactly. Independent checks recalculate all 321 HST component conversions, six catalog scale pairs, local tensor projections and species totals, and orbit-extremum summaries.

In the single-host audit for manuscript-v0.8.0.2, the host was revalidated with `run.py --reuse-solution`: this rechecks the existing numerical solution and its parameter/protocol identity, rather than repeating the external model solve. The twenty collision states, four static-boundary scenarios and four Galactic reference orbits were rerun. At that revision, all four result summaries and scientific output arrays/tables remained byte-identical to manuscript-v0.8.0.1. The cross-cluster repairs below subsequently changed local coefficient values at floating-point roundoff level and updated their provenance. Provenance records now bind the guarded implementations. The six-cluster data analysis was not rerun or relocked.

The formula review checked number/mass normalization in the scattering kernels, the diffusion factor of one half, time conversions, force/potential/tensor derivatives, rotating-frame conventions, and the analytic spatial force bound. The Galactic potential is compared with [Irrgang et al., Model I](https://arxiv.org/html/1211.4353v5); the host family and frozen-tide limitation retain [Gieles and Zocchi](https://arxiv.org/abs/1508.02120) and [Bertin and Varri](https://arxiv.org/abs/0808.2432) as their primary references. No demonstrated equation error or changed published value was found within these checks.

The eight evidence tables, eight figures, manuscript metadata, bibliography, PDF and portable archive are checked through their source manifests and the paper verifier. Markdown link and equation checks cover the active documentation rather than a fixed subset. Historical implementation details are outside the current scientific audit and are not certified by these checks.

## Interpretation retained

This is numerical and provenance verification, not external physical calibration. Smooth diffusion from a few modeled black holes, a tidally evolving host, full measurement/model covariance, remnant launch laws, impact-parameter cutoffs and physical escape criteria remain unresolved. The 8 pc control radius is conditional on one host, one reference potential, the stated force tolerance and numerically refined time sampling. It is not a measured cluster edge, escape radius, or accumulated orbit-error guarantee.

[Run the checks](../README.md#calculation-handoff) · [Manuscript verification](../paper/verification.txt) · [Current paper](../output/pdf/central-residence.pdf)

## Cross-cluster numerical checks

The six-cluster extension exposed two numerical issues beyond the original single-host checks. A positive log-cubic interpolator could reverse the nodal trend near the outer edge on an irregular mesh, creating an artificial density shell. The interpolator now checks the derivative polynomial throughout each interval and uses shape-preserving interpolation on intervals that violate the nodal trend. Both interpolants agree at their shared nodes. Two initially failed NGC 3201 pilot projections pass after revalidation; the original failure is reproduced on the same meshes in [numerical regression diagnostics](../results/cluster-comparison/numerical-regressions.json); they are not discarded as unsuitable observations.

The dense NGC 6397 standard-fit reference exposed cancellation error in a formally vanishing zero-speed drift. Antipodal quadrature directions are now exactly paired, and their contributions differ before summation. Even distributions retain zero drift exactly, while a streaming-distribution test retains its physical nonzero drift. Acceptance tolerances are unchanged. These repairs concern numerical integration and interpolation; they do not validate the diffusion closure.

The completed comparison includes six accepted references, 120 local states, 47 accepted pilot hosts out of 48 attempts, 24 numerically accepted Galactic orbits and 18 static controls. The M4 finite-extent solver failure and catalog-pericenter discrepancy remain explicit. Thirty-nine numerical and regression tests pass. A clean checkout can rebuild environmental calculations from the committed, hash-checked source rows without the ignored full-page cache. The 20-page PDF, executed notebook, source archive and current manifests have been checked together.
