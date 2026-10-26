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

The host was revalidated with `run.py --reuse-solution`: this rechecks the existing numerical solution and its parameter/protocol identity, rather than repeating the external model solve. The twenty collision states, four static-boundary scenarios and four Galactic reference orbits were rerun. All four result summaries and the scientific output arrays/tables remain byte-identical to manuscript-v0.8.0.1. Provenance records now bind the guarded implementations. The six-cluster data analysis was not rerun or relocked.

The formula review checked number/mass normalization in the scattering kernels, the diffusion factor of one half, time conversions, force/potential/tensor derivatives, rotating-frame conventions, and the analytic spatial force bound. The Galactic potential is compared with [Irrgang et al., Model I](https://arxiv.org/html/1211.4353v5); the host family and frozen-tide limitation retain [Gieles and Zocchi](https://arxiv.org/abs/1508.02120) and [Bertin and Varri](https://arxiv.org/abs/0808.2432) as their primary references. No demonstrated equation error or changed published value was found within these checks.

The seven evidence tables, seven figures, manuscript metadata, bibliography, PDF and portable archive are checked through their source manifests and the paper verifier. Markdown link and equation checks cover the active documentation rather than a fixed subset. Historical implementation details are outside the current scientific audit and are not certified by these checks.

## Interpretation retained

This is numerical and provenance verification, not external physical calibration. Smooth diffusion from a few modeled black holes, a tidally evolving host, full measurement/model covariance, remnant launch laws, impact-parameter cutoffs and physical escape criteria remain unresolved. The 8 pc control radius is conditional on one host, one reference potential, the stated force tolerance and numerically refined time sampling. It is not a measured cluster edge, escape radius, or accumulated orbit-error guarantee.

[Run the checks](../README.md#calculation-handoff) · [Manuscript verification](../paper/verification.txt) · [Current paper](../output/pdf/central-residence.pdf)
