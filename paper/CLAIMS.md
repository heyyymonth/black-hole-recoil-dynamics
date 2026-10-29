# Claim-to-evidence record

Manuscript revision 0.11.0-manuscript.2. Author: Venkata Hemanth Athota.

The paper defines a two-step transport and exposure problem and reports the completed six-cluster observation and published-posterior analysis. Data revision: `cluster-data-v0.1.0`. A detailed NGC 3201 benchmark is extended to six conditional reference hosts and 120 local states, 48 joint-row pilot attempts, 18 common static controls and 24 Galactic reference orbits. It reports no new observational fit, remnant trajectory ensemble or calibrated residence prediction.

| Statement | Basis | Boundary |
|---|---|---|
| Gravity and local scattering require a consistent mass-resolved host. | Species-number distribution definitions, mass moments, and Poisson's equation. | A potential and one velocity dispersion do not generally specify all collision coefficients. |
| Position, progenitor motion, recoil, and remnant mass belong in a joint launch law. | Initial-state kinematics and conditional probability. | No launch distribution or selection weights are inferred. |
| Drift and velocity diffusion have a declared Ito convention. | Standard kinetic theory and number-normalized Rosenbluth closure. | Local weak scattering requires calibration; strong encounters and host backreaction are unresolved. |
| The energy balance includes a diffusion contribution. | Ito's rule applied to specific kinetic plus host-potential energy. | It is not a conserved tidal escape criterion or a calibrated damping law. |
| Entry, re-entry, residence, and density exposure are distinct. | Operational region, stopping-time, and path-functional definitions. | Their distributions have not been evaluated. |
| Time-augmented occupation retains timing and velocity information. | A measure on time and phase space. | It is not a companion population or interaction law. |
| Ensemble exposure retains survival probability. | Conditional probability and interchange of nonnegative integrals. | Renormalizing survivors changes the quantity; unresolved continuation cannot be inferred. |
| The formulation supplies a calculation handoff. | Input register and acceptance requirements. | Six reference hosts and local-kernel protocols are supplied. Finite-domain BH cutoffs are declared for the local comparison; merger launches, trajectory controls and physical calibration remain open. |
| The six original Gaia field catalogs contain 621,878 rows; original HST products contain 107 radial bins. | Raw archive/file inspection; coverage and catalog-quality tables. | Field rows include nonmembers; three HST component records per bin do not represent independent measurements. |
| Twenty of 418 compiled kinematic rows require reconciliation. | Source-line flags and preserved original values. | Flagging does not establish that every value is physically impossible; no guessed corrections or approved compiled likelihood. |
| The intersection of HST median-radius spans is approximately 0.081–0.251 in projected R/rh. | Original median radii converted using one catalog's point scales. | Not exact annular overlap or a common three-dimensional domain; most bin boundaries are unavailable. |
| Mass–radius correlations are 0.900 in NGC 3201 and 0.944 in Omega Centauri. | Direct weighted covariance of original joint posterior rows; independently recomputed during manuscript review. | Conditional on published models/data; does not calibrate friction or establish exact source-fit reproduction. |
| Catalog and published-host summaries differ. | Three evidence tables and host-comparison figure generated from recorded values. | Overlapping observations and different assumptions prevent interpreting differences as evolution or independent model-rejection evidence. |
| A common interface can retain distinct conditional host distributions. | Native units, tracer metadata, posterior weights, covariance and model conventions. | Standardization does not make the clusters homologous or supply missing local species phase space. |
| A diffusion-based response predictor may quantify the residence change from finite encounters in a declared regime. | Matched encounter generators, backward occupation equation and first-response construction. | Untested at the full jump endpoint; requires matched trajectory ensembles and uncertainty intervals. |
| One joint NGC 3201 row produces mutually consistent gravity and species velocity distributions. | Host protocol, original input hashes, numerical arrays and independent quadrature checks. | Compatible historical implementation; exact original fitting versions not established; not a posterior-marginalized host. |
| Sixteen species bins reproduce mass, force, virial balance and DF moments within declared tolerances. | `results/host-reconstruction/summary.json`, density/velocity and projection checks, independent numerical tests. | Numerical consistency is not physical calibration; projection acceptance concerns the observed luminous tracer. |
| Twenty local drift and diffusion states pass quadrature and tensor checks. | `results/local-transport-kernels/`, analytic uniform-ball limits, rotational covariance and trace identity. | Per unit constant Coulomb logarithm at tested radii/speeds; no impact-parameter prescription or trajectory interpolation validated. |
| Expected 4.28 black holes supply 27.7% of central density but 95.8% of formal central diffusion at zero test velocity. | The same fixed host's species masses/densities and collision contributions. | Model expectations, not counted black holes or measured relaxation; sparse encounters require a physical validity test. |
| The formal truncation radius is 981 pc. | Selected distribution-function solution. | Not an observed Galactic tidal boundary or an adopted escape rule. |
| A circular logarithmic-field benchmark gives an extended-host saddle at 84.13 pc. | `results/outer-boundary/`, extended enclosed-mass balance and exact circular saddle checks. | Prescribed field and one host; numerical proximity to a published estimate is not independent physical validation. |
| The static logarithmic baseline has a proposed 1% isolated-force limit at 20.67 pc. | Exact maximum differential acceleration divided by the same host gravity; tolerance sensitivities. | Chosen instantaneous error budget, not a universal boundary or accumulated trajectory-error guarantee. |
| The earlier 3rh coefficient radius has a maximum omitted Galactic force of 2.962%. | Baseline field applied to the existing local-state radius. | The integral remains defined; isolated motion there exceeds the chosen 1% force budget. |
| The host puts 5.63% of its mass beyond its conditional saddle. | Enclosed-mass profile at the solved saddle. | Not escaped mass; no tidal re-equilibration or stripping evolution is simulated. |
| Central circular Hill opening speed is 16.68 km/s. | Stationary conservative effective potential and cluster central potential. | Full remnant speed, not recoil alone; not conserved with noncircular tides or collisions and not an escape probability. |

| The catalog-reference orbit has mean pericenter 8.333 kpc and apocenter 25.362 kpc. | Original orbital-table row, Model I potential, registered integration and detected turning events. | Apocenter differs from the catalog by 1.77%; not exact source reconstruction or a new Galactic fit. |
| At the selected host distance, the sufficient 1% force domain varies approximately 8.2–39.8 pc. | Analytic Hessian bound applied along the same reference orbit, checked under cadence refinement. | A sufficient instantaneous bound, not a measured edge or accumulated orbit-error guarantee. The proposed fixed 8 pc control is conditional. |
| Neglected-force bounds reach 1.73% at rh and 29.5% at 3rh. | Time-dependent field with the unchanged reconstructed enclosed mass. | Bounds rather than measured accelerations; retaining exact tides removes this omitted-force error, not other host limitations. |
| The frozen-frame scale is 62.1–183.9 pc and the Euler/effective-eigenvalue ratio reaches 0.67. | Instantaneous tensor, radial-frame angular velocity and its derivative. | Not a stationary Jacobi surface, escape rule, or stripping calculation. |
| The forcing/dynamical-timescale indicator reaches 1.16 at rh. | Fixed-frame tensor derivative; independent norm-only indicator reaches 1.14. | Not a heating rate or proof of disruption; identifies an approximation requiring validation. |
| Galactic reference motion is calculated, but the host has not been tidally re-equilibrated. | Source/model scope and explicit host evolution equation. | No fitted tidal distribution, remnant trajectories, or calibrated residence predictions are supplied. |

Analytical consistency differs from numerical verification and physical validation. Raw-file checks and independent recomputation support the stated data calculations. Build and documentation checks assess source integrity and presentation. New observational measurements, publication-level novelty, dynamical calibration, companion-acquisition rates and merger probabilities are not claimed. The manuscript cites primary source studies, not this project's experiment records or software-tool references.

[Problem statement](../docs/problem.md) · [Calculation handoff](../docs/calculation-handoff.md)

## Cross-cluster calculation comparison

Six maximum-quadrature-weight joint references repeat the validated host and 20-state local coefficient calculation. Eight weighted joint draws per cluster provide an exploratory sensitivity pilot. Released 4,000-draw derived model products supply a separate empirical population summary. The comparison resolves spatial differences in the formal BH diffusion contribution; it neither establishes diffusion failure nor predicts residence. The adopted spherical nonrotating standard-fit families define the conditioning; NGC 6397 uses its available standard-fit chain, while its distinct no-BH alternative remains unevaluated.

Evidence: `results/cluster-comparison/provenance.json`, `summary.json`, `reference-comparison.csv`, `posterior-inventory.json`; per-cluster selections and numerical checks; `results/cluster-comparison-environment/summary.json`. The source archive retains compact comparison evidence and reproduction code; full reference host meshes are in the repository, while pilot meshes are reproducible scratch.

## Residence hypothesis and analytical bridge

The primary question is when a continuous encounter approximation preserves finite-time residence and species exposure at a declared accuracy. The paper formulates a finite-angle transition measure and a diffusion generator with matching first two moments. The backward occupation equation and stopped-generator identity connect their difference in mean residence to a nonlocal residual, subject to the stated generator-domain and integrability conditions. These are analytical constructions, not computed residence results or claimed new theorems.

A separate proposed predictor averages that residual under the diffusion process. It is the first response along an auxiliary interpolation of matched generators; its prediction at the full jump endpoint is untested and need not be accurate. The hypothesis can be tested against sign, magnitude and accuracy class of the measured residence change. Bounds for means do not establish distributional agreement.

The prescribed controls separate finite transitions, integer-population correlations and reciprocal response, and evolving-host/tidal effects. Resolved BH forces replace their smooth gravitational and collision contributions. A resolved two-body channel is not also used as a strong-encounter stop. Unresolved channels retain common censoring rules. Finite-domain encounter controls are now evaluated locally. No integer-population sampling law, trajectory ensemble, cluster residence residual or accuracy classification has been evaluated.

Primary grounding: Rosenbluth et al. (1957), Tep et al. (2022), Webb et al. (2018), and Dynkin (1956). Bar-Or et al. (2013) supplies related encounter-distribution context in a different massive-black-hole nuclear setting; no numerical result is transferred from that setting.

## Finite-angle local encounter law

The six reference hosts supply 162 accepted finite-domain BH scattering settings and 600,000 independent conditional increments. Rates, drift and raw second moments use the same velocity distribution, impact annulus and exact Newtonian deflection law. Conditional moments agree with sampled increments below 2.58 estimated standard errors in the registered samples. Density and moment refinement checks address numerical integration; they do not calibrate physical encounter cutoffs or independence assumptions.

The diagnostic fraction of diffusion trace with b90 greater than the local density scale is integrated on its relative-velocity domain. It is a scale diagnostic, not a proof of nonlocal evolution or diffusion failure. Excluded close-passage flux is recorded as a prospective common unresolved-event hazard; no such event has been applied to a remnant history.

The single-transition replacement algorithm is verified against exact finite-state derivatives, numerical differentiation and independent per-path random streams. Its application to a cluster residence observable and its continuous-time convergence remain uncomputed.

Evidence: `analysis/encounters/protocol.json`, `results/finite-angle-encounters/provenance.json`, per-host tensors and samples, `analysis/encounters/test_physics.py`, `analysis/encounters/test_response.py`.
