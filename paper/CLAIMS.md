# Claim-to-evidence record

Manuscript revision 0.6.0-manuscript.1. Author: Venkata Hemanth Athota.

The paper defines a two-step transport and exposure problem and reports the completed six-cluster observation and published-posterior analysis. Data revision: `cluster-data-v0.1.0`. It also reports one compatible NGC 3201 host reconstruction and twenty evaluated local collision states. It reports no new observational fit, remnant trajectory ensemble or calibrated residence prediction.

| Statement | Basis | Boundary |
|---|---|---|
| Gravity and local scattering require a consistent mass-resolved host. | Species-number distribution definitions, mass moments, and Poisson's equation. | A potential and one velocity dispersion do not generally specify all collision coefficients. |
| Position, progenitor motion, recoil, and remnant mass belong in a joint launch law. | Initial-state kinematics and conditional probability. | No launch distribution or selection weights are inferred. |
| Drift and velocity diffusion have a declared Ito convention. | Standard kinetic theory and number-normalized Rosenbluth closure. | Local weak scattering requires calibration; strong encounters and host backreaction are unresolved. |
| The energy balance includes a diffusion contribution. | Ito's rule applied to specific kinetic plus host-potential energy. | It is not a conserved tidal escape criterion or a calibrated damping law. |
| Entry, re-entry, residence, and density exposure are distinct. | Operational region, stopping-time, and path-functional definitions. | Their distributions have not been evaluated. |
| Time-augmented occupation retains timing and velocity information. | A measure on time and phase space. | It is not a companion population or interaction law. |
| Ensemble exposure retains survival probability. | Conditional probability and interchange of nonnegative integrals. | Renormalizing survivors changes the quantity; unresolved continuation cannot be inferred. |
| The formulation supplies a calculation handoff. | Input register and acceptance requirements. | Candidate observed datasets are assessed; one host and local-kernel protocols are supplied; merger launches, physical cutoffs, trajectory solvers and calibration remain unassigned. |
| The six original Gaia field catalogs contain 621,878 rows; original HST products contain 107 radial bins. | Raw archive/file inspection; coverage and catalog-quality tables. | Field rows include nonmembers; three HST component records per bin do not represent independent measurements. |
| Twenty of 418 compiled kinematic rows require reconciliation. | Source-line flags and preserved original values. | Flagging does not establish that every value is physically impossible; no guessed corrections or approved compiled likelihood. |
| The intersection of HST median-radius spans is approximately 0.081–0.251 in projected R/rh. | Original median radii converted using one catalog's point scales. | Not exact annular overlap or a common three-dimensional domain; most bin boundaries are unavailable. |
| Mass–radius correlations are 0.900 in NGC 3201 and 0.944 in Omega Centauri. | Direct weighted covariance of original joint posterior rows; independently recomputed during manuscript review. | Conditional on published models/data; does not calibrate friction or establish exact source-fit reproduction. |
| Catalog and published-host summaries differ. | Three evidence tables and host-comparison figure generated from recorded values. | Overlapping observations and different assumptions prevent interpreting differences as evolution or independent model-rejection evidence. |
| A common interface can retain distinct conditional host distributions. | Native units, tracer metadata, posterior weights, covariance and model conventions. | Standardization does not make the clusters homologous or supply missing local species phase space. |
| Host structure and mass composition may alter normalized residence/exposure under controlled assumptions. | Proposed physical hypothesis motivated by structural differences and correlated host uncertainty. | Not tested by this data analysis; requires reconstructed host models and controlled transport calculations. |

| One joint NGC 3201 row produces mutually consistent gravity and species velocity distributions. | Host protocol, original input hashes, numerical arrays and independent quadrature checks. | Compatible historical implementation; exact original fitting versions not established; not a posterior-marginalized host. |
| Sixteen species bins reproduce mass, force, virial balance and DF moments within declared tolerances. | `results/host-reconstruction/summary.json`, density/velocity and projection checks, independent numerical tests. | Numerical consistency is not physical calibration; projection acceptance concerns the observed luminous tracer. |
| Twenty local drift and diffusion states pass quadrature and tensor checks. | `results/local-transport-kernels/`, analytic uniform-ball limits, rotational covariance and trace identity. | Per unit constant Coulomb logarithm at tested radii/speeds; no impact-parameter prescription or trajectory interpolation validated. |
| Expected 4.28 black holes supply 27.7% of central density but 95.8% of formal central diffusion at zero test velocity. | The same fixed host's species masses/densities and collision contributions. | Model expectations, not counted black holes or measured relaxation; sparse encounters require a physical validity test. |
| The formal truncation radius is 981 pc. | Selected distribution-function solution. | Not an observed Galactic tidal boundary or an adopted escape rule. |

Analytical consistency differs from numerical verification and physical validation. Raw-file checks and independent recomputation support the stated data calculations. Build and documentation checks assess source integrity and presentation. New observational measurements, publication-level novelty, dynamical calibration, companion-acquisition rates and merger probabilities are not claimed. The manuscript cites primary source studies, not this project's experiment records or software-tool references.

[Problem statement](../docs/problem.md) · [Calculation handoff](../docs/calculation-handoff.md)
