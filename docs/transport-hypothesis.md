# Encounter physics and the accuracy of remnant residence

The research question is when a continuous drift-and-diffusion description preserves the finite-time residence and environmental exposure of a kicked black-hole remnant. The published cluster distributions are the adopted calculation settings. The six reference hosts, their joint-row sensitivity calculations, local kernels and Galactic diagnostics supply the environments for this question.

## Physical hypothesis

A finite velocity change matters for residence when it changes the remnant's subsequent occupation of the chosen region. Its importance therefore depends on the encounter distribution, the remnant's position and velocity, and the remaining time. A central diffusion fraction or global expected BH count alone does not specify this effect.

A testable prediction is that a response estimate computed from the diffusion process can predict the sign, size or accuracy class of the residence change produced by finite encounters in a declared regime. The predicted response is evaluated before comparison with full encounter ensembles. It may fail when repeated finite transitions produce a large nonlinear change in the paths.

## From the adopted background to an encounter law

Each complete host realization supplies the same species number distributions to the potential, encounter calculation and exposure integral. For local independent encounters with single objects, the differential rate is

$$
d\Gamma_j=F_j(\mathbf r,\mathbf u,t)|\mathbf v-\mathbf u|\,d^3u\,b\,db\,d\varphi.
$$

The physical domain specifies impact parameters, closest approaches and excluded processes. Exact elastic Newtonian two-body scattering supplies a velocity change for each allowed encounter. Pushing the rate through that map defines the jump measure. Its first two velocity-change moments then define the matched diffusion process.

The existing per-unit-Coulomb-log kernels provide weak-scattering calculations and numerical checks. A selected finite-angle domain requires its own matched moments; multiplying the existing outputs by an arbitrary cutoff-dependent number does not establish that match. If large deflections are resolved separately, their impact-parameter range is excluded from the diffusion contribution.

The [finite-angle calculation](finite-angle-encounters.md) now specifies that local domain and evaluates the rates, moments and conditional elastic increments across the six reference hosts. It also implements and tests a single-transition replacement estimator in finite-state examples. Numerical acceptance remains recorded per local state; the matched collisional residence and physical response prediction have not been evaluated. The [collisionless controls](collisionless-residence.md) separately establish timed reference paths, finite-region occupation and path-integrated encounter budgets.

## The analytical bridge

Let `U_D(t, r, v)` be the expected remaining occupation under diffusion. Its backward equation has the occupation weight as a source and zero terminal and killed-state payoff. For residence, the weight is membership in the declared spatial region; for species exposure it is the species density.

With a common initial law and stopping convention, the difference in mean occupation between jump and diffusion processes equals the jump-process expectation of the time-integrated generator residual acting on `U_D`. That residual is the full velocity change in remaining occupation after subtracting its linear and quadratic velocity terms.

This identity is standard stopped-generator calculus applied to the observable. It requires integrability and a suitable generator domain. A smooth region weight permits an initial derivation; the sharp-region limit and absorbing boundaries require justification. Boundary-crossing jumps use the full zero-continued nonlocal expression.

The identity is not a numerical result or a hypothesis needing empirical confirmation. The separate **response predictor** replaces its jump-process expectation with a diffusion-process expectation. This is the first response along a diagnostic interpolation of the two generators. Its extrapolation to the full encounter process is the proposition tested by the experiment. The manuscript states the equations and the corresponding sufficient bound for the mean when the uniform residual bound is finite.

## Controlled comparisons

| Comparison | Held fixed | Physical information varied | Required outcome |
|---|---|---|---|
| Finite transitions versus diffusion | Host, launch, encounter domain, first two moments, external force, region, horizon and stops | Full encounter increments versus their second-order approximation | Mean and distributional residence differences; accuracy of the response predictor |
| Independent encounters versus a finite population | Conditional host realization and declared initial population ensemble | Persistent identities, repeated encounters, reciprocal BH forces and reservoir response | Change attributable to correlations and backreaction |
| Fixed versus evolving host | Launch, encounter prescription and external forcing | Consistent evolution of the species distributions and potential | Change attributable to the environment's response |
| Across adopted hosts | Declared dimensionless launch, region and time controls | Host structure and mass composition | Range of conditional regimes; no automatic universal cluster weighting |

A finite population requires an explicit integer-count and mass-sampling law; an expected count is not rounded silently. Resolved BH forces replace the smooth BH gravitational contribution, and BH collision terms already represented by those forces are removed. The remnant must act back on the resolved objects. Binaries and relativistic capture require additional interaction models.

The original fixed 60-solar-mass probe does not match the remnant-to-host or remnant-to-individual-scatterer mass ratios across clusters. Those ratios, radial concentration, encounter correlation time and relative-velocity-dependent deflection strength are separate controls.

## What the current findings contribute

- The reconstructed distributions connect gravity, scattering and exposure without mixing different host realizations.
- The radial diffusion profiles identify where launch and crossing conditions should be resolved. Compact and more extended contributions give distinct test settings.
- Joint-row variation supplies conditional parameter sensitivity. Eight draws per cluster remain a pilot.
- Zero-BH realizations provide a null comparison when only BH encounters are varied.
- Galactic-force bounds define isolated controls. A force-accuracy boundary remains distinct from an escape surface.
- Recorded unresolved calculations remain explicit; no result is filled by a replacement host or an assumed residence history.

## Criteria for interpretation

Compare normalized mean residence, its cumulative distribution, species exposure, entry times and termination probabilities. Retain censored histories and their causes. If two-body strong encounters are resolved, remove the corresponding strong-encounter stop from both comparison arms; otherwise neither arm can measure their subsequent effect.

Predeclare accuracy tolerances and reduce numerical, ensemble and smoothing errors below them. An uncertainty interval for a nonnegative discrepancy entirely below tolerance supports adequacy in that setting. An interval entirely above it establishes a material difference. An interval spanning it is inconclusive. A small mean difference alone does not establish agreement of tails or termination outcomes.

The intended theoretical result is an observable-specific accuracy criterion tested across the adopted environments. The collisionless reference paths establish the zero-collision control. No matched encounter ensemble, collisional residence distribution or evaluated residual has been produced.

## Primary foundations

[Rosenbluth et al.](https://doi.org/10.1103/PhysRev.107.1) supplies the inverse-square collision-moment framework; [Tep et al.](https://arxiv.org/abs/2201.03985) treats anisotropic cluster relaxation; [Webb et al.](https://arxiv.org/abs/1711.09100) studies kicked black holes; [Dynkin](https://doi.org/10.1137/1101004) supplies the generator framework. [Bar-Or et al.](https://arxiv.org/abs/1209.4594) studies encounter distributions in a different massive-black-hole nuclear setting; its numerical conclusions are not transferred to these hosts.

[Problem statement](problem.md) · [Completed comparison](six-cluster-comparison.md) · [Manuscript](../paper/manuscript.tex)
