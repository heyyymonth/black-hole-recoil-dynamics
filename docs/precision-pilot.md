# Precision of residence-response estimators

The question is whether changing how the same response is sampled can reduce its uncertainty at useful computational cost. The [completed response test](response-test.md) estimates a finite-step derivative and compares it with independent mixture trajectories. Its prediction precision is unresolved. This pilot examines the estimator before extending that physical comparison.

The [fixed protocol](../analysis/precision/protocol.json) defines three methods, three launch directions and four independent blocks per method and launch. It holds the conditional Omega Centauri host, encounter table, 60-solar-mass remnant, BH collision law, all-species gravity, residence region, horizon and stopping rules fixed. The earlier response protocol, source files and retained outcomes remain unchanged.

## One response target, three sampling methods

For the complete numerical transition mixture `P_a = (1-a) P_D + a P_J`, write `mu_h(a) = E_a[C_T/T]` and `S_h = mu_h'(0)`. Each transition includes gravity, occupation accumulated during the step, collisions and stopping. All three methods estimate this same finite-step `S_h`.

| Method | Gaussian factor used within a quartet | Replacement-step selection |
|---|---|---|
| `eigen_uniform` | Existing eigenbasis factor | Uniform over all 1,000 scheduled steps |
| `symmetric_uniform` | Symmetric positive-semidefinite covariance root | Uniform over all 1,000 scheduled steps |
| `symmetric_stratified` | Symmetric positive-semidefinite covariance root | Fixed proportional samples from four time strata |

At an active replacement, the encounter branch is split into a zero-event branch and a positive-count branch with their correct Poisson probabilities. Two diffusion branches use opposite Gaussian draws at that transition with unchanged drift. All branches then continue under diffusion. With full-horizon normalized rewards `Y_0`, `Y_+`, `Y_D+` and `Y_D-`, the score remains

$$
Z=N\left[p_0Y_0+p_+Y_+-\frac{Y_{D+}+Y_{D-}}{2}\right].
$$

The zero/positive probabilities and every branch's residence and termination record are retained. A path stopped before its selected insertion contributes zero and remains in the denominator. No scores are clipped or removed.

## Covariance factor and path coupling

For a diffusion tensor `D = V Lambda V^T`, the existing factor is `F = V sqrt(Lambda)`. The symmetric factor is `G = V sqrt(Lambda) V^T`. Both satisfy `F F^T = G G^T = D`. Consequently, applying either factor to a standard normal vector gives the same marginal Gaussian covariance. Drift, duration and all other operations are unchanged.

The two factors need not produce the same realized path from the same normal vector. They define different joint couplings between branches. Arbitrary eigenvector signs and basis changes can weaken the eigenbasis coupling between nearby states, whereas the unique symmetric PSD root is invariant under those representation choices. Its continuity does not imply a uniform Lipschitz bound near zero eigenvalues, nor does it make it the optimal coupling for an entire residence history. Matrix-factor and Gaussian-coupling results provide the mathematical basis for distinguishing the marginal law from its pairing. [Bhatia, Jain and Lim, *On the Bures-Wasserstein distance between positive definite matrices*](https://arxiv.org/abs/1712.01504)

The implementation loads the published transport source into an isolated module namespace and changes only that instance's Gaussian-increment function. The eigen arm retains the original update exactly. Numerical controls check full-trajectory equivalence for the eigen arm, covariance and drift preservation, repeated-eigenvalue basis invariance, stopping rewards and isolation of the published module. The same roundoff tolerance for small negative eigenvalues is used in both factors.

## Proportional insertion-time stratification

The four strata contain scheduled steps `[0,250)`, `[250,500)`, `[500,750)` and `[750,1000)`. Each has probability `w_s = 1/4` under uniform insertion. Every stratified block samples 256 independent uniform integer insertions within each stratum, for a total of 1,024 quartets. Insertion randomness has a stream separate from transition randomness.

Let `mu_s = E[Z | s]`, `sigma_s^2 = Var(Z | s)`, and `n_s` be the sample count in stratum `s`. Then

$$
S_h=\sum_s w_s\mu_s,\qquad
\widehat S_h=\sum_s w_s\overline Z_s,\qquad
\operatorname{Var}(\widehat S_h)=\sum_s\frac{w_s^2\sigma_s^2}{n_s}.
$$

For proportional allocation `n_s = n w_s`, the law of total variance gives

$$
\operatorname{Var}(\widehat S_{\rm uniform})-
\operatorname{Var}(\widehat S_{\rm stratified})
=\frac{1}{n}\sum_s w_s(\mu_s-S_h)^2\geq0.
$$

This population-variance statement compares the same symmetric-factor score law under uniform and proportionally stratified insertion. It does not compare either one with the eigenbasis score law. Estimated variances can fluctuate, and the identity does not guarantee a practically important gain or lower computing cost. No allocation is adapted to the observed score means, signs or tails.

For the stratified method, uncertainty is calculated from the within-stratum variances and fixed weights. Treating the pooled scores as independent and identically distributed would include between-stratum variation that the sampling design has already controlled.

## Fixed calculation and independence

The primary encounter table and step `0.01 t_dyn` are retained, with horizon `10 t_dyn`, residence region `r < r_rho`, and the existing spatial and velocity validity boundaries. The radial, tangential and 45-degree oblique launches retain the same initial radius and speed.

The design specifies **36 jobs**, each with 1,024 quartets: **36,864 independent quartets comprising 147,456 branch trajectories**. Four branches within a quartet are correlated. Different methods, launches and blocks have independent streams; there is no claimed paired covariance reduction between methods. Launch and method execution orders rotate between blocks to reduce systematic association with cache warm-up and machine conditions.

Seeds and the complete job inventory are fixed in the execution design before the first calculation. The lossless sampler cache is reused from the frozen response experiment, with the same physical-input namespace. Cache data are reconstructible scratch; retained arrays, metadata, protocol and source hashes constitute the calculation record.

## Reporting and interpretation

The primary precision measure is `n Var(estimated mean)`, called the equivalent score variance. It compares the sampling designs at their equal sample counts. Per-block and pooled means and standard errors are descriptive. The six independent-arm mean differences also have descriptive standard errors; they are not confidence-based equivalence tests. The common expectation target follows from the construction and numerical controls, not from observing compatible sample means.

There are six efficiency comparisons: eigen/uniform versus symmetric/uniform, and symmetric/uniform versus symmetric/stratified, for each of the three launches. Ratios are oriented as the first method's mean-estimator variance divided by the second's; values above one favor the second method. Empirical percentile bootstrap intervals use 8,192 replicates with a nominal Bonferroni family size of six. Resampling preserves each method's insertion design, including separate resampling within each fixed stratum. These are approximate intervals and do not protect against unobserved rare contributions. Invalid zero-variance denominators are reported rather than discarded.

Equal quartet counts do not imply equal CPU budgets. The cost-adjusted measure is mean-estimator variance multiplied by recorded process CPU time. CPU time includes trajectory calculation and cache access, but excludes setup and postprocessing. Timing remains conditional on the machine, cache and competing load; it is held fixed during resampling, so the cost-adjusted intervals do not include timing uncertainty. Wall time and cache loads/builds are retained separately.

Diagnostics include score extrema and tails, signed insertion-time contributions, branch-reward covariance, stopping causes and numerical continuation allowances. Branch covariance is descriptive; the response score also has state-dependent Poisson weights. Better precision does not remove numerical-domain continuation uncertainty or establish grid/time convergence.

A method supported by this pilot can be selected for a subsequent independently sampled response test. Earlier endpoint validation is not reused to certify the selected estimator's predictive accuracy. No conclusion about capture, companion acquisition, merger rates or calibrated real-cluster dynamics follows from an estimator-efficiency comparison.

[Implementation and reproduction](../analysis/precision/README.md) · [Residence-response hypothesis](transport-hypothesis.md) · [Completed response benchmark](response-test.md)

## Completed findings

All 36 planned ensembles completed and passed independent record checks. The pre-execution protocol and implementation were published in commit `7e3b32d`; every retained sample belongs to the fixed design. Ratios below are estimated variance of the first method divided by that of the second. Brackets are the prespecified exploratory six-family bootstrap intervals.

| Launch | E/S variance ratio | S/T variance ratio | E/S including recorded CPU cost |
|---|---:|---:|---:|
| Radial | 112.50 [33.35, 517.46] | 0.876 [0.129, 6.379] | 110.04 |
| Tangential | 1133.58 [639.71, 2367.00] | 1.099 [0.410, 3.017] | 1183.56 |
| Oblique | 1885.25 [1073.74, 3855.74] | 0.978 [0.409, 2.517] | 1928.84 |

All three E/S intervals lie above one. At equal sample count, the standard-error reduction factors are 10.61, 33.67 and 43.42. All twelve individual batch E/S point estimates also favor S, ranging from 45.42 to 3256.83; that batch diagnostic does not carry an additional confidence claim. CPU-adjusted ratios retain the large gain under the measured cache and machine conditions.

All three S/T intervals include one. There is no resolved additional benefit from the fixed four insertion-time strata. Estimated ratios below one reflect sampling variation and do not contradict the population identity. The result supports symmetric/uniform as a straightforward reference for another independently sampled response test; it does not establish uniform sampling as superior to proportional stratification.

| Launch | Symmetric/uniform response mean | Descriptive standard error | Finite-sample continuation allowance |
|---|---:|---:|---:|
| Radial | −0.021604 | 0.018500 | 0 |
| Tangential | 0.005397 | 0.019142 | 0.00003849 |
| Oblique | −0.011634 | 0.013035 | 0 |

These errors still exceed the earlier 0.0005 simultaneous-half-width target before any multiplicity adjustment. The means do not resolve the response sign or endpoint accuracy. Zero observed continuation allowance does not bound its population contribution. The tangential stratified arm has an allowance of 0.0002162, above the earlier 0.0002 sample budget. A precision gain does not remove continuation or refinement work.

The prespecified diagnostic shows that the largest one percent of absolute scores carries 33–41% of total absolute score among the symmetric-factor ensembles. A further descriptive diagnostic, specified after the pilot, examines the centered-square contributions to estimated variance. For symmetric/uniform, one radial score supplies 61.1% of that sample variance; the largest 41 contributions (approximately one percent) supply 92.9%, 89.6% and 84.0% in radial, tangential and oblique order. These are concentrations in the observed sample, not population-tail estimates. No observations are removed, and neither the fixed analysis nor its confidence family is changed. This is why the efficiency ratios cannot guarantee future sample sizes or coverage against unseen rare scores.

The experimental conclusion is a large, reproducible improvement in estimation for the same finite-step model. It is not confirmation of the response hypothesis, continuous-time convergence, or calibrated cluster dynamics. Independent response predictions and validation with explicit refinement and continuation controls remain necessary.

[Recorded summary](../results/precision/summary.json) · [Comparison figure](../results/precision/comparison.png) · [Post-pilot variance diagnostic](../results/precision/tail-diagnostic.json)

The subsequent [independent response test](independent-response-test.md) uses this selected estimator with fresh predictions and validation; it keeps the pilot samples separate.
