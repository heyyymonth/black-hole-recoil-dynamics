# Testing the residence-response prediction

The question is whether the change in mean residence caused by finite encounters can be predicted from one encounter-transition replacement in otherwise diffusive histories. The earlier [matched comparison](matched-residence.md) measures the two numerical endpoints directly. This experiment evaluates the first response independently, then compares its prediction with new mixture trajectories. It uses the same conditional Omega Centauri host, BH encounter law, 60-solar-mass remnant, all-species smooth gravity, region and stopping conventions.

The frozen protocol is [analysis/response/protocol.json](../analysis/response/protocol.json). Its tolerances and sampling limits are inputs to the test, not quantities selected after seeing validation results.

## Quantity tested

Let `Y = C/T` be residence divided by the full horizon, including stopped paths. For the complete transition mixture `P_a = (1-a) P_D + a P_J`, define `mu_h(a) = E[Y]` and `S_h = d mu_h(a)/da` at `a=0`. The index `h` retains the numerical table and step size.

The prediction is `mu_h(a) - mu_h(0) ≈ a S_h`. The residual is `R_h(a) = mu_h(a) - mu_h(0) - a S_h`. Differentiating the finite product of transition operators proves the single-replacement identity for `S_h`; it does not prove that the first derivative predicts the full endpoint `a=1`.

Every transition includes both gravity half-flows, occupation accumulated during the step, the collision update, and all stopping outcomes. Replacing only a velocity increment after otherwise completing the step would test a different quantity.

## Response estimator

Choose one scheduled step uniformly from all `N` steps, independently for each sample. Follow diffusion to that step. A sample already stopped contributes zero and remains in the denominator.

At the selected transition, condition on its excluded-event clock. With the resulting duration `delta`, the allowed count has mean `lambda = Gamma * delta`. Split its law into zero events with weight `p0 = exp(-lambda)` and a positive Poisson count with weight `pplus = 1-p0`. The positive branch permits more than one event. Two diffusion branches use opposite normal increments with the same drift. After this transition, all four branches continue with diffusion.

The sample score is

$$
\widehat S_i=N\left[p_0Y_0+p_+Y_+-\frac{Y_{D+}+Y_{D-}}{2}\right].
$$

Aligned normal and exponential random numbers couple the branches; each separated branch still evaluates its own local coefficients and hazard. The score has the same expectation as the direct two-branch replacement estimator. Conditioning and antithetic normal increments aim to reduce its variance, but their efficiency must be measured.

An independent pilot compares 1024 samples of each estimator for each launch. The weighted-score variance alone determines the primary production count, between 8192 and 65,536 samples according to the fixed rule. Pilot samples do not enter the production prediction. A capped count does not establish the requested precision.

## Comparisons and controls

| Control | Registered choice | Purpose |
|---|---|---|
| Launch | Radial outward, tangential, and 45-degree oblique; same radius and speed | Separate orbital direction from initial energy |
| Primary resolution | Refined table; step `0.01 t_dyn` | Main prediction and mixture comparison |
| Numerical controls | Coarse and refined tables; step `0.02 t_dyn` | Separate table and step sensitivity |
| Horizon and region | `T = 10 t_dyn`; `r < r_rho` | Preserve the earlier observable |
| Primary mixture values | `0.05, 0.1, 0.25, 0.5, 1` | Examine small changes and the full endpoint |
| Control mixture values | `1` | Check endpoint sensitivity |
| Validation sampling | 8192 mixture/control pairs per setting | Estimate the shift using paired differences |

All prediction arrays and their hash register must be saved before any validation trajectories are generated. Validation uses independent streams, with a separate paired diffusion control for every mixture value. The radial and tangential endpoint conditions were seen in the earlier pilot, so fresh samples there are replication. The oblique launch is a previously unused direction; it does not provide a new host or physical calibration.

## Decision and uncertainty

At `a=1`, a residual interval wholly within `[-0.001, 0.001]` supports the specified finite-resolution prediction accuracy. An interval wholly outside establishes disagreement. Overlap is inconclusive. At `a=0.05` and `0.1`, report `R_h(a)/a`; these are finite-change checks, not direct observations of a derivative limit.

The normal intervals use a predeclared Bonferroni family size of 64. The 63 comparisons consist of nine response means, 21 validation shifts, 21 residuals, six response-refinement contrasts and six endpoint-refinement contrasts. The same response estimate may appear in several predictions; independence among those reported errors is not required by Bonferroni. Independence between the prediction and validation batches is required by the variance calculation. Coverage remains approximate, especially when rare scores dominate.

The numerical-sensitivity budget is `0.0005`; the numerical missing-residence budget is `0.0002`. A numerical-domain stop permits up to `(T-stop_time)/T` omitted residence. For the four-branch response, its pathwise allowance is the sum of these branch bounds multiplied by the absolute score weights, including `N`. Excluded close-passage stops are part of the modeled killed process and are kept separate from numerical-domain continuation uncertainty.

Signed score contributions are grouped by the predeclared insertion-time and radius bins, using all samples in each denominator. These decompose the sampled score along the diffusion history; unresolved response estimates do not identify physical sensitivity without a separate uncertainty analysis. Selecting only survivors or conditioning on large contributions would change the quantity.

Finite grid/time contrasts do not supply rigorous continuum bounds. Small sampling error alone does not establish physical accuracy. The result concerns a mean-residence predictor; it does not establish agreement of residence distributions, capture probabilities, companion-acquisition rates or merger probabilities.

## Completed observations

All six pilot, nine prediction and 21 validation ensembles are complete. The retained evidence contains 423,936 independent pairs or quartets and 1,345,536 branch trajectories; branches within each group are correlated. The prediction register predates every validation run. Pilot samples are excluded from the production estimates.

At the primary resolution, each response estimate uses 65,536 weighted samples and each endpoint shift uses 8192 independent validation pairs. The table reports estimates and approximate family-adjusted intervals for normalized residence. The error is the measured endpoint shift minus the response prediction.

| Launch | Response slope | Endpoint shift | Prediction error | Accuracy test |
|---|---|---|---|---|
| Radial | -0.03680 [-0.19323, +0.11962] | -0.00782 [-0.01241, -0.00323] | +0.02899 [-0.12751, +0.18548] | Inconclusive |
| Tangential | -0.14211 [-0.71474, +0.43052] | -0.00671 [-0.01219, -0.00122] | +0.13540 [-0.43725, +0.70806] | Inconclusive |
| Oblique | -0.06169 [-0.61073, +0.48735] | -0.00432 [-0.00919, +0.00056] | +0.05737 [-0.49169, +0.60644] | Inconclusive |

The radial and tangential endpoint shifts exclude zero in these intervals; the oblique interval includes zero. These statements concern the fixed numerical process. All response intervals are much wider than the measured endpoint shifts and the requested precision. Their half-widths are 0.15643, 0.57263 and 0.54904, compared with the target 0.0005. All three endpoint accuracy tests are therefore inconclusive: the results establish neither agreement nor failure of the response hypothesis.

All six small-mixture scaled-residual intervals include zero. That does not demonstrate a derivative limit. None of the twelve grid/time refinement contrasts has its whole interval within the 0.0005 budget. Primary response continuation allowances are 0.11345, 0.05044 and 0.04497, exceeding the 0.0002 budget. These conservative finite-sample bounds can be amplified by the response weight and are not measured physical missing residence or population confidence bounds.

The pilot's weighted/direct variance reduction factors are 1.96, 1.27 and 1.21. The variance-only count rule reached its cap for every primary launch. This makes estimation precision the immediate obstacle; repeating the same method without addressing its variance would require substantially greater sampling.

## A specific precision question

A deterministic diagnostic, specified after examining the pilot, compares covariance factors at neighboring states. Eigenbasis factors and symmetric covariance roots both reproduce the correct marginal covariance. Eigenvector signs and basis changes can nevertheless make the eigenbasis factor discontinuous under shared noise. The retained radial-state example has a relative tensor change of about 4.13e-10, yet normalized shared-noise distance squared of 1.26; the symmetric root gives 4.33e-20 for the same local comparison.

This identifies a controlled follow-up: compare trajectory-estimator variance under the two covariance factors while keeping the marginal law, cutoffs, rewards and stops fixed. The separate [fixed precision comparison](precision-pilot.md) now evaluates that trajectory-estimator effect. The local diagnostic alone does not establish the cause of aggregate variance or guarantee the precision target. The completed experiment retains its original factor and protocol. Its lossless disk storage adapter changes only sampler storage and records equivalence tests and provenance.

The result leaves the physical response hypothesis open. Resolving the estimator, numerical refinement and boundary continuation is required before extending an accuracy claim to other hosts or to calibrated cluster dynamics.

[Reproduction and stored records](../analysis/response/README.md) · [Hypothesis](transport-hypothesis.md) · [Manuscript](../paper/manuscript.tex)
