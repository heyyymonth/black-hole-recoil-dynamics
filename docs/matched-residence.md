# Matched encounters and remnant residence

The question is whether replacing finite BH encounters by their drift and velocity diffusion preserves finite-time core occupation. Local moment matching does not by itself answer that question: future excursions depend on the full velocity-change distribution and on when and where the encounters occur.

This calculation uses the conditional Omega Centauri reference. It was selected before the collisional outcomes because its collisionless launch paths had the largest integrated allowed BH encounter rate among the six reference hosts. The experiment connects the [local encounter law](finite-angle-encounters.md) and [collisionless controls](collisionless-residence.md) to actual stochastic residence. It does not fit new cluster data or supply a calibrated merger population.

## Fixed physical inputs

- Frozen gravity from all physical species in the same host realization.
- A 60-solar-mass remnant and BH collisions only. Lighter-species collisions are disabled identically in both arms; their smooth gravity remains.
- Initial radius `0.5 r_rho`, with radial-outward or tangential velocity at equal energy. The speed corresponds to a collisionless radial apocenter of `2 r_rho`.
- Residence inside `r < r_rho`, with every visit counted over `T = 10 t_dyn`. No remnant is brought to rest on return.
- The existing allowed elastic encounter domain, `b_max = 0.1 L`, and the same excluded close-passage flux. An excluded event ends modeled accumulation without assigning a capture outcome.
- Numerical limits `r < 4 r_rho`, `speed < 4 v_dyn`, inside the inherited isolated-force domain. Crossing these limits is censoring, not escape.

The scales are `r_rho = 0.0247824 pc`, `v_dyn = 11.6746 km/s`, `t_dyn = 0.00207562 Myr`. The diagnostic launch speed is `19.8835 km/s`. These are properties of the adopted conditional host, not measured merger parameters.

## A shared numerical collision law

The coarse grid contains 550 radius/speed/pitch vertices; bisecting every interval gives 3591 refined vertices. Positive interpolation mixes whole vertex rate measures. An event chooses vertex `i` with probability proportional to `w_i Gamma_i`. The Gaussian arm uses the first moment and **raw second moment per unit time** of exactly that same mixture. Both arms share velocity quadrature, local symmetries and coordinate transformations.

Every vertex passes its local refinement requirements. At 72 local checks the maximum relative diffusion-moment discrepancy falls from 9.24% to 1.97% on grid refinement, and the maximum scaled drift discrepancy falls from 4.57% to 0.963%. These finite checks are not a uniform error bound. Independent sampling at six states per grid, with 40,000 increments each, passes the six-standard-error consistency diagnostic for means and raw second moments. This verifies the sampled numerical law; it does not demonstrate continuum or residence convergence.

At each step the local collision law is frozen. An excluded-event clock may shorten the duration. Two gravity half-flows surround either one Gaussian increment or the sum of a Poisson number of finite increments. Multiple finite events are allowed and their sum is applied at the midpoint. Their individual event times are not resolved. The first two moments agree for this frozen velocity update, while the changing trajectories require time-step refinement.

Core occupation is exact along the solver's drift chords. Density exposure uses midpoint quadrature along those same chords. Spatial exits are located on the chords; speed exits are checked at resolved endpoints. All stopped paths retain the full declared horizon denominator.

## Numerical controls and inference

The protocol fixes two grids, two steps (`0.04` and `0.02 t_dyn`), two launch directions and two collision arms, with 4096 independent paths per setting. All sixteen settings are retained. The statistical intervals are descriptive and marginal; paths are not paired across settings.

- Independent-ensemble normal intervals describe mean residence and exposure contrasts.
- Union-DKW bounds describe the population CDF-supremum discrepancy at each numerical setting. Their half-width is about 0.0463 for two samples of 4096.
- Wilson intervals describe stop probabilities, including a nonzero upper bound when no stop was observed.
- A numerical stop at time `tau` leaves at most `(T - tau)/T` unmodeled normalized residence. Both its empirical mean and a conservative probability-based upper bound are retained.

The predeclared targets are a mean residence difference of 0.01, a numerical mean sensitivity budget of 0.005, a CDF distance of 0.05 and a missing-residence budget of 0.002. Finite grid/time changes are sensitivity tests with sampling uncertainty, not rigorous truncation-error bounds.

Independent null controls reproduce the collisionless residence and exposure to better than `1.30e-5` in normalized units and improve with step refinement. Gravity agrees with an independent all-species density integral to relative discrepancy `5.46e-4`, below its `1e-3` criterion. Analytical tests cover anisotropic moments, multiple events, exponential killing, competing stops and occupation accounting.

## Measured outcomes

At the refined grid and `dt = 0.02 t_dyn`:

| Launch | Jump mean residence | Diffusion mean residence | Jump minus diffusion, 95% interval | CDF distance, 95% union-DKW interval |
|---|---:|---:|---|---|
| Radial outward | 0.30400 | 0.31158 | -0.00757 [-0.01126, -0.00389] | 0.12646 [0.08021, 0.17272] |
| Tangential | 0.32426 | 0.32924 | -0.00499 [-0.00952, -0.00045] | 0.10547 [0.05921, 0.15173] |

Both arms exceed the collisionless residence controls, 0.27194 and 0.28866. This is the combined effect of the specified collisions; drift and random deflections are not separately varied.

Mean BH exposures are 0.40837 versus 0.41317 for radial launches and 0.41461 versus 0.41661 for tangential launches. Their signed contrast intervals are [-0.00735, -0.00225] and [-0.00483, 0.00083]. The tangential exposure contrast is unresolved at this sampling precision.

Across all eight arm comparisons, empirical residence CDF distances span 0.09692–0.12646. Each marginal DKW interval lies above the registered 0.05 threshold at its numerical resolution. No simultaneous claim is made. The means alone would miss this distributional difference.

A descriptive diagnostic finds 88 radial and 130 tangential zero-event, full-horizon jump paths in the refined-step settings (2.15% and 3.17%). They preserve the numerical collisionless residence. This discrete component is smaller than the observed CDF separation and does not explain that separation by itself.

**The continuum accuracy classification remains unresolved.** Some absolute mean-refinement intervals extend beyond the 0.005 numerical budget, reaching 0.00912. In the refined-step tangential comparison, the one diffusion spatial stop gives an empirical missing-residence bound of `5.39e-5`; the jump arm has no numerical stop. The sum of their marginal Wilson upper endpoints is 0.00232, above the 0.002 missing-residence budget. That sum is not a joint 95% bound. All sixteen ensembles retain ten spatial stops and one excluded-channel stop, with no speed-domain stops.

The supported finding is that the two **matched numerical laws** give different residence distributions in this conditional setting. It is not yet a converged continuum result, an observed property of Omega Centauri, or a validation of the response predictor.

## Connection to the hypothesis

The theoretical hypothesis concerns an **observable-specific response**, not a universal threshold in cluster BH count. This pilot evaluates the full finite-encounter endpoint and its matched diffusion endpoint. The diffusion-based first-response predictor has not been evaluated, so these results cannot confirm its accuracy.

A difference in the residence distribution can coexist with a smaller mean difference. The occupation residual in the paper predicts a mean; validating it would not certify the entire distribution. Distributional comparisons and stopping probabilities therefore remain separate observables. Grid/time precision must be resolved before assigning continuum accuracy, and independently resolved dynamics are needed for physical calibration.

[Protocol and reproduction](../analysis/matched/README.md) · [Retained summary](../results/matched-residence/summary.json) · [Numerical controls](../results/matched-residence/controls.json) · [Sampler validation](../results/matched-residence/sampling-validation.json) · [Research hypothesis](transport-hypothesis.md)
