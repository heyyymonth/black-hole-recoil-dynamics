# Independent test of the residence response

The [precision comparison](precision-pilot.md) supports symmetric-root coupling as a substantially more precise estimator of the same finite-step response. The scientific question remains whether that response predicts the change in finite-time residence when diffusion transitions are replaced by finite encounters. This calculation supplies independent predictions and validation observations, while retaining numerical-resolution and continuation uncertainties.

## Fixed physical comparison

The conditional Omega Centauri host, 60-solar-mass remnant, BH encounter law, all-species gravity, three launch states, residence region and horizon are identical to the [earlier response test](response-test.md). The numerical domain remains radius and speed below four in their declared units. Existing kernel tables do not support silently extending that domain. Excluded-channel killing is part of the modeled process, distinct from numerical-domain termination and from capture.

For the complete numerical transition mixture `P_a=(1-a)P_D+aP_J`, including rewards and stopping, define `Y=C_T/T`, `mu_h(a)=E_a[Y]`, and `S_h=mu_h'(0)`. Prediction uses the exact single-insertion identity with the existing conditioned quartet. The symmetric PSD covariance root is used in every diffusion update. It preserves the marginal transition law and changes the coupling used to estimate differences.

The validation statistic is

$$
R_h(a)=\left[\mu_h(a)-\mu_h(0)\right]-aS_h.
$$

Independent prediction and validation batches give estimated sampling variance

$$
\widehat{\operatorname{Var}}(\widehat R_h(a))
=\frac{s_{\Delta Y}^2}{n_v}
+a^2\frac{s_S^2}{n_p}.
$$

Mixture and diffusion rewards are paired within each validation row; their covariance is already included in `s_DeltaY`. Different validation batches and every prediction batch have disjoint random streams. Reusing a prediction at several mixture values makes their residuals dependent; the declared family adjustment does not require residual independence.

## Sample allocation and registration

| Setting | Prediction quartets per launch | Validation pairs per launch and mixture value | Mixture values |
|---|---:|---:|---|
| Coarse table, step 0.02 | 8,192 | 8,192 | 1 |
| Fine table, step 0.02 | 8,192 | 8,192 | 1 |
| Fine table, step 0.01 | 65,536 | 8,192 | 0.05, 0.1, 1 |

Steps are in the common dynamical-time unit; the horizon is ten. There are nine prediction ensembles and fifteen validation ensembles: 245,760 independent response quartets and 122,880 independent validation pairs, containing 1,228,800 correlated branch trajectories. All nine predictions must be complete and hash-registered before any validation begins. No previous samples enter these estimates.

Counts are fixed before execution. The prior symmetric/uniform variances imply optimistic primary half-widths near 0.0164, 0.0170 and 0.0116 under the present confidence allocation. They are much larger than the unchanged 0.0005 precision target, and unseen rare scores can invalidate those projections. The allocation is a resource limit, not a promise of decisive inference. No count, target, cutoff or acceptance threshold is changed after outcomes.

Two independent workers use private cache directories. Existing immutable sampler arrays are hard-linked into each worker's cache, with the same physical-input namespace and checksums. Each directory has one writer. Execution order can overlap across workers, while the prediction-before-validation gate remains mandatory. Memory and disk caching change storage, not numerical sampler values or random-number laws.

## Sampling and resolution checks

The normal-approximation family contains 51 statements: nine responses, fifteen measured shifts, fifteen prediction residuals, six response-refinement contrasts and six endpoint-refinement contrasts. A conservative family size of 64 receives failure budget 0.025. Small-mixture residuals divided by `a` reuse the original interval event; they are finite-`a` checks, not an exact derivative-limit test.

The endpoint error tolerance remains 0.001. An interval wholly inside its closed band supports agreement for the declared stopped numerical model; an interval wholly outside resolves disagreement at that tolerance; crossing a tolerance boundary is inconclusive. Sampling precision is assessed separately. Coarse/fine table and 0.02/0.01 time-step contrasts retain the 0.0005 sensitivity budget. A finite contrast is not a rigorous continuum truncation bound.

All scores are retained. Extrema, absolute-score concentration and centered-square variance concentration remain visible. Normal intervals may be unreliable for rare influential contributions even at a large sample count.

## Population uncertainty from numerical termination

Let `B_i` be the retained row-level upper allowance for missing normalized residence. It sums the two validation-branch allowances, or weights the four response allowances and multiplies by the scheduled-step count `N`. Thus `0 <= B_i <= M`, where `M=2` for a validation pair and `M=2N` for a response quartet. Branches within a row are correlated; the independent statistical units are rows.

A finite-sample mean allowance does not control its population expectation. If `k` of `n` independent rows have positive allowance, the one-sided Clopper–Pearson occurrence limit gives

$$
u=\operatorname{Beta}^{-1}(1-\delta;k+1,n-k),\qquad
\mathbb E[B]\leq M\Pr(B>0)\leq Mu.
$$

For `k=n`, set `u=1`. For no observed events, use `u=1-delta^(1/n)`, which is positive. Each of the 24 continuation bounds receives `delta=0.025/24`. The binomial coverage is exact under fixed-count iid assumptions; numerical quantiles use float64 evaluation with an outward guard. The response/validation interval family remains approximate, so the combined nominal error allocation is not an exact joint normal-interval guarantee. [Clopper and Pearson, 1934](https://doi.org/10.1093/biomet/26.4.404)

Even with no exits, this allocation gives primary prediction upper allowance about 0.210, control prediction about 0.838, and validation about 0.00168. These exceed the unchanged 0.0002 continuation budget. The worst-case range bound is conservative, but it makes the unresolved population uncertainty explicit. A smaller observed mean cannot override it.

For the response residual, a conservative continuation allowance is the validation upper bound plus `a` times the prediction upper bound. This bounds possible changes under bounded-reward continuations agreeing with the numerical process up to termination. It does not specify or validate dynamics outside the table domain. Physical calibration and continuous-time convergence remain separate requirements.

[Fixed protocol](../analysis/retest/protocol.json) · [Implementation and verification](../analysis/retest/README.md)

## Recorded findings

The fixed protocol and implementation were published in commit `487ac443` before production trajectories. All nine predictions were then saved and hash-registered before the fifteen independently sampled validations. The completed design contains **368,640 independent groups and 1,228,800 correlated branch trajectories**. No earlier samples are pooled with this test.

| Primary launch | Response estimate | Standard error | Response interval half-width | Endpoint residual interval | Accuracy |
|---|---:|---:|---:|---|---|
| Radial outward | -0.011685 | 0.003543 | 0.012565 | [-0.010392, 0.016622] | Inconclusive |
| Tangential | -0.002609 | 0.005037 | 0.017861 | [-0.019211, 0.018409] | Inconclusive |
| Oblique | -0.011350 | 0.003984 | 0.014128 | [-0.010782, 0.019199] | Inconclusive |

The response half-width target remains 0.0005 and the endpoint error tolerance remains 0.001. Every primary residual interval includes both acceptable and unacceptable errors. None establishes agreement or failure of the prediction. All six small-mixture scaled-residual intervals contain zero; that does not establish the derivative limit.

At equal primary sample count, the observed standard errors are smaller than those of the earlier eigenbasis test by descriptive factors 13.1, 33.8 and 41.0. These ratios have no additional confidence interval, and earlier validation is not used to certify current accuracy. The independently measured endpoint shifts are −0.008569, −0.003010 and −0.007141. Radial and oblique shift intervals exclude zero under the approximate family adjustment; the tangential interval includes zero. An encounter effect and successful prediction of that effect are different claims.

The largest one percent of centered-square score contributions accounts for 89.9%, 90.6% and 91.2% of the three primary sample variances. The largest individual shares are 2.84%, 2.49% and 8.76%. All scores are retained. These observed concentration diagnostics do not establish normal coverage or account for unseen tails.

**All twelve refinement classifications remain inconclusive. None of the 24 population continuation bounds meets the 0.0002 budget.**

| Primary launch | Positive-allowance prediction rows | Empirical response allowance | Population response upper bound | Empirical endpoint allowance | Population endpoint upper bound |
|---|---:|---:|---:|---:|---:|
| Radial outward | 27 | 0.125621 | 1.438373 | 0.00007609 | 0.002728 |
| Tangential | 21 | 0.035061 | 1.198915 | 0.00000000 | 0.001676 |
| Oblique | 25 | 0.121685 | 1.359293 | 0.00011472 | 0.003174 |

Response bounds can exceed one because they concern a step-weighted derivative score, not a probability. A large bound is not a measured physical continuation error. All continuation-expanded endpoint residual intervals remain inconclusive. Their sampling component is still a normal approximation; the exact binomial construction does not turn the combined envelope into an exact joint confidence interval.

This test makes the statistical comparison more precise and identifies the remaining uncertainty separately. Further work needs a more precise response estimate, resolved transition refinement, and a sharper justified continuation bound or an extended verified dynamical domain. Additional histories alone cannot supply missing physical dynamics. The present calculation does not establish continuum response accuracy, physical-cluster calibration, or a merger rate.

[Recorded summary](../results/retest/summary.json) · [Comparison figure](../results/retest/comparison.png) · [Prediction register](../results/retest/prediction-register.json)
