# Scientific relevance review

Initial baseline review: 18 October 2026. E003 relevance and protocol review: 19 October 2026, before E003 execution.

## E001: relevant as a numerical benchmark

The first-return calculation establishes a checked orbit solution and exposes the distinction between binding and return time. It is a derivation in a standard potential, not a claim of new astrophysical physics. Independent orbit integration and analytical limits are appropriate validation.

## E002: relevant as a controlled approximation test

Central residence and density exposure answer a specific question that first-return time alone cannot answer: what fraction of a finite window is actually spent in a concentrated environment? The continuous trajectory and the two deliberately simplified comparisons make that question measurable without inventing a physical merger efficiency.

The comparison must preserve both possible directions of error. Waiting until the first return removes genuine early exposure. Assuming permanent central residence thereafter eventually overestimates exposure, because the undamped orbit continues to leave the center.

## Relationship to existing work

[Islam (2026), sections V.1.3-V.1.4](https://arxiv.org/html/2604.04546v2), already investigates remnant displacement, dynamical-friction recentering, and the suppression of repeated mergers. Its phenomenological probability depends on a frictional return time. Our first ballistic center crossing is a different quantity. These experiments neither reproduce nor falsify that probability prescription.

[Torniamenti et al. (2024)](https://arxiv.org/abs/2401.14837v2) already includes host evolution and uncertain pairing in hierarchical-merger predictions. A static-potential benchmark does not establish its applicability to those evolving systems.

This focused check supports publication of the experiments as reproducible benchmarks. It does not establish publication-level novelty or an observational discriminator.

## E003: relevant as a test of the density-to-encounter approximation

The [pre-execution protocol](e003-protocol.md) fixes an explicit target population, close-passage definition, common time window, and 10% approximation screen. These are the missing ingredients needed to test whether the E002 density proxy represents passage opportunities. The [result](e003-encounter-accessibility.md) identifies a velocity-dependent discrepancy and its geometric/focusing contributions. This is a meaningful physical extension, whereas a denser deterministic kick grid would only refine the same earlier calculation.

The apocenter admission criterion is motivated by Islam's stated choice of five half-mass radii. It admits E003's eta=0.9 orbit but excludes E002's eta=0.99 near-escape example. This does not invalidate the isolated-orbit benchmark; it limits its environmental interpretation. E003 retains an untruncated equilibrium Plummer bath and does not implement a self-consistent tidal boundary.

The rate ingredients are standard. [Leigh et al.](https://arxiv.org/abs/1711.10494) already studies encounters and binary hardening, and [the galpy distribution reference](https://docs.galpy.org/en/v1.12.0/reference/dfplummer.html) supplies the exact Plummer distribution used for the analytical velocity calculation. We have checked relevance and correctness of this limited extension, not established a novel assembly model.

The completed screen rejects the density proxy for two selected larger kicks. Absolute passage counts remain tiny for the imposed bath. Neither the discrepancy nor passing numerical checks establishes an observable merger-rate effect. Physical calibration of the bath and trajectories, and conditional encounter outcomes, remain outside this experiment.

## E004: relevant as an external assumption check

Reviewed 20 October 2026. The [executed comparison](e004-cluster-validation.md) tests the assumed target population and host potential using a separately evolved cluster model. It identifies sparse BH–BH targets and a return-time discrepancy much larger than shell-compression error. The external simulation adds physical context, but a single frozen snapshot cannot validate dissipative motion or infer a merger rate. The [literature review](cluster-dynamics-context.md) distinguishes established mass segregation and encounter frameworks from two unexecuted experiments needed to calibrate this model.

## E005: relevant as a target-availability and data-sufficiency audit

Reviewed before execution on 21 October 2026 in [the protocol](e005-protocol.md). E004's two-target state motivates checking the scope of that finding before fitting a continuous companion model. The [executed audit](e005-target-availability.md) extends the comparison to the available late-time summary release, with one age-selected state per model and explicit missing coverage. It finds small instantaneous BH–BH populations throughout the inspected saved states, but does not invalidate a continuous ensemble description or identify distinct binaries across time.

The next physical comparison must measure changing target availability and joint properties conditional on merger events. Equal-weight old-cluster states are not that population. This is a relevant empirical constraint on our modeling choices; it is not evidence of a novel formation channel, a calibrated encounter rate, or validated live recoil motion.
