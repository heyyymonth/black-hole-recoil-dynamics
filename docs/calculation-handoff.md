# Calculation handoff: recoil transport and environmental exposure

**Scope:** Step 1 supplies remnant motion; Step 2 measures residence and environmental exposure. This is the general trajectory specification. The [host and local-coefficient calculation](local-transport-calculations.md) has separate explicit configurations and numerical checks. It supplies neither a merger-conditioned launch law nor a calibrated path or residence distribution.

## Scientific contract

Both steps must refer to the same merger selection, host, launch law, physical horizon, and stopping convention. Null entries in the [input register](../formulation/inputs.json) mean unresolved inputs, not zero values, absent effects, or permission to choose silent defaults.

| Input group | Required specification | Physical consistency requirement |
|---|---|---|
| Units and frame | Physical units; translating nonrotating cluster axes; external reference orbit | Positions, velocities, recoils, potentials, densities, and times share conventions |
| Host matter | Species masses and number distribution functions; streaming, anisotropy, evolution or frozen-host validity | Density and velocity laws describe the same species, host, and epoch |
| Gravity | Potential, boundary conditions, external potential and reference acceleration | Internal gravity matches the mass inventory; remnant self-gravity is excluded |
| Launch ensemble | Joint merger position, progenitor velocity, recoil, mass, and selection weights | Rotate recoil into the host frame; preserve host/launch correlations |
| Collision model | Drift/diffusion prescription, applicable Coulomb cutoffs, validity criteria | Consistent number/mass normalization, Ito convention, and positive-semidefinite diffusion |
| Strong events | Detector, resolved objects, or justified hazard; weak/strong partition | An absent detector is not evidence of absent encounters; avoid double-counting |
| Termination | Physical horizon, escape and validity rules, simultaneous-event ordering | Separate escape, unresolved motion, numerical failure, and censoring |
| Region | Level-set function, boundary evolution, boundary-launch rule | Projected core radius does not directly define a spatial region |
| Exposure | Species weights and positive reference densities where needed | Use the same host densities, units, and normalization conventions |
| Numerical controls | Solver, temporal resolution, error controls, ensemble precision and seeds | Separate numerical and sampling errors from physical uncertainty |
| Evidence | Source versions, unit conversions, uncertainty and calibration domain | Distinguish observations, inferred model products, and simulations |

No row is fulfilled by naming a model family alone. Fitting luminous stars does not directly measure every remnant species or validate friction.

## Outputs from Step 1

Each trajectory must carry its realization, host and launch identifiers, remnant-mass label, original ensemble weight, physical time, position, velocity, frame, and units. Preserve enough model metadata to reconstruct coefficient evaluation, numerical provenance, and random seeds if used.

If a sample is selected on future survival, preserve the conditioning or reweighting of the full path law, not only the launch weights. Include any stochastic stopping-event randomness in ensemble averages.

Record the horizon, modeled endpoint, stopping cause, validity flags, and crossing events in distinct fields. A path may re-enter the core and later escape. Missing intervals must not be filled as rest. Strong events cannot be continued without a transition description.

A transport-density solution must declare its conditioning, mass coordinate or fixed-mass label, initial and boundary conditions, and absorbed or killed probability. The alive density integrates to survival probability. History-dependent events require an augmented state or a justified treatment rather than an assumed closed phase-space equation.

## Outputs from Step 2

| Output | Required interpretation |
|---|---|
| First entry, exit, and re-entry | Sign-changing finite-region events with explicit launch-side and moving-boundary treatment |
| Residence | Accumulated time in the region, with endpoint and denominator convention |
| Physical species exposure | Mass density multiplied by physical time along the path |
| Normalized species exposure | Time integral of a declared density ratio; not generally bounded by the horizon |
| Time-augmented occupation | Association of elapsed time, position, and velocity, with resolution if discretized |
| Ensemble summaries | Original weights, survival probability, unresolved/censored outcomes, and uncertainty labels |

A full-horizon physical conclusion requires coverage of that horizon. A stopped integral is an accumulated modeled quantity. Its stopping status must remain attached in every comparison.

## Decisions required before execution

1. Define the conditional population and physical interval.
2. Identify which host inputs are observations, inferred functions, simulations, or assumptions.
3. Specify the approximations being tested and their validity domains.
4. Declare the observables, independent reference, and accuracy criteria.
5. Define the treatment of correlated uncertainty, stochastic variation, numerical errors, and incomplete histories.

These decisions belong in a versioned calculation specification. The present register keeps them open.

## Acceptance conditions — unexecuted

- Check dimensions, mass/number normalization, force–potential consistency, and diffusion symmetry and positivity.
- Recover appropriate conservative and independently justified analytical stochastic limits, including the Ito energy balance.
- Test tangencies, inside/outside/boundary launches, moving regions, competing stops, and probability accounting.
- Compare path-accumulated exposure with survival-weighted density integrals where both representations exist.
- Establish convergence of the relevant distributions, crossings, residence, and exposure under a declared refinement procedure.
- Compare with independently resolved dynamics appropriate to the physical regime before claiming calibrated predictions.

This general trajectory specification assigns no trajectory solver, sample count or pass outcome. The separate host and coefficient protocols assign their own numerical checks and tolerances. Internal numerical agreement alone does not establish physical adequacy.

## Interface beyond the present scope

Timed positions and velocities can later be combined with an independently justified companion population and interaction law. Neither is supplied by the two-step outputs. This handoff contains no encounter-rate, capture, exchange, binary-evolution, or merger-probability calculation.

[Problem statement](problem.md) · [Step 1](e001-ballistic-return.md) · [Step 2](e002-core-exposure.md)
