# Host structure and data requirements

Transport and exposure require a consistent environment. Potential, mass density, velocity distributions, stellar species, and their evolution must refer to compatible hosts and epochs.

## What data constrain

Surface brightness, stellar counts, line-of-sight velocities, proper motions, and stellar mass information constrain different projections of cluster structure. Obtaining a three-dimensional mass distribution or potential requires a model, distance and selection treatment, and uncertainty propagation. An inferred potential is not a directly observed force field.

The scattering population can include faint stars and remnants poorly constrained by luminous tracers. Total mass alone does not determine species-resolved phase space. A present-day host fit also need not describe the environment at a historical merger.

The [observed-cluster assessment](cluster-data-assessment.md) audits original public products for six candidates and identifies what can be put into common units. Its [host schema](cluster-host-schema.md) retains distinct joint parameter distributions, tracer selections and source conventions. This data comparison does not yet reconstruct the local distributions needed for drift and diffusion.

## Host representations and tradeoffs

| Representation | Information available | Further requirements |
|---|---|---|
| Distribution-function family fitted to observations | Joint density/velocity description with declared mass components and anisotropy | Fit uncertainty, species content, model adequacy, and evolution over the recoil interval |
| Empirical or inferred profiles | Mass, potential, or velocity moments on a measured domain | Compatible phase space, geometry, boundaries, and treatment outside that domain |
| Simulation snapshots | Species inventories, phase-space samples, model-dependent structure | Time conventions, sampling effects, identities, resolution, and consistency across epochs |
| Resolved dynamical histories | Correlated motion and encounters | Numerical accuracy, appropriate masses and relaxation, and physical applicability |

[Gieles & Zocchi (2015)](https://arxiv.org/abs/1508.02120) supplies a multimass lowered-isothermal family with truncation and anisotropy. It is a possible mathematical description, not a selected fit. [Tep, Fouvry & Pichon (2022)](https://arxiv.org/html/2201.03985) addresses anisotropy and local diffusion. [Webb et al. (2018)](https://arxiv.org/abs/1711.09100) motivates validating kicked-remnant drift and scattering instead of treating a drag law as exact.

Independently resampled snapshot positions are not necessarily a continuous object trajectory. Data products must document the meaning of positions, velocities, and time before they support return or residence.

## Connection to the two steps

Step 1 consumes host distributions, gravity, a joint launch law, collision coefficients, and stopping rules. Step 2 consumes those histories and the same host densities at the same times, together with a region and exposure convention.

The [input register](../formulation/inputs.json) chooses no representation or cluster. The [handoff](calculation-handoff.md) states the evidence required. Companion populations and interaction laws are separate information for later questions.
