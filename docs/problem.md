# Recoil transport and environmental exposure of black-hole merger remnants

**Researcher: Venkata Hemanth Athota**

## Problem statement

A black-hole merger can leave a remnant moving through a stellar cluster under the combined effects of its recoil, the progenitor's existing motion, the host gravitational field, and gravitational encounters. Binding to the host does not by itself specify the remnant's location, speed, return to a dense region, or residence within that region.

The physical question is:

> How do merger-conditioned initial states and a consistent stellar environment determine the probability distribution of recoil trajectories and the environmental exposure accumulated over a finite physical interval?

The study has two components. [Step 1](e001-ballistic-return.md) defines the phase-space evolution. [Step 2](e002-core-exposure.md) defines the functionals of that evolution that measure entry, residence, and species-resolved exposure. An [observational analysis of six clusters](cluster-data-assessment.md) supplies numerical comparisons of candidate host inputs, and a [fixed-host reconstruction with local collision calculations](local-transport-calculations.md) supplies the first numerical connection between those inputs and the transport operator. A [static boundary calculation](outer-boundary.md) separates a conditional tidal opening from the accuracy domain of neglecting the Galactic field. A [Galactic-orbit calculation](galactic-orbit.md) now supplies phase-dependent external forcing and tests the fixed host's force and response-timescale domains. Remnant trajectories and residence outcomes remain uncalculated.

## Observational findings and their role

The completed analysis compares original projected kinematic profiles and published host posteriors using common variables and units. It retains source selection, shared distance uncertainty, radial coverage and joint parameter dependence. Its findings support one host interface with a distinct conditional distribution for each cluster; they do not establish a universal density profile or calibrated local scattering background.

The [paper](../output/pdf/central-residence.pdf) presents the six-cluster coverage, structural scales, posterior summaries and model comparisons as findings of this study. Parameter dependence and differences between structural estimates identify assumptions that a transport calculation must preserve or test. Neither observation directly predicts how long a kicked remnant returns or resides.

A testable physical hypothesis is that host structure and mass composition affect consistently normalized residence and species exposure under fixed dimensionless launch, remnant-mass, external-field, collision and event conventions. Testing it requires reconstructed, observationally checked host realizations and controlled transport calculations. The current results motivate the comparison without assuming its outcome.

## Quantities to determine

- The joint distribution of position and velocity, conditional on remnant mass, host, and merger selection.
- Survival within a declared modeling domain and the probabilities of distinct stopping causes.
- Entry or re-entry into a finite region, with explicit treatment of launches outside or on its boundary.
- Residence and density exposure over a common time interval, including incomplete trajectories.
- The time-resolved occupation measure needed to evaluate environmental kernels without discarding velocity or temporal ordering.

A scalar return time or a surviving-sample average cannot stand in for all of these quantities.

## Mathematical formulation

Let the state be position, velocity, and a fixed remnant-mass label between discrete transitions. The host supplies species number distributions, a compatible potential and external reference orbit, and an explicit scattering closure. The initial state follows a joint merger-conditioned law. No factorization among launch variables, host properties, and merger selection is assumed.

The transport problem requires an initial distribution, differential operator, spatial/velocity domain, boundary conditions, stopping rules, and validity interval. The exposure problem requires a physical time horizon, an operational region, species weight functions, and an accounting convention for stopped paths. These constitute a fully specified calculation only when the [input register](../formulation/inputs.json) is populated and its physical consistency justified.

## Scope of the physical assumptions

The equations use nonrelativistic cluster motion in translating, nonrotating coordinates. A relativistic merger model supplies the recoil and remnant mass. The remnant is a test object in a supplied background until a declared transition or validity cutoff. A local diffusion closure assumes weak, uncorrelated encounters on appropriate scales. The host's evolution and its backreaction are not solved implicitly by Poisson's equation.

The host need not be spherical or isotropic. Its time dependence, rotation, anisotropy, mass spectrum, and species concentration must be represented consistently where retained. Leaving a dependence symbolic does not establish observational support or numerical accuracy.

## Research objective and evidence standard

The objective is to determine which launch and environmental information must be retained to predict finite-time residence and exposure at a declared accuracy. Useful future comparisons would hold the host, launch ensemble, time horizon, and observable definitions fixed while testing a specified simplification. The comparison, accuracy requirement, and uncertainty treatment must be declared before interpreting its outcomes.

A justified conclusion could identify a regime where a reduction preserves the observables, or a reproducible discrepancy showing why additional information is required. Neither outcome is assumed. Analytical identities establish internal consistency; numerical convergence establishes implementation accuracy; comparison with independently resolved dynamics is needed for physical calibration.

## Connection to subsequent physics

The trajectory record retains times and velocities because a future companion-access model would need to match the remnant's motion to an independently justified target population. Density exposure alone is not an encounter count. Companion acquisition, binary evolution, merger probabilities, and inference from gravitational-wave observations require further models and are outside steps 1 and 2.

[Step 1](e001-ballistic-return.md) · [Step 2](e002-core-exposure.md) · [Calculation handoff](calculation-handoff.md) · [Scientific relevance](relevance.md)
