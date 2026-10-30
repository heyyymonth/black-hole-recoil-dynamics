# Recoil transport and environmental exposure of black-hole merger remnants

Independent physics research · Stellar dynamics and gravitational-wave source environments

**Contributed by Hemanth.**

[![Research verification](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml/badge.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[Problem statement](docs/problem.md) · [Step 1](docs/step1-transport.md) · [Step 2](docs/step2-exposure.md) · [Calculation handoff](docs/calculation-handoff.md)

## Abstract

Black-hole mergers in stellar clusters leave remnants whose recoil can carry them repeatedly through regions of different stellar density and composition. Their subsequent availability for dynamical interactions depends on the time spent in those environments, as well as on whether they remain bound. We study how gravitational encounters determine finite-time residence and species-resolved exposure. Published multimass distributions for six Galactic clusters provide the adopted potentials and scattering backgrounds. Reconstructions, local anisotropic collision calculations and Galactic environmental controls connect these backgrounds to a merger-conditioned transport formulation. The calculations reveal distinct radial profiles of the black-hole contribution to velocity diffusion. Finite-angle encounter calculations specify a common collision law and separate event frequency from individual deflection strength. We formulate the hypothesis that the accuracy of continuous transport for residence depends on the encounter increments sampled along the orbit and on how those increments change the expected remaining occupation. A matched comparison retains either the full local encounter transitions or their first two velocity-change moments. A backward occupation equation relates their difference in mean residence to a nonlocal generator residual, defining a route from local scattering physics to an observable accuracy criterion. Finite-population response and time-dependent host evolution enter as separate controlled comparisons. Collisionless controls connect these backgrounds to finite-time residence and encounter budgets along equal-energy launches. The matched collisional residence predictions and their response estimate remain to be evaluated.

## Research question

When does a continuous drift-and-diffusion description preserve a kicked remnant's finite-time residence and exposure, and which encounter or host-response effects must be retained to reach a declared accuracy?

The [transport hypothesis](docs/transport-hypothesis.md) connects the completed local calculations to a matched encounter comparison and a backward equation for expected remaining occupation.

The active scope consists of two linked steps:

| Step | Physical task | Output required |
|---|---|---|
| [1. Recoil transport](docs/step1-transport.md) | Evolve a joint merger-conditioned launch under consistent host gravity, external acceleration, collisional drift, and velocity diffusion. | Timed position–velocity histories or a transport density, with mass labels, uncertainty, survival probability, and stopping causes. |
| [2. Residence and exposure](docs/step2-exposure.md) | Define finite-region entry and re-entry, residence, and species-resolved density exposure on those histories. | Event distributions, accumulated exposure, and a time-resolved occupation measure with explicit censoring. |

The [six-cluster calculation comparison](docs/six-cluster-comparison.md) repeats the reference host and local-transport experiments across the shortlist, with a weighted joint-row pilot, explicit [data-suitability decisions](docs/six-cluster-data-quality.md), and [environmental controls](docs/six-cluster-environment.md). An [executed comparison notebook](analysis/cluster_comparison/cluster_comparison.ipynb) presents the recorded results.

The [finite-angle encounter calculation](docs/finite-angle-encounters.md) evaluates 162 local settings across the six hosts and checks 600,000 sampled elastic increments against their matched moments. It also supplies a tested numerical first-response estimator; matched collisional residence predictions remain unevaluated.

The [collisionless residence controls](docs/collisionless-residence.md) integrate two equal-energy launch directions in each host. They verify residence through independent radial and Cartesian calculations and accumulate BH encounter budgets along the paths. Collisions are disabled in these reference trajectories.

The detailed single-host benchmark constrains one candidate environment. A completed NGC 3201 reconstruction now supplies consistent species distributions, and twenty evaluated local states test its anisotropic drift and diffusion. Four Galactic reference orbits now supply time-dependent tides and a [host-support audit](docs/galactic-orbit.md). Collisional remnant evolution and residence distributions remain unevaluated; the sparse black-hole component and the host in a Galactic tide require physical validation. A [boundary calculation](docs/outer-boundary.md) gives a conditional tidal-opening scale and a separate accuracy limit for isolated controls. No central launch, radial orbit, permanent rest after return, or universal cluster profile is imposed.

## Physical foundation

The host's species distributions determine both gravitational mass density and the local scattering environment. Initial remnant motion includes the progenitor's center-of-mass velocity and the recoil vector in a common coordinate system. Mean drift and random deflections retain the remnant's continuing motion through and beyond a declared core region.

The local weak-scattering closure has explicit applicability assumptions. It does not resolve strong encounters, evolve the host by itself, or guarantee calibrated recoil trajectories. [Mathematical reading order](docs/recoil-transport-theory.md) · [Host and data requirements](docs/cluster-dynamics-context.md) · [Primary literature](docs/references.md).

Timed position and velocity records can support later companion-access calculations. Local elastic BH encounter rates are now evaluated; binary acquisition, capture outcomes and repeated-merger predictions require additional interaction physics.

## Observed-cluster data

The [six-cluster assessment](docs/cluster-data-assessment.md) examines original public observations and published candidate host posteriors for 47 Tuc, NGC 1851, NGC 3201, Omega Centauri, M4 and NGC 6397. The [data package](data/observed-clusters/README.md) includes source hashes, common-unit tables, quality flags, three comparison figures and an [executed notebook](analysis/cluster_data/cluster_data_audit.ipynb).

A [common host interface](docs/cluster-host-schema.md) preserves a separate weighted joint parameter distribution for each cluster. Unit harmonization does not establish a shared stellar population, a calibrated potential or a local scattering background. The [fixed-host calculation](docs/local-transport-calculations.md) reconstructs one NGC 3201 joint realization; exact historical fit reproduction and physical calibration are not claimed. The manuscript presents the data findings and the completed host and collision calculations in ten tables and ten vector figures. It records 107 original HST radial bins, 621,878 Gaia field-catalog rows, 20 flagged compilation rows and strong mass–radius dependence in several published host posteriors. Field-catalog rows include nonmembers, and the display uncertainties condition on catalog point scales.

## Calculation handoff

The general [trajectory input register](formulation/inputs.json) leaves the launch law, globally adopted collision prescription, stopping rules and trajectory controls unassigned. The finite-angle local protocol separately declares its conditional cutoffs. The completed host and local-coefficient calculations have separate versioned protocols, linked from that register. Missing inputs are not zero-valued defaults. The [handoff specification](docs/calculation-handoff.md) defines their units, dependencies, required output fields, and acceptance checks before a calculation can be interpreted.

From the repository root, create a Python 3.11 environment and install the active checks:

```sh
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pytest -q
```

The default test suite covers all active scientific stages and regression checks. Reconstructing a host uses its separate pinned [scientific environment](analysis/host_reconstruction/README.md). Integrity checks perform no physical calculations:

```sh
python scripts/verify_foundation.py
python paper/verify.py
python analysis/cluster_data/verify.py
python analysis/host_reconstruction/verify.py
python analysis/boundary/verify.py
python analysis/orbit/verify.py
python analysis/cluster_comparison/verify.py
python analysis/encounters/verify.py
python analysis/residence/verify.py
```

The [research paper](output/pdf/central-residence.pdf), [LaTeX source](paper/manuscript.tex), and [portable source archive](output/pdf/manuscript-source.zip) combine the two-step formulation with the observed-cluster findings. It is an unreviewed manuscript reporting data analysis, six conditional reference hosts, joint-row sensitivity, local collision kernels and environmental controls. It also reports conditional collisionless trajectories and their occupation, without claiming calibrated collisional residence predictions. [Paper status and build instructions](paper/README.md) · [Calculation audit](docs/repository-audit.md).

## Repository contents

| Path | Purpose |
|---|---|
| `docs/problem.md` | Research question, scope, and mathematical objectives |
| `docs/step1-transport.md` | Step 1: recoil transport |
| `docs/step2-exposure.md` | Step 2: residence and environmental exposure |
| `docs/calculation-handoff.md` | Required inputs, outputs, and validation conditions |
| `formulation/` | General trajectory inputs and links to defined calculation stages |
| `data/observed-clusters/` | Public-source manifest, candidate host records and audited comparison products |
| `analysis/cluster_data/` | Source retrieval, reproducible data audit, figures and executed notebook |
| `analysis/host_reconstruction/`, `analysis/transport/` | Reproducible host and local collision calculations with independent checks |
| `results/host-reconstruction/`, `results/local-transport-kernels/` | Completed calculation records and shared provenance manifest |
| `analysis/boundary/`, `results/outer-boundary/` | Static tidal and validity diagnostics with primary-source provenance |
| `analysis/orbit/`, `data/orbit/`, `results/orbit/` | Galactic reference orbits, original phase-point inputs and time-dependent host diagnostics |
| `analysis/cluster_comparison/`, `results/cluster-comparison/` | Six conditional references, joint-row sensitivity, data-quality decisions and comparisons |
| `data/cluster-comparison-environment/`, `results/cluster-comparison-environment/` | Per-cluster source phase points, Galactic orbits and conditional boundary diagnostics |
| `analysis/residence/`, `results/collisionless-residence/` | Collisionless reference trajectories, independent occupation checks and path-integrated encounter budgets |
| `scripts/` | Structural documentation and manuscript integrity checks |
| `paper/` | Formulation and observational findings, evidence-bound tables and figures, bibliography, and typesetting build |
| `tests/` | Regression checks for input validity, configuration handling and evidence consistency |
| [Historical records](docs/computational-records.md) | Immutable source and results for earlier studies |

## Citation and reuse

The [BibTeX citation](CITATION.bib) and [CITATION.cff](CITATION.cff) identify the [versioned research formulation](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/formulation-v1), separately from an executed study or a journal publication. No DOI has been assigned.

Original code and documentation are available under the [MIT license](LICENSE). Referenced publications, dependencies, and external data retain their own rights. [Contribution guidance](CONTRIBUTING.md) · [Author and contact](AUTHORS.md).
