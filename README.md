# Recoil transport and environmental exposure of black-hole merger remnants

Independent physics research · Stellar dynamics and gravitational-wave source environments

**Contributed by Hemanth.**

[![Research verification](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml/badge.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[Problem statement](docs/problem.md) · [Step 1](docs/e001-ballistic-return.md) · [Step 2](docs/e002-core-exposure.md) · [Calculation handoff](docs/calculation-handoff.md)

## Abstract

Anisotropic gravitational-wave emission can impart a recoil velocity to the black hole formed by a merger. In a stellar cluster, the resulting motion depends on the merger location, the progenitor's center-of-mass velocity, the cluster's mass and velocity distributions, and gravitational scattering. Remaining bound does not determine when a remnant enters a dense region or how long it occupies that region. We formulate a two-stage description of this problem. A merger-conditioned phase-space transport model combines smooth gravity, external differential acceleration, mean collisional drift, and stochastic velocity diffusion. Finite-region entry, residence, and species-resolved environmental exposure are then defined on the resulting trajectories, with explicit treatment of termination, censoring, and host uncertainty. The formulation specifies the physical inputs and trajectory information needed for subsequent numerical investigation without selecting a cluster model or claiming calibrated transport predictions. Its immediate objective is to determine how recoil and the stellar environment jointly control the remnant's time-dependent occupation of dense cluster regions.

## Research question

How do a merger's initial conditions and the host's mass-resolved phase-space structure determine the remnant's motion and finite-time environmental exposure?

The active scope consists of two linked steps:

| Step | Physical task | Output required |
|---|---|---|
| [1. Recoil transport](docs/e001-ballistic-return.md) | Evolve a joint merger-conditioned launch under consistent host gravity, external acceleration, collisional drift, and velocity diffusion. | Timed position–velocity histories or a transport density, with mass labels, uncertainty, survival probability, and stopping causes. |
| [2. Residence and exposure](docs/e002-core-exposure.md) | Define finite-region entry and re-entry, residence, and species-resolved density exposure on those histories. | Event distributions, accumulated exposure, and a time-resolved occupation measure with explicit censoring. |

These are mathematical specifications. A numerical solver, host selection, parameter calibration, and results for this formulation are pending. No central launch, radial orbit, permanent rest after return, or universal cluster profile is imposed.

## Physical foundation

The host's species distributions determine both gravitational mass density and the local scattering environment. Initial remnant motion includes the progenitor's center-of-mass velocity and the recoil vector in a common coordinate system. Mean drift and random deflections retain the remnant's continuing motion through and beyond a declared core region.

The local weak-scattering closure has explicit applicability assumptions. It does not resolve strong encounters, evolve the host by itself, or guarantee calibrated recoil trajectories. [Mathematical reading order](docs/recoil-transport-theory.md) · [Host and data requirements](docs/cluster-dynamics-context.md) · [Primary literature](docs/references.md).

Timed position and velocity records can support later companion-access calculations, but binary acquisition, encounter rates, and repeated-merger predictions are outside the current two-step scope.

## Observed-cluster data

The [six-cluster assessment](docs/cluster-data-assessment.md) examines original public observations and published candidate host posteriors for 47 Tuc, NGC 1851, NGC 3201, Omega Centauri, M4 and NGC 6397. The [data package](data/observed-clusters/README.md) includes source hashes, common-unit tables, quality flags, three comparison figures and an [executed notebook](analysis/cluster_data/cluster_data_audit.ipynb).

A [common host interface](docs/cluster-host-schema.md) preserves a separate weighted joint parameter distribution for each cluster. Unit harmonization does not establish a shared stellar population, a calibrated potential or a local scattering background. Exact source-fit reproduction and reconstruction of mass-resolved local phase space remain necessary before transport. The data audit leaves the theoretical manuscript and physical experiments unchanged.

## Calculation handoff

The [input register](formulation/inputs.json) leaves host functions, launch distributions, collision cutoffs, stopping rules, and numerical controls unassigned. Missing inputs are not zero-valued defaults. The [handoff specification](docs/calculation-handoff.md) defines their units, dependencies, required output fields, and acceptance checks before a calculation can be interpreted.

Documentation checks perform no physical calculations:

```sh
python scripts/verify_foundation.py
python paper/verify.py
python analysis/cluster_data/verify.py
```

The [research paper](output/pdf/central-residence.pdf), [LaTeX source](paper/manuscript.tex), and [portable source archive](output/pdf/manuscript-source.zip) present the same two-step formulation. It is an unreviewed theoretical manuscript, with no numerical transport results. [Paper status and build instructions](paper/README.md).

## Repository contents

| Path | Purpose |
|---|---|
| `docs/problem.md` | Research question, scope, and mathematical objectives |
| `docs/e001-ballistic-return.md` | Step 1: recoil transport |
| `docs/e002-core-exposure.md` | Step 2: residence and environmental exposure |
| `docs/calculation-handoff.md` | Required inputs, outputs, and validation conditions |
| `formulation/` | Explicitly unassigned calculation inputs |
| `data/observed-clusters/` | Public-source manifest, candidate host records and audited comparison products |
| `analysis/cluster_data/` | Source retrieval, reproducible data audit, figures and executed notebook |
| `scripts/` | Structural documentation and manuscript integrity checks |
| `paper/` | Theory manuscript, bibliography, and calculation-free build |
| `src/`, `experiments/`, `tests/`, `results/` | [Computational reference records](docs/computational-records.md); these do not implement or validate the active transport formulation |

## Citation and reuse

The [BibTeX citation](CITATION.bib) and [CITATION.cff](CITATION.cff) identify the [versioned research formulation](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/formulation-v1), separately from an executed study or a journal publication. No DOI has been assigned.

Original code and documentation are available under the [MIT license](LICENSE). Referenced publications, dependencies, and external data retain their own rights. [Contribution guidance](CONTRIBUTING.md) · [Author and contact](AUTHORS.md).
