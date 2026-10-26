# Recoil transport and environmental exposure of black-hole merger remnants

Independent physics research · Stellar dynamics and gravitational-wave source environments

**Contributed by Hemanth.**

[![Research verification](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml/badge.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[Problem statement](docs/problem.md) · [Step 1](docs/step1-transport.md) · [Step 2](docs/step2-exposure.md) · [Calculation handoff](docs/calculation-handoff.md)

## Abstract

Black-hole mergers in stellar clusters can leave remnants moving through their hosts after receiving gravitational-wave recoil. Remaining gravitationally bound does not determine how soon a remnant enters a dense region or how much time it spends there. Predicting this motion requires a stellar environment that consistently supplies gravity, collisional drift, and velocity diffusion. We formulate merger-conditioned transport and define finite-region residence and species-resolved exposure on its trajectories. Public observations and published multimass-model posteriors for six Galactic clusters establish the available constraints and their shared uncertainties. For one joint realization of NGC 3201, we verify the mass distribution, potential, velocity moments, and projected tracer profiles, then evaluate anisotropic scattering from the same species distributions. A sparse black-hole component dominates the formal central diffusion, identifying a regime where the smooth closure requires a discrete-encounter test. Galactic reference orbits constructed from public phase-space coordinates supply the time-dependent external field. An analytic force bound separates a sufficient accuracy domain for isolated controls from a diagnostic tidal scale. Disk passages restrict this domain and produce rapid changes in forcing, identifying phases where the frozen-host approximation needs additional validation. These calculations establish conditional host, scattering, and environmental diagnostics; consistent host evolution and calibrated remnant residence distributions remain unresolved.

## Research question

How do a merger's initial conditions and the host's mass-resolved phase-space structure determine the remnant's motion and finite-time environmental exposure?

The active scope consists of two linked steps:

| Step | Physical task | Output required |
|---|---|---|
| [1. Recoil transport](docs/step1-transport.md) | Evolve a joint merger-conditioned launch under consistent host gravity, external acceleration, collisional drift, and velocity diffusion. | Timed position–velocity histories or a transport density, with mass labels, uncertainty, survival probability, and stopping causes. |
| [2. Residence and exposure](docs/step2-exposure.md) | Define finite-region entry and re-entry, residence, and species-resolved density exposure on those histories. | Event distributions, accumulated exposure, and a time-resolved occupation measure with explicit censoring. |

The observational analysis constrains candidate hosts. A completed NGC 3201 reconstruction now supplies consistent species distributions, and twenty evaluated local states test its anisotropic drift and diffusion. Four Galactic reference orbits now supply time-dependent tides and a [host-support audit](docs/galactic-orbit.md). Remnant trajectory evolution and residence distributions remain unevaluated; the sparse black-hole component and the host in a Galactic tide require physical validation. A [boundary calculation](docs/outer-boundary.md) gives a conditional tidal-opening scale and a separate accuracy limit for isolated controls. No central launch, radial orbit, permanent rest after return, or universal cluster profile is imposed.

## Physical foundation

The host's species distributions determine both gravitational mass density and the local scattering environment. Initial remnant motion includes the progenitor's center-of-mass velocity and the recoil vector in a common coordinate system. Mean drift and random deflections retain the remnant's continuing motion through and beyond a declared core region.

The local weak-scattering closure has explicit applicability assumptions. It does not resolve strong encounters, evolve the host by itself, or guarantee calibrated recoil trajectories. [Mathematical reading order](docs/recoil-transport-theory.md) · [Host and data requirements](docs/cluster-dynamics-context.md) · [Primary literature](docs/references.md).

Timed position and velocity records can support later companion-access calculations, but binary acquisition, encounter rates, and repeated-merger predictions are outside the current two-step scope.

## Observed-cluster data

The [six-cluster assessment](docs/cluster-data-assessment.md) examines original public observations and published candidate host posteriors for 47 Tuc, NGC 1851, NGC 3201, Omega Centauri, M4 and NGC 6397. The [data package](data/observed-clusters/README.md) includes source hashes, common-unit tables, quality flags, three comparison figures and an [executed notebook](analysis/cluster_data/cluster_data_audit.ipynb).

A [common host interface](docs/cluster-host-schema.md) preserves a separate weighted joint parameter distribution for each cluster. Unit harmonization does not establish a shared stellar population, a calibrated potential or a local scattering background. The [fixed-host calculation](docs/local-transport-calculations.md) reconstructs one NGC 3201 joint realization; exact historical fit reproduction and physical calibration are not claimed. The manuscript presents the data findings and the completed host and collision calculations in seven tables and seven vector figures. It records 107 original HST radial bins, 621,878 Gaia field-catalog rows, 20 flagged compilation rows and strong mass–radius dependence in several published host posteriors. Field-catalog rows include nonmembers, and the display uncertainties condition on catalog point scales.

## Calculation handoff

The general [trajectory input register](formulation/inputs.json) leaves the launch law, physical collision cutoffs, stopping rules and trajectory controls unassigned. The completed host and local-coefficient calculations have separate versioned protocols, linked from that register. Missing inputs are not zero-valued defaults. The [handoff specification](docs/calculation-handoff.md) defines their units, dependencies, required output fields, and acceptance checks before a calculation can be interpreted.

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
```

The [research paper](output/pdf/central-residence.pdf), [LaTeX source](paper/manuscript.tex), and [portable source archive](output/pdf/manuscript-source.zip) combine the two-step formulation with the observed-cluster findings. It is an unreviewed manuscript reporting data analysis, a reconstructed host, and local collision kernels. It does not report remnant trajectories or calibrated residence predictions. [Paper status and build instructions](paper/README.md) · [Calculation audit](docs/repository-audit.md).

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
| `scripts/` | Structural documentation and manuscript integrity checks |
| `paper/` | Formulation and observational findings, evidence-bound tables and figures, bibliography, and typesetting build |
| `tests/` | Regression checks for input validity, configuration handling and evidence consistency |
| [Historical records](docs/computational-records.md) | Immutable source and results for earlier studies |

## Citation and reuse

The [BibTeX citation](CITATION.bib) and [CITATION.cff](CITATION.cff) identify the [versioned research formulation](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/formulation-v1), separately from an executed study or a journal publication. No DOI has been assigned.

Original code and documentation are available under the [MIT license](LICENSE). Referenced publications, dependencies, and external data retain their own rights. [Contribution guidance](CONTRIBUTING.md) · [Author and contact](AUTHORS.md).
