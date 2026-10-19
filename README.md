# Black-hole recoil dynamics

Independent physics research · Stellar dynamics and gravitational-wave source environments

**Contributed by Hemanth.**

[![Numerical verification](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml/badge.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-0.1.0-14677d.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.0)

[Research paper](output/pdf/central-residence.pdf) · [Paper source](paper/manuscript.md) · [Cite this work](CITATION.cff) · [Release](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.0)

## Research question

How accurately do simplified retention and return prescriptions represent a recoiling remnant's exposure to the dense center of a stellar cluster?

The experiments follow a test particle kicked radially from the center of a static, isolated Plummer potential. They measure first central return, finite-time residence inside declared spheres, and density integrated along the orbit. This isolates the orbital part of the black-hole assembly problem.

## Published experiments

| Experiment | Calculation | Research record |
|---|---|---|
| E001 | First central return as a function of kick speed | [Method and results](docs/e001-ballistic-return.md) · [Table](results/e001/return_times.csv) |
| E002 | Central residence and density exposure over finite windows | [Method and results](docs/e002-core-exposure.md) · [Tables](results/e002/) |

At 99% of the central escape speed, the bound orbit samples a mean density of **0.130%** of the central value and spends **0.193%** of its radial period inside one Plummer scale radius. A prescription that makes the object permanently central after its first return misses early exposure and eventually overestimates later exposure.

![Density exposure compared with simplified prescriptions](results/e002/exposure_comparison.png)

These are controlled numerical benchmarks. The model omits dynamical friction, encounters, tidal boundaries, and host evolution. It does not infer encounter rates, merger probabilities, or constraints from observed gravitational-wave events. The [technical note](output/pdf/central-residence.pdf) records the derivation, results, independent checks, and limitations; it is not peer reviewed.

[Physical assumptions](docs/problem.md) · [Scientific relevance](docs/relevance.md) · [References](docs/references.md)

## Reproduce the results

Python **3.11** is the supported reference environment. From the repository root:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pytest -q
python experiments/verify_results.py
```

The tests compare the production quadrature with independently integrated trajectories, conservation laws, analytic limits, and convergence checks. The verifier separately reproduces **146 recorded table rows** and checks source/output hashes.

To regenerate every table and figure:

```sh
python experiments/run.py
python experiments/verify_results.py
```

The runner replaces the corresponding generated outputs. No network access is needed after dependency installation. [Configuration](experiments/config.json) fixes the grids and tolerances; the [run manifest](results/run_all.json) records source revision, file hashes, and numerical-library versions. [Paper build instructions](paper/README.md) cover the separate PDF dependencies and source.

## Repository structure

| Path | Purpose |
|---|---|
| `src/recoil_dynamics/` | Reusable dimensionless orbital calculations |
| `experiments/` | Configuration, execution, result verification, and record template |
| `tests/` | Independent numerical and analytical validation |
| `results/` | Generated tables, figures, provenance, and verification record |
| `docs/` | Physical problem, experiment records, and literature context |
| `paper/` | Editable technical note and PDF build script |
| `output/pdf/` | Published paper and build manifest |

[Contribution guidance](CONTRIBUTING.md) describes how to add a separately identified experiment with its own assumptions, configuration, validation, and results.

## Citation and reuse

Athota, Venkata Hemanth (2026). *Black-hole recoil dynamics: central residence and density exposure in a Plummer potential*. Version 0.1.0. [Versioned research release](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.0).

Use [CITATION.cff](CITATION.cff), [BibTeX](CITATION.bib), or GitHub's **Cite this repository** control. Cite the version used for the calculation. This release has no DOI; the versioned release identifies the research snapshot.

Original code, technical documentation, and generated research artifacts are provided under the [MIT license](LICENSE). Preserve its copyright and permission notice when reusing substantial portions. Academic citation is requested separately from the license conditions. Referenced publications and dependencies retain their own rights and licenses.

[Author and contact](AUTHORS.md) · [Release record](CHANGELOG.md)
