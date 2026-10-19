# Black-hole recoil dynamics

Independent physics research · Stellar dynamics and gravitational-wave source environments

**Contributed by Hemanth.**

[![Numerical verification](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml/badge.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-0.1.1-14677d.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.1)

[E001/E002 paper](output/pdf/central-residence.pdf) · [LaTeX source](paper/manuscript.tex) · [Source archive](output/pdf/manuscript-source.zip) · [E003 research record](docs/e003-encounter-accessibility.md) · [Baseline release](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.1)

## Research question

How does a retained black-hole remnant's motion affect its access to potential companions, and how accurately does density exposure represent that access?

The experiments follow a test particle kicked radially from the center of a fixed Plummer potential. E001/E002 measure return, residence, and density exposure. E003 adds an explicit tracer-binary population and velocity-dependent close-passage cross section. This tests the connection between orbital motion and encounter opportunity, an ingredient in the larger black-hole assembly problem.

## Experiments

| Experiment | Calculation | Research record |
|---|---|---|
| E001 | First central return as a function of kick speed | [Method and results](docs/e001-ballistic-return.md) · [Table](results/e001/return_times.csv) |
| E002 | Central residence and density exposure over finite windows | [Method and results](docs/e002-core-exposure.md) · [Tables](results/e002/) |
| E003 | Density proxy versus velocity-weighted close passages over a common time window | [Protocol](docs/e003-protocol.md) · [Method and results](docs/e003-encounter-accessibility.md) · [Table](results/e003/encounters.csv) |

In E003, density exposure overestimates normalized close-passage opportunities by **4.4%, 53.1%, and 98.6%** at kick fractions 0.3, 0.8, and 0.9. Two cases fail the predeclared 10% approximation screen. Gravitational focusing dominates this specified configuration; faster relative motion reduces its contribution. Absolute expected passage counts are tiny: 1.61e-6 for the stationary control and 4.71e-8 for the 0.9 kick during the common 0.149 Myr window.

![Density proxy compared with velocity-weighted close passages](results/e003/proxy_comparison.png)

These are controlled numerical benchmarks. E003 calculates an expected close-passage count in an imposed bath along an unperturbed orbit; it does not simulate encounter outcomes or infer mergers. All three experiments omit dynamical friction and host evolution. E003 screens remnant apocenters against a declared outer radius without truncating the host. The former E002 near-escape example remains an isolated-model limit and does not pass that screen. No observed gravitational-wave data are used.

The [unreviewed manuscript](output/pdf/central-residence.pdf) and v0.1.1 release cover E001/E002. E003 is documented separately on `main`; it has not been added to that paper or release.

[Physical assumptions](docs/problem.md) · [Scientific relevance](docs/relevance.md) · [References](docs/references.md)

## Reproduce the results

Python **3.11** is the supported reference environment. From the repository root:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pytest -q
python experiments/verify_results.py
python experiments/verify_e003.py
```

The **62 tests** compare production calculations with independent trajectories, numerical velocity integrals, conservation laws, analytical limits, and scaling checks. The verifiers separately reproduce **146 baseline rows and 4 E003 rows**, and check recorded source/output hashes. These rows are model calculations, not independent astronomical observations.

To regenerate E001/E002 tables and figures:

```sh
python experiments/run.py
python experiments/verify_results.py
```

For E003, commit its source inputs, then run:

```sh
python experiments/run_e003.py
python experiments/verify_e003.py
```

Each runner replaces its corresponding outputs. No network access is needed after dependency installation. The [baseline configuration](experiments/config.json) and [E003 configuration](experiments/e003.json) fix inputs; [baseline](results/run_all.json) and [E003](results/run_e003.json) manifests record source revisions, hashes, and library versions. [Paper build instructions](paper/README.md) cover the separate PDF dependencies and source.

## Repository structure

| Path | Purpose |
|---|---|
| `src/recoil_dynamics/` | Reusable orbit, velocity-moment, and close-passage calculations |
| `experiments/` | Configuration, execution, result verification, and record template |
| `tests/` | Independent numerical and analytical validation |
| `results/` | Generated tables, figures, provenance, and verification record |
| `docs/` | Physical problem, experiment records, and literature context |
| `paper/` | LaTeX manuscript, BibTeX, vector figures, and build script |
| `output/pdf/` | Published paper and build manifest |

[Contribution guidance](CONTRIBUTING.md) describes how to add a separately identified experiment with its own assumptions, configuration, validation, and results.

## Citation and reuse

Athota, Venkata Hemanth (2026). *Black-hole recoil dynamics: central residence and density exposure in a Plummer potential*. Version 0.1.1. [Versioned research release](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.1).

Use [CITATION.cff](CITATION.cff), [BibTeX](CITATION.bib), or GitHub's **Cite this repository** control. Cite the version used for the calculation. This release has no DOI; the versioned release identifies the research snapshot.

Those citation files describe the E001/E002 v0.1.1 baseline. For E003, use its [BibTeX citation](results/e003/CITATION.bib), which identifies [results snapshot 7944b82](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/7944b82b3c99e070b7c312011020e4ee53b4255d). The [E003 manifest](results/run_e003.json) also identifies the generating source commit. This commit-based citation has no DOI; E003 is not part of v0.1.1.

Original code, technical documentation, and generated research artifacts are provided under the [MIT license](LICENSE). Preserve its copyright and permission notice when reusing substantial portions. Academic citation is requested separately from the license conditions. Referenced publications and dependencies retain their own rights and licenses.

[Author and contact](AUTHORS.md) · [Release record](CHANGELOG.md)
