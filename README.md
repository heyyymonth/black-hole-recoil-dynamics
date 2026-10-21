# Black-hole recoil dynamics

Independent physics research · Stellar dynamics and gravitational-wave source environments

**Contributed by Hemanth.**

[![Numerical verification](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml/badge.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/Version-0.1.1-14677d.svg)](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.1)

[E001–E004 paper](output/pdf/central-residence.pdf) · [LaTeX source](paper/manuscript.tex) · [Source archive](output/pdf/manuscript-source.zip) · [E004 cluster comparison](docs/e004-cluster-validation.md) · [Research context](docs/cluster-dynamics-context.md) · [Baseline release](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.1)

## Research question

How does a retained black-hole remnant's motion affect its access to potential companions, and how accurately does density exposure represent that access?

E001/E002 follow a test particle in a fixed Plummer potential and measure return, residence, and density exposure. E003 adds a prescribed tracer-binary population and velocity-dependent close-passage cross section. E004 compares those population and potential assumptions with an externally evolved CMC cluster snapshot. Together these test ingredients connecting recoil motion to companion access; they do not yet predict repeated-merger rates.

## Experiments

| Experiment | Calculation | Research record |
|---|---|---|
| E001 | First central return as a function of kick speed | [Method and results](docs/e001-ballistic-return.md) · [Table](results/e001/return_times.csv) |
| E002 | Central residence and density exposure over finite windows | [Method and results](docs/e002-core-exposure.md) · [Tables](results/e002/) |
| E003 | Density proxy versus velocity-weighted close passages over a common time window | [Protocol](docs/e003-protocol.md) · [Method and results](docs/e003-encounter-accessibility.md) · [Table](results/e003/encounters.csv) |
| E004 | Companion inventory and ballistic returns in an evolved cluster snapshot | [Protocol](docs/e004-protocol.md) · [Method and results](docs/e004-cluster-validation.md) · [Tables](results/e004/) |

**Latest result:** E004 finds only **two BH–BH binaries** among 682,583 systems in one evolved cluster snapshot: insufficient to calibrate a continuous companion distribution. At identical physical kicks, a Plummer model matched to current host mass and half-mass radius predicts returns **2.96, 2.81, and 4.45 times as long** as the frozen snapshot potential. Direct full-shell checks place compression errors below 0.16%. These results diagnose assumptions; **live recoil trajectories remain unvalidated**. [Results and limitations](docs/e004-cluster-validation.md) · [Specific next experiments and existing research](docs/cluster-dynamics-context.md).

![External cluster comparison](results/e004/snapshot_comparison.png)

In E003, density exposure overestimates normalized close-passage opportunities by **4.4%, 53.1%, and 98.6%** at kick fractions 0.3, 0.8, and 0.9. Two cases fail the predeclared 10% approximation screen. Gravitational focusing dominates this specified configuration; faster relative motion reduces its contribution. Absolute expected passage counts are tiny: 1.61e-6 for the stationary control and 4.71e-8 for the 0.9 kick during the common 0.149 Myr window.

![Density proxy compared with velocity-weighted close passages](results/e003/proxy_comparison.png)

These are controlled numerical benchmarks. E003 calculates an expected close-passage count in an imposed bath along an unperturbed orbit; it does not simulate encounter outcomes or infer mergers. All four experiments omit dynamical friction and host evolution. E003 screens remnant apocenters against a declared outer radius without truncating the host. The former E002 near-escape example remains an isolated-model limit and does not pass that screen. No observed gravitational-wave data are used.

The [unreviewed manuscript](output/pdf/central-residence.pdf), **From recoil motion to companion access: controlled benchmarks and a cluster-snapshot comparison**, now integrates **E001–E004**. It develops the connection from orbital exposure to velocity-weighted passages and then assesses transfer of the assumed population and potential to an external cluster. The [claim-to-evidence record](paper/CLAIMS.md) explains the resulting scope. It makes no calibrated merger-rate or live-trajectory claim. The archived v0.1.1 release remains the earlier E001/E002 paper.

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
python experiments/verify_e004.py
```

The **108 tests** compare production calculations with independent trajectories, numerical velocity integrals, conservation laws, analytical limits, and scaling checks. The verifiers separately reproduce **146 baseline rows, 4 E003 rows, and E004's 5 population rows, 12 trajectory cases, and 160 profile points**, and check recorded source/output hashes. These rows are model calculations, not independent astronomical observations.

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

For E004, `python experiments/fetch_e004.py` downloads approximately 46 MB of pinned external inputs to ignored local cache; the verifier also fetches missing inputs. After committing source changes, use `python experiments/run_e004.py` to regenerate. The [E004 manifest](results/run_e004.json) records external and local input hashes. Subsequent checks work offline with the populated cache.

Each runner replaces its corresponding outputs. E001–E003 need no network after dependency installation. The [baseline configuration](experiments/config.json) and [E003 configuration](experiments/e003.json) fix inputs; [baseline](results/run_all.json) and [E003](results/run_e003.json) manifests record source revisions, hashes, and library versions. [Paper build instructions](paper/README.md) cover the separate PDF dependencies and source.

## Repository structure

| Path | Purpose |
|---|---|
| `src/recoil_dynamics/` | Reusable orbit, close-passage, snapshot-parsing, and shell-potential calculations |
| `experiments/` | Configuration, execution, result verification, and record template |
| `tests/` | Independent numerical and analytical validation |
| `results/` | Generated tables, figures, provenance, and verification record |
| `docs/` | Physical problem, experiment records, and literature context |
| `paper/` | LaTeX manuscript, BibTeX, vector figures, and build script |
| `output/pdf/` | Working paper, source archive, and build manifest |

[Contribution guidance](CONTRIBUTING.md) describes how to add a separately identified experiment with its own assumptions, configuration, validation, and results.

## Citation and reuse

For the current **E001–E004 manuscript**, use the [paper-specific BibTeX citation](paper/CITATION.bib), identifying [immutable manuscript snapshot 0b2767e](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/0b2767ebe455eb891bffddb03a7e65cc442539a7). Revision `0.2.0-manuscript.1` includes the expanded eight-page paper and supporting source; it is an unreviewed manuscript with no DOI. The root citation files below continue to identify the earlier baseline software release.

Athota, Venkata Hemanth (2026). *Black-hole recoil dynamics: central residence and density exposure in a Plummer potential*. Version 0.1.1. [Versioned research release](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.1).

Use [CITATION.cff](CITATION.cff), [BibTeX](CITATION.bib), or GitHub's **Cite this repository** control. Cite the version used for the calculation. This release has no DOI; the versioned release identifies the research snapshot.

Those citation files describe the E001/E002 v0.1.1 baseline. For E003, use its [BibTeX citation](results/e003/CITATION.bib), which identifies [results snapshot 7944b82](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/7944b82b3c99e070b7c312011020e4ee53b4255d). The [E003 manifest](results/run_e003.json) also identifies the generating source commit. This commit-based citation has no DOI; E003 is not part of v0.1.1.

For E004, use its [BibTeX citation](results/e004/CITATION.bib), identifying [results snapshot 14d9f49](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/14d9f493e6eeae25e51a244dfb60789753992db6). Its [manifest](results/run_e004.json) distinguishes generating code from external simulation provenance. Cite the upstream CMC/cmctoolkit sources as appropriate when using those data. This snapshot also contains the working abstract revision; it has no DOI and does not replace v0.1.1.

Original code, technical documentation, and generated research artifacts are provided under the [MIT license](LICENSE). Preserve its copyright and permission notice when reusing substantial portions. Academic citation is requested separately from the license conditions. Referenced publications, dependencies, and external simulation files retain their own rights and licenses. The fetched CMC snapshot and upstream summaries are not redistributed or relicensed by this repository.

[Author and contact](AUTHORS.md) · [Release record](CHANGELOG.md)
