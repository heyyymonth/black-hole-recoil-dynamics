# Black-hole recoil dynamics

**Hemanth AV Research** · Researcher: **Hemanth AV**

How well do a retained/not-retained label and a first-return time represent a recoiling remnant's access to the dense central region of a stellar cluster?

We study this question in a controlled model: a test particle kicked radially from the center of a static, isolated Plummer potential. We calculate the first central return, finite-time residence inside specified spheres, and density integrated along the orbit. These are orbital quantities; no companion-acquisition rate or repeated-merger probability is inferred.

| Experiment | Question | Record |
|---|---|---|
| E001: ballistic return | How does the first central-return time depend on kick speed? | [Method and results](docs/e001-ballistic-return.md) · [Numerical table](results/e001/return_times.csv) |
| E002: central residence and density exposure | What information is lost by replacing continuous orbital motion with an assumed central residence or a first-return gate? | [Method and results](docs/e002-core-exposure.md) · [Numerical tables](results/e002/) |

![Finite-horizon density exposure](results/e002/exposure_comparison.png)

The comparison uses two explicitly defined approximations: an object assumed always central, and an object assumed permanently central after its first return. The latter is not a dynamical-friction model. In the collisionless calculation, the remnant continues oscillating after each return.

[Problem and physical assumptions](docs/problem.md) · [Scientific relevance review](docs/relevance.md) · [References](docs/references.md)

## Reproduce

Python 3.11 is the reference environment. Run from the repository root:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pytest -q
python experiments/run.py
python experiments/verify_results.py
```

The experiments require no network access after dependency installation. Direct dependencies are specified in `requirements.txt`; `requirements-lock.txt` pins the tested environment including transitive packages. The run produces CSV tables, PNG/SVG figures, and `results/run_all.json`, which records parameters through a configuration hash, source revision and hashes, output hashes, and numerical-library versions. Running the script replaces the corresponding generated results.

The tests compare energy quadrature against independently integrated equations of motion and analytical limits. The result verifier checks the numerical tables against the current implementation and detects changed input or output files. The repository's verification workflow runs both checks.

## Scope

The model excludes dynamical friction, encounters, tidal truncation, and evolution of the host potential. Its density weight follows a Plummer profile; a real black-hole population can have a different spatial distribution. Near-escape excursions are particularly sensitive to omitted host boundaries. These experiments are reproducible benchmarks, not observational constraints or a claim of a new formation channel.

This repository contains original code, generated experiment results and source citations. Hemanth AV Research identifies this personal research collection; it does not denote an institutional affiliation.
