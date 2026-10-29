# Research manuscript

**Recoil transport and environmental exposure of black-hole merger remnants**

Venkata Hemanth Athota · Independent researcher<br/>
bharadwaj.hemanth.002@gmail.com

[PDF](../output/pdf/central-residence.pdf) · [LaTeX source](manuscript.tex) · [Portable source](../output/pdf/manuscript-source.zip) · [Citation](CITATION.bib)

The manuscript studies when continuous encounter dynamics preserve merger-remnant residence and environmental exposure. It connects this hypothesis to calculations in six adopted cluster environments. It reports source coverage, native-to-physical conversions, joint host-posterior summaries, structural-estimate comparisons and their implications for the transport inputs. Nine tables and nine vector figures present the data audit, a fixed NGC 3201 host reconstruction, twenty local anisotropic collision states, four static boundary scenarios, and four Galactic reference orbits with time-dependent host-support diagnostics. A six-cluster comparison adds six reference hosts, 120 local states, a 48-row sensitivity pilot and per-cluster environmental controls, with model-suitability decisions retained separately. The central hypothesis concerns observable accuracy: a matched finite-encounter and diffusion comparison is connected to residence through a backward occupation equation and its nonlocal generator residual. A separate diffusion-based response predictor is proposed for testing. Finite-angle calculations add 162 local settings and 600,000 independently sampled increments, with an explicit encounter domain and a discrete single-transition response estimator verified against exact finite-state calculations. Finite-population response and host evolution have distinct controls. No new cluster fit, calibrated friction coefficient, remnant trajectory or residence result is claimed.

Revision **0.11.0-manuscript.1** is an unreviewed formulation, data-analysis and local-calculation manuscript. No submission, acceptance, or DOI is claimed. The [claim-to-evidence record](CLAIMS.md) distinguishes established theory, completed data calculations, proposed physical tests and remaining calibration requirements. [working.json](working.json) records the title, abstract, data revision and scope.

## Build and verify

The article uses the two-column APS/PRD REVTeX format with its standard typography. See [FORMAT.md](FORMAT.md).

With Python 3.11, the paper dependencies, and Tectonic 0.16.9 installed:

```sh
python -m pip install -r paper/requirements.txt
python scripts/verify_foundation.py
python paper/data_assets.py --verify
python analysis/host_reconstruction/verify.py
python paper/calculation_assets.py --verify
python analysis/boundary/verify.py
python paper/boundary_assets.py --verify
python analysis/orbit/verify.py
python paper/orbit_assets.py --verify
python analysis/cluster_comparison/verify.py
python paper/comparison_assets.py --verify
python analysis/encounters/verify.py
python paper/encounter_assets.py --verify
python paper/build.py
python paper/verify.py
```

The builder checks the recorded data revision and evidence assets, then typesets LaTeX and BibTeX. It does not retrieve observations, rerun their analysis, generate plots, fit hosts, run experiments or evaluate transport equations. It writes the PDF, portable source ZIP, and build manifest under `output/pdf/`.

To explicitly regenerate the presentation assets, use the Python environment described in the [data package](../data/observed-clusters/README.md) and run `python paper/data_assets.py`. This reads existing tables and posterior diagnostics, generates native LaTeX tables and publication-sized vector figures, and records their source/output hashes. The read-only `--verify` mode requires only the Python standard library. The immutable data audit and its original comparison figures remain separate from these paper layouts.

Hashes bind the article, bibliography, formulation, handoff, source audit, evidence tables, figure-generation inputs and build metadata. Verification checks citations and archive consistency, page dimensions, embedded fonts, and unresolved references or overfull equations. The archive includes the tables, figures, supporting evidence and pinned unmodified REVTeX distribution; other standard TeX packages remain external dependencies.

Render and inspect every page after a layout change. Source hashes and successful compilation do not substitute for visual inspection or scientific validation.

Use [manuscript-v0.11.0.1](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/manuscript-v0.11.0.1) for an exact citation. The data analysis is pinned to [cluster-data-v0.1.0](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/cluster-data-v0.1.0). Computational reference records have separate [immutable citations](../docs/computational-records.md).

The [host and local-kernel reproduction instructions](../analysis/host_reconstruction/README.md) specify the scientific environment. `python paper/calculation_assets.py` regenerates only the presentation of recorded calculations; it never reruns the model or collision integrals.

The [outer-boundary calculation](../analysis/boundary/README.md) has its own source and numerical record. `python paper/boundary_assets.py` presents those recorded results without evolving any orbit.

The [Galactic reference-orbit stage](../analysis/orbit/README.md) supplies external forcing and host-validity diagnostics. It does not evolve a merger remnant or re-equilibrate the stellar host. `python paper/orbit_assets.py` generates the presentation from its recorded results.

The [six-cluster comparison](../analysis/cluster_comparison/README.md) retains 47 accepted pilot hosts and one unresolved M4 reconstruction. `python paper/comparison_assets.py` presents the recorded comparison; it does not rerun the pilot.
