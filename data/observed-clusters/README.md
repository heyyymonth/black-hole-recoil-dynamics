# Observed-cluster data package

Venkata Hemanth Athota · Data revision `cluster-data-v0.1.0` · 24 October 2026

Original public source products and an executed audit of six candidate stellar-cluster environments. The package standardizes records and units, describes existing joint host posteriors, and identifies unresolved inputs for recoil transport. It contains no fitted host, evaluated potential, friction coefficient, remnant trajectory or residence calculation.

[Scientific assessment](../../docs/cluster-data-assessment.md) · [Common host schema](../../docs/cluster-host-schema.md) · [Executed notebook](../../analysis/cluster_data/cluster_data_audit.ipynb) · [Candidate records](candidate-hosts.json)

## Source capture and attribution

The [source manifest](sources.json) records 51 original files (269,365,248 bytes), retrieval URLs, byte counts, SHA256 hashes, capture time, roles and source rights. Large originals are retrieved into `data/observed-clusters/raw/`, which is excluded from Git. The repository contains derived tables, figures, an executed notebook and the manifest needed to verify the captured originals.

| Source | Products used | Attribution and rights |
|---|---|---|
| [Baumgardt cluster database](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/parameter.html) and [profile definitions](https://people.smp.uq.edu.au/HolgerBaumgardt/globular/newdata/indivdata/velsurf.html) | Two separate structural-table paths; six compiled kinematic and six surface profiles | Cite the database and the studies identified for each series. Provider rights apply; no blanket MIT relicensing. Compiled quantities remain attributed to their original studies. |
| [MAST HACKS](https://archive.stsci.edu/hlsp/hacks), DOI [10.17909/jpfd-2m08](https://doi.org/10.17909/jpfd-2m08) | Original HST dispersion bins for five clusters | Libralato et al. (2022), DOI [10.3847/1538-4357/ac7727](https://doi.org/10.3847/1538-4357/ac7727); CC BY 4.0. |
| [oMEGACat VI release](https://zenodo.org/records/14978551) | Omega Centauri HST dispersion bins, LOS profile and mass-dependent supporting products | Häberle et al. (2025), [source study](https://arxiv.org/abs/2503.04903), DOI 10.5281/zenodo.14978551; CC BY 4.0. |
| [Gaia EDR3 v2 release](https://zenodo.org/records/4891252) | Original stellar catalogs, membership/quality fields and fitted dispersion/rotation curves | Vasiliev & Baumgardt (2021), DOI [10.1093/mnras/stab1475](https://doi.org/10.1093/mnras/stab1475), data DOI 10.5281/zenodo.4891252; CC BY 4.0. Provider archive MD5 is checked in addition to our SHA256. |
| [GCfit results at the captured revision](https://github.com/nmdickson/GCfit-results/tree/c45a03b1681570b6ee96be76f362894c6628e1f1) | Six standard-fit chains and summaries, Gaia reanalysis products, NGC 6397 alternative summary | Dickson et al., [source study](https://arxiv.org/abs/2308.13037). Preserve release attribution and provider terms. No alternative NGC 6397 no-BH chain was found. |
| [Historical GCfit observation resources](https://github.com/nmdickson/GCfit/tree/ca63e5b465d01bb5b006dce7d40e28d0a95643ed/gcfit/resources) | Six candidate observation bundles, including mass-function information | Preserve per-dataset citations and repository terms. Exact pairing to the old posterior runs is unverified. |

Changes made to source products are limited to extraction, source/quality annotation, declared unit conversions, weighted posterior summaries, and visualization. CC BY source attribution is retained in this file and the source manifest. The extracted `gaia_profile_README.txt` is the original archive documentation. The repository MIT license covers our original analysis code and documentation; external material retains its own rights.

Some database URLs are mutable. Re-fetching succeeds only when the original bytes still match the captured hashes. If an upstream product changes or disappears, a fresh download cannot reproduce this capture without the matching original; do not weaken the hash check or silently update this data revision. GitHub sources use immutable commit URLs and Zenodo sources use version-specific record identifiers.

## Reproduce the audit

Run from the repository root. The recorded environment is Python 3.11 on macOS; the lock file pins the analysis and notebook dependencies. `curl` is required for retrieval. Raw source capture needs about 270 MB, in addition to the Python environment.

```sh
python3.11 -m venv .venv-data
.venv-data/bin/python -m pip install -r analysis/cluster_data/requirements-lock.txt
.venv-data/bin/python analysis/cluster_data/fetch.py
.venv-data/bin/python analysis/cluster_data/run_notebook.py
.venv-data/bin/python analysis/cluster_data/verify.py
```

The notebook executes `analyze.py` and `figures.py`, verifies every raw-file SHA256, checks source schemas and keys, and recreates the derived products. It uses the interpreter invoking `run_notebook.py`. The analysis and figure scripts can also run separately. None imports or runs the recoil experiments.

For a check of the versioned artifacts without downloading sources or installing analysis dependencies:

```sh
python3 analysis/cluster_data/verify.py
```

The verification script checks the released artifact hashes, coverage totals, source links, candidate-record consistency, execution-disabled state and completed notebook outputs. It is an integrity and consistency check, not an independent reproduction of the astrophysical models. Re-execution with another platform or library set can change serialization, plot fonts or floating-point details. Review any differences before creating an explicitly revised lock with `verify.py --write-lock`; that option does not certify their scientific correctness.

## Records and meanings

| Product | Record grain and contents | Limits |
|---|---|---|
| `derived/catalog_scalars.csv` | One row per cluster; one selected UQ table path | Model estimates, not direct measurements. Radius and escape-speed errors are missing. |
| `derived/root_catalog_comparison.csv` | One row per cluster from the other UQ table path | Kept separate; not combined with selected values or interpreted as evolution. |
| `derived/hst_native_and_common_units.csv` | Cluster × source × radial bin × component; 321 rows from 107 bins | Combined, radial and tangential quantities share stars. Combined dispersion is not a third independent component. |
| `derived/compiled_kinematics_audit.csv` | 418 original compilation rows with source line and quality flags | Audit only. Twenty rows are quarantined; unflagged rows are not automatically approved likelihood inputs. |
| `derived/surface_profiles.csv` | 422 native projected surface-profile records | Stitched count/brightness normalization is not a 3D mass density. |
| `derived/gaia_correlated_curves.csv` | Six fitted profiles × 101 grid points | Correlated fitted curves; no full covariance. Column 12 is retained without interpretation. |
| `derived/gaia_star_catalog_quality.csv` | Six summaries of 621,878 original field-catalog rows | Field rows include nonmembers. High-probability counts are diagnostics, not a fitted selection. |
| `derived/posterior_parameters.csv` | Six clusters × 13 parameters | Direct weighted empirical moments and marginal quantiles of published standard fits; not newly fitted values. |
| `derived/posterior_diagnostics.json` | Sample/weight diagnostics, prior and exclusion metadata, covariance and correlation matrices | Covariances use explicitly labeled native units; mass needs a factor of one million per mass axis to match the solar-mass table. |
| `derived/model_product_inventory.json` and `observation_bundle_inventory.json` | Available derived fields and observation resources | Inventories do not establish exact posterior/data/software compatibility. |
| `derived/coverage.csv` and `quality.json` | Sampling coverage, flags and audit outcomes | Median-radius spans do not supply missing annular edges or an approved common likelihood. |
| `candidate-hosts.json` | One common interface, six separate candidate host records | Local species distributions and all launch/transport choices remain unassigned; execution is disabled. |

`gaia_archive_members.json` records the exact extracted curve members and hashes; the catalog-quality table records catalog-member hashes. Figure pairs in `derived/` provide PNG previews and vector PDF copies.

## Units and uncertainty conventions

- Canonical identifiers retain four digits, for example `NGC0104`. Native radii are angular; proper motions use mas/yr. Radial and tangential PM refer to projected components, not 3D velocity components.
- Physical display views use pc, km/s, solar masses and Myr. For distance D in kpc, transverse speed is `4.740470463533349 × D × proper_motion_masyr`. Angular radii in arcsec are converted by `1000 × D / 206264.80624709636`.
- `measurement_error_*` fields and plotted bars condition on the chosen distance. Shared distance uncertainty and distance derivatives are retained separately; the common contribution to covariance has not been added independently to each row.
- Display coordinates `R/rh` and `sigma/sqrt(G M/rh)` condition on one catalog's point estimates. They do not propagate the unavailable joint mass/radius/distance uncertainty or establish homologous hosts. The time scale `rh/sqrt(G M/rh)` is converted to Myr; it is not a dynamical-friction time.
- HACKS errors are the authors' symmetric summaries; H25 retains lower and upper offsets separately. Missing HACKS annular boundaries are empty fields, not zero radii.
- Posterior quantiles use linear interpolation in the sorted weighted cumulative distribution. The reported covariance is the weighted empirical covariance, without the finite-weight unbiased-estimator correction. Printed source summaries average statistics from resampled runs, so rounding-identical moments are not expected.
- Complete posterior rows retain their normalized original weights. No independent Gaussian parameter draws, posterior resampling, or model evaluation is performed. NGC 6397's standard fit is not a substitute for its unavailable alternative chain.
- Shared distance, calibration, tracer selection and reused observations can create dependence across records. `covariance_group` marks same-bin shared stars; it does not imply different bins have zero covariance.

The observations constrain projected tracer populations. A physically consistent host reconstruction must still provide mass-resolved local density and velocity distributions before the Step 1 drift/diffusion terms or Step 2 species exposure can be evaluated. See the [assessment](../../docs/cluster-data-assessment.md) for each cluster's role and unresolved assumptions.
