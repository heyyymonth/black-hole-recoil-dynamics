# E004: companion population and frozen-potential validation

**Researcher:** Venkata Hemanth Athota. Results recorded 20 October 2026.

**Result:** this external cluster snapshot does not justify treating E003's imposed companion bath as calibrated. Only two binary black holes are present. A Plummer model matched to the current host mass and half-mass radius also fails the declared 10% ballistic-return screen at every tested kick. Neither comparison validates live remnant trajectories.

## Question and provenance

E003 supplied the density, velocities, masses, and separations of its targets by assumption. E004 asks which of those assumptions an evolved cluster actually supports, and how much the assumed potential changes an otherwise identical recoil calculation. The [protocol](e004-protocol.md) was committed before calculation as `50cf824`; its dated addendum identifies two subsequent validation checks. The final [source revision](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/7f7bc0eea4569ff1c4f5e2db443eb29daa6c7374) and all input/output hashes appear in the [run manifest](../results/run_e004.json).

We use `initial.snap0382.dat.gz` for model `N8-RV1-RG8-Z0.1`, distributed by [Nicholas Rui's cmctoolkit](https://github.com/NicholasRui/cmctoolkit/tree/a4876e7706f7b23a14c63a5be7df63a07f2dccd9/N8-RV1-RG8-Z0.1). The snapshot, conversion file, and authors' separately tabulated `PARAM0382.txt` are pinned to revision `a4876e7706f7b23a14c63a5be7df63a07f2dccd9` with SHA-256 hashes in the [configuration](../experiments/e004.json). Selection used the available full evolved example, not a search for a large discrepancy. This is one simulated cluster state, not an observed cluster or a representative population.

The snapshot age is 12.7323 Gyr, current mass is 224,633 solar masses, and three-dimensional half-mass radius is 4.79771 pc. A matched Plummer model has scale radius 3.67706 pc. These differ from E003's illustrative host: E004 tests transfer of the assumptions to this specific external state, not E003's internal arithmetic.

## Companion inventory

| Population | Systems | Median radius / host half-mass radius |
|---|---:|---:|
| All systems | 682,583 | 1.2468 |
| All binaries | 31,122 | 1.0991 |
| Single black holes | 21 | 0.04018 |
| Binaries with at least one black hole | 6 | 0.03851 |
| Binary black holes (both components) | 2 | 0.01014 |

The final two categories overlap; they must not be added. There are 29 black-hole objects in total. All binaries and all black holes are inappropriate substitutes for the specific BH–BH targets assumed in E003.

The two BH–BH system masses are 26.672 and 29.821 solar masses; their semimajor axes span 0.23958–12.40947 AU. Their median radius is 0.04864 pc. E003's fixed 60-solar-mass, 0.1-AU targets and unsegregated distribution are therefore not calibrated by this snapshot. Two systems fail the predeclared minimum-count screen of 30; that threshold is an operational rule, not a theorem that 30 would suffice. No continuous distribution is fitted, and repeated synthetic draws would not add independent information.

The [population table](../results/e004/population.csv) also reports descriptive radial and speed-fraction CDF distances. Host-mass comparisons use exact step CDFs on their combined jump locations. Plummer comparisons use the continuous analytical distributions. These distances have no independent-sample significance interpretation. Speeds are scaled by the *matched Plummer* local escape speed: discrepancies mix potential and tracer-population differences and cannot alone establish anisotropy or equipartition failure. Speeds above this reference escape are not evidence of escape from the simulated cluster.

## Ballistic sensitivity to potential shape

For spherical shells, with potential zero at infinity,

```math
\Phi(r)=-G\left[\frac{M(\lt r)}{r}+\sum_{r_i\geq r}\frac{m_i}{r_i}\right].
```

Exterior shells contribute to binding energy even though they exert no interior force. The snapshot central escape speed is 26.547 km/s, compared with 22.924 km/s for the matched Plummer model. The latter matches total mass and half-mass radius but is shallower centrally.

Both models receive the **same physical kick**. The first return is a center-to-center ballistic excursion; it is not dissipative settling. Values below use the declared 2048-shell calculation.

| Kick / Plummer escape | Kick (km/s) | Snapshot return (Myr) | Plummer return (Myr) | Plummer / snapshot |
|---|---:|---:|---:|---:|
| 0.3 | 6.8771 | 0.26250 | 0.77588 | 2.9557 |
| 0.8 | 18.3388 | 0.89636 | 2.52230 | 2.8139 |
| 0.9 | 20.6312 | 1.39333 | 6.20464 | 4.4531 |

All three differences exceed the declared 10% screen. Apocenters and all four shell resolutions are retained in [trajectories.csv](../results/e004/trajectories.csv). The agreement of compressed and full-shell results below leaves this conclusion unchanged.

![Population, potential, and ballistic-return comparison](../results/e004/snapshot_comparison.png)

## Verification and what it establishes

Eight derived quantities agree with the authors' tabulated summary: age, total mass, black-hole mass, system count, object count, binary count, black-hole count, and BH–BH binary count. Counts match exactly; masses use relative tolerance 1e-7 to accommodate rounded binary center-of-mass entries; age uses 1e-12. This checks parsing against a separately processed version of the same simulation, not a second physical dataset. Snapshot masses are already in solar masses. Velocity uses the N-body time conversion, whereas age uses the distinct CMC relaxation-time conversion. Conversion scripts are parsed, never executed. See the [CMC output conventions](https://clustermontecarlo.github.io/CMC-COSMIC/output/index.html).

Compression uses mass-weighted harmonic shell radii, preserving total mass and central potential. The final two resolutions differ by less than 0.005% in return time and 0.023% in apocenter. This is not an absolute error bound: comparison with **all 680,110 unique shell radii** gives return errors of 0.1555%, 0.00131%, and 0.000568%. All satisfy the 1% screen. The finite-shell central cavity is force-free, so arbitrarily small kicks cannot be interpreted using the smooth-core harmonic limit.

The full suite has 108 passing tests, including direct potential sums, force derivatives, analytical one-shell returns, independent adaptive quadrature, smooth-Plummer recovery, unit conversion, binary accounting, and discrete-CDF examples. Recomputing stored tables is a separate reproducibility check. [Verification output](../results/e004/verification.txt) records the commands and results.

**Live trajectories remain unvalidated.** CMC randomly samples orbital phases during Monte Carlo evolution, so joining successive snapshot radii would not reconstruct a remnant path. The frozen spherical calculation omits relaxation, friction, strong encounters, tides, and evolving host structure. Its isolated-potential escape criterion is distinct from CMC's tidal-removal treatment. These distinctions follow the [CMC method, sections II.5–II.7](https://arxiv.org/html/2106.02643).

## Reproduction and output conventions

After installing the root README dependencies:

```sh
python experiments/fetch_e004.py
python experiments/verify_e004.py
```

The first fetch requires network access and approximately 46 MB of local storage. Subsequent checks reuse checksum-verified files in ignored `tmp/cmc/`. To regenerate, first commit changed source inputs, then run `python experiments/run_e004.py`. The runner requires committed inputs and records their hashes. Raw upstream files are neither redistributed nor relicensed as MIT.

- `population.csv`: five overlapping cohort summaries; lengths in pc except binary semimajor axes in AU, masses in solar masses, counts integral, CDF distances and speed fractions dimensionless. Empty binary columns mean not applicable.
- `trajectories.csv`: twelve numerical cases (four resolutions, three kicks), not twelve independent clusters. Kicks are km/s, lengths pc, return times Myr. `eta_snapshot` uses snapshot escape; `eta_plummer` defines the identical physical kick. Relative return error is Plummer/snapshot minus one. Central potential is in units of GM/b. Bound flags refer only to the isolated potential.
- `profiles.csv`: 160 plotting-grid locations; radius in units of host half-mass radius, potentials in GM/b, cumulative fractions dimensionless. Grid points are not independent observations.
- `summary.json`: scales, accounting checks, count screen, adjacent-resolution and full-shell checks, and explicit calibration/trajectory status. `velocity_unit_kms` is the matched-model scale sqrt(GM/b), not the file's N-body conversion factor.

The [literature context](cluster-dynamics-context.md) specifies the two missing comparisons: calibrating a joint companion population across cluster histories, and controlled live-kick experiments measuring return and encounter accessibility. E004 diagnoses those needs; it does not yet answer how often remnants acquire companions or merge again.
