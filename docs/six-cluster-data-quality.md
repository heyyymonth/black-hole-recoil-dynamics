# Suitability of the six cluster datasets

The comparison uses six distinct published standard-fit host distributions. No shared stellar population is fitted, no missing input is imputed, and no cluster is accepted merely because a numerical model can be constructed. The original source audit remains pinned to its release; the following decisions govern the cross-cluster calculation.

## Source integrity and common inputs

All six have a complete 13-parameter joint nested-sampling product and a candidate historical observation bundle. The 51 source files pass their recorded SHA-256 checks. Each reconstruction preserves an actual joint row, including its distance, stellar-population parameters and covariances. Individual marginal medians are not assembled into a synthetic host.

The released derived-model products contain 4,000 realizations per cluster. Their black-hole quantities are summarized separately from direct weighted chain parameters. The compatible historical GCfit code creates these products from equal-weight chains, but the exact generating seed and mapping to original nested rows are absent. Consequently, derived rows are never joined to raw posterior rows by row number. The [inventory](../results/cluster-comparison/posterior-inventory.json) records the source hashes, units and zero-population handling.

## Scientific suitability

| Cluster | Use in this comparison | Limitation requiring separate interpretation |
|---|---|---|
| 47 Tuc | Conditional multimass reference and posterior sensitivity | Spherical nonrotating model; published mass-segregation/dynamical-age consistency concern. Numerical success does not validate its geometry or encounter closure. |
| NGC 1851 | Conditional structural and scattering comparison | Four suspect compiled kinematic entries remain audit-only. They are not repaired or inserted into the historical fit. |
| NGC 3201 | Sparse-population and model-sensitivity example | The source authors specifically caution that this fit may struggle with the cluster's dynamical regime; alternative methods infer larger BH populations. It is not a clean observational anchor for a few-BH thesis. |
| Omega Centauri | Massive-reservoir comparison within the stated model | Spherical nonrotating treatment omits rotation and complex population structure. Competing central-mass descriptions are not tested here. |
| M4 | Low-count comparison | Sparse original HST radial coverage and model dependence remain. One of eight pilot reconstructions fails to converge to finite extent and remains unresolved. The catalog-distance orbit passes numerical checks but differs from the catalog mean pericenter by 20.6%. Zero-BH realizations, where present, are retained. |
| NGC 6397 | **Standard-fit sensitivity only** | Core-collapsed system. The preferred no-BH alternative has a summary and figures but no corresponding chain in the inspected pinned release. The standard-fit calculation cannot substitute for that missing alternative. |

The source study identifies NGC 3201 as its strongest example of a suspicious high-segregation, dynamically young fit; 47 Tuc and Omega Centauri also occupy the flagged regime, although their inferred BH contents agree better with other work. Its discussion of NGC 6397 emphasizes the role of white dwarfs and little or no BH content. These are limitations of physical interpretation, not reasons to erase conditional calculations. [Dickson et al., sections 3 and 5.3](https://arxiv.org/html/2308.13037)

## Quality gates and missing information

1. A source checksum, missing required joint parameter, invalid weight, or nonfinite model input blocks that calculation.
2. Independent mass, force, potential, virial, distribution-function moment and luminous-projection checks determine numerical eligibility for scattering calculations. Failed rows are retained with their failure status; no replacement row is drawn to conceal them.
3. Local coefficient refinement, symmetry, positivity and trace checks determine numerical acceptance of the evaluated states. A failed quadrature is an implementation-resolution issue, not evidence that a real cluster is unsuitable.
4. Physical model suitability is reported separately from numerical acceptance. These snapshot models do not establish a persistent bath, integer BH count distribution, binary population, strong-encounter law or a merger-conditioned launch.
5. A zero-BH realization is a valid model outcome. Its BH diffusion contribution is zero; its BH half-mass radius and probe-to-BH-reservoir mass ratio are undefined.

The 20 flagged rows in the external compilation remain outside reconstruction inputs. Gaia field counts are not membership counts, fitted profile points are correlated, and the different HST components are not independent experiments. Full phase-space covariance and Galactic-potential uncertainty are unavailable for the environmental comparison; distance variants are sensitivity cases, not joint posterior intervals.

Original observations and source definitions: [data assessment](cluster-data-assessment.md), [source manifest](../data/observed-clusters/sources.json), [pinned model release](https://github.com/nmdickson/GCfit-results/tree/c45a03b1681570b6ee96be76f362894c6628e1f1).
