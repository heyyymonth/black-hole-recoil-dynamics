# Suitability of the six cluster datasets

The comparison adopts six distinct published standard-fit host distributions as conditional dynamical backgrounds. Joint parameter records preserve each cluster's inferred structure and stellar population. Numerical acceptance and source-specific applicability are recorded separately. The original source audit remains pinned to its release; the following decisions govern the cross-cluster calculation.

## Source integrity and common inputs

All six have a complete 13-parameter joint nested-sampling product and a candidate historical observation bundle. The 51 source files pass their recorded SHA-256 checks. Each reconstruction preserves an actual joint row, including its distance, stellar-population parameters and covariances. Individual marginal medians are not assembled into a synthetic host.

The released derived-model products contain 4,000 realizations per cluster. Their black-hole quantities are summarized separately from direct weighted chain parameters. The compatible historical GCfit code creates these products from equal-weight chains, but the exact generating seed and mapping to original nested rows are absent. Consequently, derived rows are never joined to raw posterior rows by row number. The [inventory](../results/cluster-comparison/posterior-inventory.json) records the source hashes, units and zero-population handling.

## Scientific suitability

| Cluster | Use in this comparison | Limitation requiring separate interpretation |
|---|---|---|
| 47 Tuc | Conditional multimass reference and posterior sensitivity | Spherical nonrotating model; published mass-segregation/dynamical-age consistency concern. Numerical success does not validate its geometry or encounter closure. |
| NGC 1851 | Conditional structural and scattering comparison | Four suspect compiled kinematic entries remain audit-only. They are not repaired or inserted into the historical fit. |
| NGC 3201 | Conditional sparse-population reference and joint-parameter sensitivity | The source authors flag this fit's mass-segregation and dynamical-age regime; alternative methods infer larger BH populations. The adopted standard fit defines the experiment's background rather than a unique measured BH population. |
| Omega Centauri | Massive-reservoir comparison within the stated model | Spherical nonrotating treatment omits rotation and complex population structure. Competing central-mass descriptions are not tested here. |
| M4 | Low-count comparison | Sparse original HST radial coverage and model dependence remain. One of eight pilot reconstructions fails to converge to finite extent and remains unresolved. The catalog-distance orbit passes numerical checks but differs from the catalog mean pericenter by 20.6%. Zero-BH realizations, where present, are retained. |
| NGC 6397 | Adopted standard-fit reference and sensitivity | Core-collapsed system. A distinct no-BH alternative has a summary and figures but no corresponding chain in the inspected pinned release. It remains unevaluated; the standard-fit results are conditional on their own host prescription. |

The source study identifies NGC 3201 as its strongest example of a high-segregation, dynamically young regime requiring additional model assessment; 47 Tuc and Omega Centauri also occupy that regime, although their inferred BH contents agree better with other work. Its discussion of NGC 6397 emphasizes the role of white dwarfs and little or no BH content. These source-specific statements condition the astrophysical interpretation of the adopted backgrounds. [Dickson et al., sections 3 and 5.3](https://arxiv.org/html/2308.13037)

The inspected [NGC 6397 no-BH directory](https://github.com/nmdickson/GCfit-results/tree/c45a03b1681570b6ee96be76f362894c6628e1f1/clusters/NGC6397/no_BHs) contains six figures and a summary, with no sampler chain. This is an availability statement, not a ranking of the two host prescriptions. The inventory records the adopted standard-fit chain as available for every cluster and records this inspected alternative only for NGC 6397; it makes no claim about alternative-chain availability for the other clusters.

## Quality gates and missing information

1. A source checksum, missing required joint parameter, invalid weight, or nonfinite model input blocks that calculation.
2. Independent mass, force, potential, virial, distribution-function moment and luminous-projection checks determine numerical eligibility for scattering calculations. Failed rows are retained with their failure status; no replacement row is drawn to conceal them.
3. Local coefficient refinement, symmetry, positivity and trace checks determine numerical acceptance of the evaluated states. A failed quadrature is an implementation-resolution issue, not evidence that a real cluster is unsuitable.
4. Physical model suitability is reported separately from numerical acceptance. These snapshot models do not establish a persistent bath, integer BH count distribution, binary population, strong-encounter law or a merger-conditioned launch.
5. A zero-BH realization is a valid model outcome. Its BH diffusion contribution is zero; its BH half-mass radius and probe-to-BH-reservoir mass ratio are undefined.

The 20 flagged rows in the external compilation remain outside reconstruction inputs. Gaia field counts are not membership counts, fitted profile points are correlated, and the different HST components are not independent experiments. Full phase-space covariance and Galactic-potential uncertainty are unavailable for the environmental comparison; distance variants are sensitivity cases, not joint posterior intervals.

Original observations and source definitions: [data assessment](cluster-data-assessment.md), [source manifest](../data/observed-clusters/sources.json), [pinned model release](https://github.com/nmdickson/GCfit-results/tree/c45a03b1681570b6ee96be76f362894c6628e1f1).
