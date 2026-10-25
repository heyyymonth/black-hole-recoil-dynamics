# Common cluster host schema

**Scope:** preparation and comparison of public observations and fitted-model products for NGC 0104 (47 Tuc), NGC 1851, NGC 3201, NGC 5139 (Omega Centauri), NGC 6121 (M4), and NGC 6397. A common input structure is proposed below. It does not establish a compatible observation/posterior/software pairing, approve a cluster model, or assign recoil experiment parameters. No host reconstruction or remnant dynamics is performed here.

## One structure, six conditional hosts

Each cluster has its own observations, uncertainties, source selection and conditional host distribution. The common structure standardizes names, units and meanings; it preserves physical differences between clusters. It neither pools all observed stars into one cluster nor treats the six deliberately selected clusters as a representative population.

For cluster c, the intended statistical object is

```math
p(\Theta_c\mid\mathcal D_c,\mathcal M_c),
```

where the data collection, model family and parameter vector remain separately identified. A shared model family is a hypothesis to assess against each cluster, not a consequence of storing all six in the same table. Priors, posterior samples and their original weights must remain associated with the model and observations that produced them.

| Layer | Common fields | Meaning and use |
|---|---|---|
| Source record | Cluster identifier and aliases; release; immutable URL/revision; file hash; citation; reuse terms | Identifies the exact bytes and their provenance |
| Observation series | Observable; projected radius and bin bounds; native unit; value; asymmetric errors; tracer selection; instrument/source; footprint; covariance or covariance status; quality flags | Holds measured or published observational summaries without changing their physical interpretation |
| Fitted host | Model family and implementation revision; data selection; parameter ordering; fixed values; priors; joint samples and weights; convergence diagnostics | Describes an inference conditional on a declared data/model pairing |
| Reconstructed background | Species masses and number distributions, potential, streaming and velocity structure, boundary and reference frame | Requires a compatible model implementation and validation; is not supplied by a scalar catalogue row |
| Recoil calculation | Merger-conditioned launch, remnant mass, scattering prescription, termination and observation interval | Remains a separate, unassigned configuration |

Retain proper motions in mas/yr and projected radii in angular units in the observation layer. Distance-dependent conversions belong in a model-aware view: transverse velocity is 4.74047 D times proper motion when D is in kpc and velocity is in km/s. The same uncertain distance affects all converted measurements for a cluster. Line-of-sight dispersions already expressed in km/s must not receive this conversion. Radial proper motion, tangential proper motion and their combined dispersion from the same stars are correlated views, not three independent datasets.

Projected radius R is not three-dimensional radius r. Surface brightness and number density are not spatial mass density. A published mean mass, half-mass radius or central escape speed does not uniquely determine a potential or a friction coefficient. Observations of visible tracers do not directly measure the black-hole velocity distribution. The [GCfit model description](https://arxiv.org/html/2308.13037v2) explains how multiple observables enter a conditional multimass inference.

## Common posterior parameter pattern

The older [released GCfit results](https://github.com/nmdickson/GCfit-results/tree/c45a03b1681570b6ee96be76f362894c6628e1f1) expose the same thirteen sampled parameter names for the six standard fits. Their conventions are documented by the [historical model implementation](https://github.com/nmdickson/GCfit/blob/e4154e8e3e2505d2c49248a71c1eebcea5844442/gcfit/core/data.py). This historical implementation is a convention reference; its exact compatibility with each released chain remains to be established.

| Released name | Canonical interpretation | Required handling |
|---|---|---|
| W0 | Dimensionless central-potential parameter of the specified distribution-function model | Preserve the model's potential and mass-scale convention; it is not a measured central escape speed |
| M | Total model mass in units of one million solar masses | Multiply by one million for a solar-mass view; the direct Model constructor uses solar masses |
| rh | Three-dimensional half-mass radius in pc | Keep distinct from projected half-light, core and truncation radii |
| ra | Base-10 logarithm of the model's dimensionless anisotropy radius | Exponentiate only when passing to the matching model; do not label the stored value as pc |
| g | Dimensionless truncation-shape parameter | Preserve the selected distribution-function family |
| delta | Dimensionless dependence of the velocity scale on stellar mass | Preserve the mass-component prescription; do not equate it directly with a measured dispersion |
| a1, a2, a3 | Slopes of the model's piecewise stellar initial mass function | Retain mass breaks, stellar evolution and remnant prescriptions from the matching implementation |
| BHret | Black-hole retention parameter in percent under the historical evolution prescription | Divide by 100 only where that prescription requires a fraction; it is not the present-day fraction of cluster mass in black holes |
| d | Distance in kpc | Propagate jointly with structural parameters and angular observations |
| s2 | Additional number-density uncertainty nuisance parameter | Retain its original units and likelihood role; do not use it as a dynamical input |
| F | Multiplicative stellar-mass-function uncertainty nuisance parameter | Retain its likelihood role; do not interpret it as a mass or velocity factor |

The posterior is a distribution over complete parameter rows. It is not thirteen independent normal distributions assembled from means and errors. Original parameter order, fixed values, priors and sampling weights must accompany every row. Posterior moments are diagnostics, not substitutes for the joint distribution. In particular, the released summaries distinguish posterior standard deviations from the numerical uncertainty in their estimates.

For nested-sampling output, normalize posterior weights from the stored log weights, retaining sample alignment. The [historical analysis code](https://github.com/nmdickson/GCfit/blob/e4154e8e3e2505d2c49248a71c1eebcea5844442/gcfit/analysis/runs.py) uses log weights and final log evidence when computing parameter moments. Raw nested-sampling rows are not equally weighted. Numerical inspection of these outputs does not require resampling or solving the cluster model.

The nuisance coordinates s2 and F remain in the archived posterior and likelihood record. A physical-host view marginalizes over them through the original weighted joint samples; it does not reset them and refit implicitly. Likewise, an excluded likelihood remains excluded when reproducing a published fit. The standard NGC 0104 summary excludes pulsar components. A posterior must not be multiplied by the likelihood of its own source observations as if those were independent new evidence.

## Common scales without a universal cluster

For a compatible individual host realization, use its own total mass and three-dimensional half-mass radius to define

```math
r_0=r_h,\qquad
v_0=\sqrt{\frac{GM}{r_h}},\qquad
t_0=\frac{r_h}{v_0}.
```

Store physical r_h in pc, M in solar masses and v_0 in km/s. Convert t_0 from pc/(km/s) to the declared physical time unit, such as Myr, explicitly. The scale v_0 is a convenient gravitational scale; it is neither the measured central velocity dispersion nor the escape speed. No numerical scale values are assigned here.

Dimensionless coordinates and fields are

```math
\widetilde{\mathbf r}=\frac{\mathbf r}{r_h},\qquad
\widetilde{\mathbf v}=\frac{\mathbf v}{v_0},\qquad
\widetilde t=\frac{t}{t_0},\qquad
\widetilde\Phi=\frac{\Phi}{v_0^2},\qquad
\widetilde\rho_j=\frac{\rho_j r_h^3}{M}.
```

A potential zero and outer boundary must still be specified. These scalings allow comparison at corresponding dimensionless radii and times, while retaining each host's concentration, anisotropy, truncation, mass spectrum and component fractions. They do not make different clusters homologous or justify using a common Plummer profile. The [lowered-isothermal family](https://arxiv.org/abs/1508.02120) supplies examples of structural freedom that mass and radius alone do not fix.

Scales must be computed from each complete posterior draw when propagating host uncertainty. Combining a mass from one catalogue, a radius from another fit and a distance from an unrelated posterior would create a parameter combination without the original joint support. A plot of projected R/r_h remains projected and cannot be used directly as a spatial density profile.

## Role of each candidate

| Cluster | Proposed assessment role | Condition before a shared host family is defensible |
|---|---|---|
| NGC 3201 | Initial baseline candidate with public photometric, kinematic and mass-function information | Reproduce the source fit and check tracer selection, anisotropy and the radial domain affected by tides |
| NGC 0104 / 47 Tuc | Initial baseline candidate with extensive public observations | Assess rotation and its treatment in the dispersion data; a spherical nonrotating model is an approximation to test |
| NGC 1851 | Additional concentrated-host comparison | Reconcile flagged compilation rows against original series and establish the correct posterior/data pairing |
| NGC 6121 / M4 | Additional nearby-host comparison | Check the competing interpretations of central dark mass and the tracer/mass-function selections |
| NGC 6397 | Core-collapsed comparison with greater model sensitivity | Keep the standard and zero-black-hole alternatives distinct; the released no-BH directory contains summaries and figures but no sampler chain in the inspected revision |
| NGC 5139 / Omega Centauri | Rotation, flattening and tracer-selection comparison | Assess an axisymmetric or rotating representation where required; reconcile changing LOS/PM releases and tracer coverage before combining them |

These roles describe a data assessment, not an approved six-cluster experiment. Individual-cluster limitations and the treatment of core-collapsed systems are discussed in [Dickson et al.](https://arxiv.org/html/2308.13037v2). The public results summaries were updated in January 2024. Omega Centauri's GCfit LOS resource changed [later that January](https://github.com/nmdickson/GCfit/commit/30dac602598aee6af1d0576615c4365833db2f5b); its proper-motion source changed again in the [April 2026 resource update](https://github.com/nmdickson/GCfit/commit/679d0cdc0e39af848107edd8679035bb7e52dff5). The latest software resource is therefore not automatically the observation file used by an older posterior. Replacing data requires an explicit fit or justified inference update, not a file substitution.

## Connection to the two steps

[Step 1](e001-ballistic-return.md) requires the species-resolved number distribution F_j(r,u,t), not only total density or a single central dispersion. The same distribution supplies mass density, streaming and velocity structure for smooth gravity and the selected drift/diffusion closure. A fitted equilibrium host can be a candidate source for these functions after reconstruction and observational checks. Its inferred dark components remain conditional on evolutionary and structural assumptions. Rotation, anisotropy and non-Gaussian velocity distributions cannot be repaired merely by changing units.

[Step 2](e002-core-exposure.md) must use the same host realization, physical epoch and reference frame as the motion that produced a path. Its region boundary and species densities then define residence and exposure. Density curves from a different fitted host cannot be attached to that path without an explicit comparison design. Host uncertainty should remain associated with all quantities passed between the two steps.

The following remain unassigned in the [calculation handoff](calculation-handoff.md) and [input register](../formulation/inputs.json):

- Merger epoch and any host evolution over the observation interval. Present-day observations do not identify a historical merger environment.
- Merger-conditioned position and progenitor center-of-mass velocity, recoil magnitude/direction, and remnant mass. A visible-star distribution is not automatically the merger launch distribution.
- Scattering closure, impact-parameter cutoffs, Coulomb logarithms and the separation of weak scattering from resolved strong encounters.
- External potential and host orbit, coordinate and potential-zero conventions, escape boundary, strong-encounter treatment and model-validity cutoff.
- Observation horizon, finite-region boundary, exposure normalization and termination/censoring rules.
- Numerical representation, time step and convergence criteria.

The data-analysis outcome needed before those choices is a traceable per-cluster observation inventory, quality and overlap flags, a verified joint posterior description where available, and a documented account of which source pairings remain unresolved. Common fields make these objects comparable; observational adequacy and dynamical validity remain separate questions.
