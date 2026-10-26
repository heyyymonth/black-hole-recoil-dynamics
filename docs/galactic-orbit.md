# Galactic orbit and time-dependent host validity

This calculation supplies NGC 3201's external reference motion for Step 1 and tests the spatial and temporal assumptions that any subsequent residence calculation must inherit. It integrates the **cluster around the Galaxy**, not a merger remnant inside the cluster. The stellar distribution is held fixed for a diagnostic; no host refit or tidal re-equilibration is claimed.

## Sources and defined cases

The [source record](../data/orbit/ngc3201-source.json) and original [row extract](../data/orbit/ngc3201-row.html) preserve the database's `/newdata/orbits.html` values. The published Cartesian point is (7.59, −4.65, 0.71) kpc and (−250.21, −200.78, 149.10) km/s. X points toward the Sun, Y in the direction of solar motion, and Z north. We use those coordinates directly, avoiding an invented solar-frame reconstruction.

The native distance is 4.74 ± 0.04 kpc. Proper motions are 8.350 and −1.976 mas/yr, both with quoted 0.005 mas/yr errors and correlation 0.18. The radial velocity is 495.38 ± 0.06 km/s. The table is a compiled source associated with [Vasiliev and Baumgardt (2021)](https://doi.org/10.1093/mnras/stab1475) and [Baumgardt and Vasiliev (2021)](https://doi.org/10.1093/mnras/stab1474). Quoted errors are retained; a full covariance is unavailable.

We use [Irrgang et al. (2013), Model I](https://doi.org/10.1051/0004-6361/201220540), a fixed axisymmetric bulge, disk and spherical halo potential. This reference model is chosen for comparison with the catalog, not asserted to be the latest or uniquely correct Galactic potential. Its exact physical normalization, parameters and source version are specified in the [protocol](../analysis/orbit/protocol.json).

Four cases use the catalog distance, the selected host's 4.54335 kpc distance, and the catalog distance shifted by its quoted ±0.04 kpc error. A distance change adds ΔD times the line-of-sight unit vector to position and 4.74047 ΔD times the proper-motion tangent vector to velocity. The catalog solar-frame anchor is preserved. These are conditional sensitivities, not a combined posterior or independent likelihood multiplication.

The reference orbit obeys

```math
\ddot{\mathbf R}=-\nabla\Phi_G(\mathbf R),\qquad
\mathbf a_{\rm ext}(\mathbf r,t)=\mathbf a_G(\mathbf R(t)+\mathbf r)-\mathbf a_G(\mathbf R(t)).
```

Integration spans 1 Gyr forward in the specified potential. This duration samples repeated orbital phases; it does not assume that the real cluster's internal structure stays unchanged for a Gyr.

## A sufficient force-accuracy domain

Let T = −Hessian(Φ_G), in fixed nonrotating axes. If B(R,r) bounds the operator norm of T everywhere in the radius-r ball around R, then the fundamental theorem of calculus gives

```math
|\mathbf a_{\rm ext}|\le r B(\mathbf R,r),\qquad
\eta_B=\frac{B(\mathbf R,r)r^3}{G M(<r)}.
```

We solve η_B = 0.01 at each orbital phase. The bound is analytic, not an angular sampling estimate. The sum of component bounds applies to the full potential. For the disk Hessian, its diagonal positive part has minimum eigenvalue at least 1 in its dimensionless matrix, while its subtracted rank-one term has norm at most 3. This yields an operator-norm bound proportional to max(2, Jzz)/d³. Minimum cylindrical radius and absolute height in the spatial ball give the conservative supremum. Bulge and halo bounds follow their spherical radial/tangential eigenvalues. The [implementation](../analysis/orbit/galaxy.py) states the expressions and units.

This bounds the instantaneous error from omitting the Galactic force, conditional on the supplied host. It does not establish small cumulative orbit error or a valid stellar distribution function.

## Frozen-frame scale and response times

For comparison only, define the minimal radial-frame angular velocity and effective tensor:

```math
\boldsymbol\omega=\frac{\mathbf R\times\dot{\mathbf R}}{|\mathbf R|^2},\qquad
\mathsf T_f=\mathsf T+|\boldsymbol\omega|^2\mathsf I-\boldsymbol\omega\boldsymbol\omega^{\mathsf T}.
```

The diagnostic radius r_f satisfies G M(<r_f)/r_f³ = λ_max(T_f). On the eccentric inclined orbit this is **not a stationary Jacobi radius or an escape boundary**. The omitted Euler term is −dot(ω) × r; its coefficient can be comparable to the extensive effective eigenvalue. Coriolis acceleration is also needed if dynamics is actually evaluated in that rotating frame. The nonrotating transport equation continues to use the exact differential force above.

We also calculate t_dyn ||dot(T)||₂ / ||T||₂, where t_dyn = sqrt(r³/[G M(<r)]). It compares host dynamical time with changes in forcing, including orientation in the nonrotating axes. A separate indicator uses variation of the tensor norm alone. These are timescale diagnostics, not heating rates, relaxation predictions or calibrated adiabatic corrections.

## Results and consequences

The catalog-reference orbit has mean pericenter 8.333 kpc and apocenter 25.362 kpc, compared with catalog estimates 8.33 and 25.82 kpc. The apocenter differs by 1.77%; source rounding, averaging and exact model configuration are not reconstructed. This is a scale check rather than exact reproduction or a new observational fit.

The selected-host-distance case gives:

| Diagnostic | Result over the defined orbit |
|---|---:|
| Mean pericenter / apocenter | 8.215 / 23.810 kpc |
| Sufficient 1% isolated-force radius | approximately 8.2–39.8 pc |
| Refined minimum sufficient radius | 8.173 pc |
| Proposed fixed isolated-control radius | 8 pc |
| Frozen-frame radius, only a diagnostic | 62.1–183.9 pc |
| Maximum omitted-force bound at r_h / 3r_h | 1.73% / 29.5% |
| Maximum forcing/dynamical-timescale ratio at r_h | 1.16 (1.14 for norm variation alone) |
| Maximum Euler-coefficient / extensive-eigenvalue ratio | 0.67 |
| Fixed host mass beyond r_f | 1.37–8.74%, not an escaped fraction |

Disk crossings tighten the force domain and make the forcing vary rapidly. Thus the earlier 20.67 pc static-logarithmic limit is not transferable unchanged to the orbit. The proposed 8 pc value rounds the refined minimum downward for this conditional isolated-control experiment. An orbit with exact Galactic forcing is not limited by this omitted-force criterion, but still needs a valid host and scattering treatment.

A remnant launched outside a chosen control domain is invalid at time zero. An interior launch is censored at the first validity exit. Neither event means irreversible escape, and Step 2 cannot infer complete future exposure from a stopped trajectory.

The host is not globally validated by the orbit. [Bertin and Varri (2008)](https://doi.org/10.1086/592684) construct self-consistent nonspherical models for a circular orbit in a frozen tidal field. That result does not authorize replacing a changing eccentric-orbit tide with one energy cutoff. Our calculation identifies the supplied host's required evolution problem: consistent species distributions and gravity under the time-dependent tide, with any escape or re-entry treatment explicit. It does not solve that problem by deleting exterior stars. The sparse central black-hole population remains a separate unresolved collision-closure issue.

## Validation and remaining uncertainty

Force–potential and tensor derivatives, disk symmetry, halo matching, all-direction force bounds and coordinate transformations pass independent checks. Energy and axial angular momentum are conserved to better than the registered 1e-8 tolerances; refined integration and time reversal pass. Halving output spacing changes the minimum sufficient radius by 0.0674%, and the peak timescale ratio by less than 4.1e-8. These numerical errors are distinct from missing Galactic-model, solar-frame, host-posterior and measurement covariance uncertainty.

[Reproduce the calculation](../analysis/orbit/README.md) · [Protocol](../analysis/orbit/protocol.json) · [Summary](../results/orbit/summary.json) · [Orbit and field records](../results/orbit/tracks.csv) · [Turning and disk-crossing events](../results/orbit/events.csv)
