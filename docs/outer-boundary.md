# Outer validity and tidal boundaries

The boundary question is part of Step 1: how far can the reconstructed environment support a remnant calculation, and what would count as loss of validity or escape? Step 2 must preserve that distinction when accumulating residence and exposure. A stopped calculation is not evidence that the remnant has escaped or can no longer return.

## Physical basis

Three radii have different meanings. A fitted limiting radius describes a density model, a Jacobi scale describes a tidal opening under specified dynamics, and an accuracy boundary limits a particular approximation. The classical empirical profile of [King (1962)](https://articles.adsabs.harvard.edu/pdf/1962AJ.....67..471K) motivates the first; the tidal-tensor treatment of [Renaud, Gieles and Boily (2011)](https://doi.org/10.1111/j.1365-2966.2011.19531.x) defines the second for appropriate stationary configurations. [Fukushige and Heggie (2000)](https://doi.org/10.1046/j.1365-8711.2000.03811.x) demonstrate why opening an escape channel does not imply immediate departure. These are foundational or directly relevant studies, without an asserted citation ranking.

For NGC 3201, [Wan et al. (2021)](https://doi.org/10.1093/mnras/stab306) distinguish a 36.1 pc King limit from a cited 83.46 pc Jacobi estimate. The captured structural catalog lists 85.86 pc, with its header specifying Eq. 8 of [Webb et al. (2013)](https://doi.org/10.1088/0004-637X/764/2/124). That equation corrects a pericentric radius using eccentricity and orbital phase. We retain the catalog value as a comparison, not an instantaneous tidal tensor or an independent measurement to average with the other estimates.

## A specified external field

The fixed reconstructed host retains its original mass, potential and species. For a controlled boundary diagnostic, adopt a spherical logarithmic Galactic potential with constant circular speed:

```math
\Phi_G(\mathbf R)=v_c^2\ln(|\mathbf R|/R_0),\qquad
\mathbf a_{\rm ext}=-\nabla\Phi_G(\mathbf R+\mathbf r)+\nabla\Phi_G(\mathbf R).
```

The baseline uses catalog Galactocentric distance 8.93 kpc and a 220 km/s literature benchmark. The external distance is a separately fixed environmental setting, not a joint draw with the host's inferred heliocentric distance. Varying the circular speed by ±10% and using a separate 13 kpc circular surrogate tests assumptions; these choices do not define confidence intervals. The latter follows a comparison discussed by Wan et al., and is not the cluster's present position or a reconstructed orbit. No eccentric Galactic orbit is integrated here.

In translating nonrotating axes, the linear external tensor has eigenvalues proportional to (1, −1, −1). For a circular orbit only, adding centrifugal acceleration gives effective eigenvalues proportional to (2, 0, −1). The common scale is v_c²/R². Coriolis acceleration acts on velocity in the rotating frame and is absent from the nonrotating transport equation. Mixing these tensors would insert an incorrect force into Step 1.

The circular Hill approximation gives an extended-host saddle radius through

```math
\frac{G M(<r_J)}{r_J^3}=\frac{2v_c^2}{R^2}.
```

This uses enclosed mass rather than treating all material outside the saddle as interior mass. Exact circular inward and outward radial force balances provide a separate check of the local tidal approximation. The host is not trimmed or re-equilibrated during this diagnostic.

## A proposed accuracy limit

For this spherical Galactic field and any direction at radius 0 < r < R, the maximum differential acceleration occurs toward the Galactic center. Its exact value is v_c² r/[R(R−r)]. The maximum relative omitted force is therefore

```math
\eta(r)=\frac{v_c^2 r^3}{G M(<r) R(R-r)},\qquad
r_\epsilon=\inf\{r>0:\eta(r)>\epsilon\}.
```

We propose ε = 0.01 for an **isolated-host control calculation**, with 0.03 and 0.1 as explicit tolerance comparisons. This is a chosen local force budget, not a universal physical law or a bound on accumulated trajectory error. For an interior launch, the first outward crossing would be a validity stop, not an escape event. An initial point already outside is invalid at time zero; a boundary launch needs an explicit side convention. A calculation including the external force has a different validity domain and must not inherit this control's stop automatically.

## Computed findings

For the baseline:

| Quantity | Result | Meaning |
|---|---:|---|
| Extended-host circular saddle radius | 84.13 pc | Conditional tidal opening |
| Point-mass approximation | 85.77 pc | Uses all reconstructed mass |
| Exact circular inward / outward saddles | 84.00 / 84.27 pc | Checks the Hill approximation |
| 1% / 3% / 10% force limits | 20.67 / 31.12 / 47.94 pc | Approximation tolerances |
| Omitted force at r_h / 3r_h | 0.175% / 2.962% | Connects to existing local coefficient radii |
| Reconstructed mass outside the circular saddle | 5.63% | Host-consistency diagnostic, not escaped mass |
| Central Hill opening speed | 16.68 km/s | Circular conservative energy diagnostic, not a recoil-retention threshold |

The formal model truncation is 981 pc. Its disagreement with the tidal scale cannot be resolved merely by deleting exterior stars: changing the distribution function and host equilibrium requires a separate calculation. The numerical proximity to a published 83.46 pc estimate is a scale check under related assumptions, not a validation of the whole host or of its outer velocity distribution.

Across the ±10% speed cases, the saddle ranges from 78.79 to 90.44 pc and the 1% force limit from 19.22 to 22.39 pc. At the separate 13 kpc surrogate the corresponding values are 108.79 and 27.40 pc. These sensitivities are conditional on the same single host realization.

The calculation preserves the numerical meaning of the earlier local kernels. It restricts their interpretation: the state at 3r_h lies beyond the proposed 1% isolated-force domain, and the sparse central black-hole component still needs discrete-encounter validation. Neither issue is fixed by convergence of an integral. Extending motion to the outskirts needs a Galactic trajectory, time-dependent external forces and an appropriate host description, followed by explicit tracking of departures and possible returns. The escape-opening energy is conserved only in the stationary conservative circular problem; stochastic scattering and noncircular tides invalidate that conservation.

## Numerical evidence

Root residuals are below 9.5e-13. Independent integration of density agrees with enclosed mass at tested boundary radii to 1.2e-7, and thinning the host radial mesh changes radii by less than 1.5e-7 in relative terms. Tests also cover analytic point-mass and extended Plummer saddles, exact angular force bounds, the correct frame tensor, central mass scaling and force–potential consistency. These are numerical checks, not observational uncertainties.

[Protocol and inputs](../analysis/boundary/protocol.json) · [Source register](../analysis/boundary/sources.json) · [Reproduction](../analysis/boundary/README.md) · [Scenario results](../results/outer-boundary/scenarios.csv) · [Numerical summary](../results/outer-boundary/summary.json)
