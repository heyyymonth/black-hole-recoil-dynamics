# Collisionless residence and encounter budgets

The collisionless control connects the adopted cluster potentials to timed remnant paths and occupation of a finite region. It establishes the zero-collision reference for the matched encounter comparison. All physical species contribute to smooth gravity; collisional drift and random deflections are disabled for every species. Evaluating encounter coefficients along a reference path does not apply those coefficients to its motion.

## Conditional launches and region

Each of the six verified reference hosts supplies its own BH density radius `r_rho`, where BH mass density falls to half its central value. The region is the spatial sphere `r < r_rho`. Define

```math
v_{\rm dyn}=\sqrt{-r_\rho\psi'(r_\rho)},
\qquad t_{\rm dyn}=\frac{r_\rho}{v_{\rm dyn}},
\qquad H=10t_{\rm dyn}.
```

The time scale is converted from pc/(km/s) to Myr. Both launch directions start at `r0=0.5 r_rho` with speed

```math
v_0=\sqrt{2[\psi(r_0)-\psi(2r_\rho)]}.
```

One velocity points radially outward and the other is tangential. They have the same energy and different angular momenta. The radial apocenter is `2 r_rho`; the tangential turning points follow from the same potential and its nonzero angular momentum. These are diagnostic total velocities, not inferred recoil velocities or a merger-conditioned population. No cluster weighting is assigned.

The finite horizon is ten dynamical times, not ten orbital periods. Residence includes every visit to the region. A return does not bring the remnant to permanent rest. Results at this horizon retain the launch phase; they are not cycle-averaged residence fractions.

## Evolution and stopping

The Cartesian equations are integrated in dimensionless coordinates, with scales `r_rho`, `v_dyn` and `t_dyn`:

```math
\frac{d\mathbf r}{dt}=\mathbf v,
\qquad
\frac{d\mathbf v}{dt}=\psi'(r)\frac{\mathbf r}{r}.
```

The force is the derivative of the same reconstructed potential used for the launch energy. The exact center has zero force. The original potential interpolation is retained at every positive radius. Its small central derivative residual is tested through orbit refinement and energy conservation.

The control ends at the horizon or the first outward crossing of the smaller of the host support radius and that host's recorded minimum one-percent isolated-force domain. This boundary is a validity stop, not escape. The prior force diagnostic bounds an instantaneous neglected acceleration; it is not an accumulated orbit-error guarantee. The frozen host remains an adopted control, without tidal re-equilibration or population response.

## Occupation and path diagnostics

The trajectory gives the stopped residence and normalized BH density exposure,

```math
C_0(H)=\int_0^{H\wedge\tau} \mathbf{1}_{r(t)<r_\rho}\,dt,
\qquad
X_0(H)=\frac{1}{H}\int_0^{H\wedge\tau}
\frac{\rho_{\rm BH}(r(t))}{\rho_{\rm BH}(0)}\,dt.
```

The declared horizon remains the denominator when a path stops early. All entry and exit roots, radial turning points and stopping causes are retained.

At each quadrature point the actual velocity is rotated into the local radial/tangential frame. The original anisotropic BH number distributions and the [finite-angle encounter domain](finite-angle-encounters.md) then supply the rates, drift and diffusion. The baseline uses `b_max=0.1 L(r)` and the same closest-approach floor. Recorded diagnostics are

```math
\mathcal H_{\rm el}=\int\Gamma_{\rm el}(r(t),v(t))\,dt,
\qquad
\mathcal H_{\rm exc}=\int\Gamma_{\rm exc}(r(t),v(t))\,dt,
```

```math
B_A=\frac{1}{v_0}\int\|\mathbf A(r(t),v(t))\|\,dt,
\qquad
B_D=\frac{1}{v_0^2}\int\operatorname{tr}\mathsf D(r(t),v(t))\,dt.
```

These are budgets on a prescribed collisionless path. They are not measured event counts, capture probabilities, velocity changes or corrections to residence. In particular, the drift norm integral discards vector cancellation, and the diffusion trace integral does not include the changed positions and velocities of a stochastic path. A small budget alone does not establish accurate residence near a region boundary.

## Independent checks

The calculation refines the orbit step and tolerances, measures energy and angular-momentum conservation, and verifies residence independently from the spherical radial equation

```math
\dot r^2=2[E+\psi(r)]-\frac{\ell^2}{r^2}.
```

Numerically grazing contacts use the radial curvature to distinguish an inner apocenter touch from an outer pericenter touch; stationary and circular paths on the strict boundary contribute no residence. The boundary tolerance is refined with the orbit controls. Analytical checks include contacts at and just before or after the terminal time.

The independent integral uses the actual initial phase and finite horizon. A sine-squared change of radial variable resolves the turning-point endpoints. It does not use the Cartesian orbit's crossing times.

Encounter-budget quadrature splits at region crossings and radial turning points, including center passages of radial trajectories. The velocity quadrature is centered on relative velocity to resolve the close-passage flux near zero relative speed. Its integration ball contains the full shifted field distribution, which is evaluated at its original field velocity and vanishes outside the host support. Separate refinement checks the time quadrature, velocity quadrature, number density and diffusion tensor. The allowed plus excluded rate is checked against the incoming parent flux. The numerical thresholds are declared in [the protocol](../analysis/residence/protocol.json); unresolved records remain distinguishable from accepted calculations.

## Recorded findings

All twelve paths pass the registered checks and reach their horizon without a validity stop. The largest independent residence discrepancy divided by the horizon is below `8e-10`; the largest energy drift divided by `v_dyn²` is below `8e-9`. These measure numerical consistency, not physical host uncertainty.

| Host | Horizon (Myr) | Residence: radial / tangential | Allowed rate integral: radial / tangential |
|---|---:|---:|---:|
| 47 Tuc | 0.19 | 0.3025 / 0.2783 | 0.998 / 0.834 |
| NGC 1851 | 0.0818 | 0.2627 / 0.2753 | 0.142 / 0.137 |
| NGC 3201 | 1.67 | 0.3161 / 0.2910 | 0.0728 / 0.0599 |
| Omega Centauri | 0.0208 | 0.2719 / 0.2887 | 3.53 / 3.44 |
| M4 | 0.242 | 0.2625 / 0.2750 | 0.0488 / 0.0469 |
| NGC 6397 (standard) | 0.0224 | 0.2636 / 0.2765 | 0.0258 / 0.0248 |

The launch-direction ordering changes with host and finite-horizon phase. Similar residence fractions accompany substantially different encounter integrals. The shared horizon in dynamical units is not a shared physical duration, and the rate integral includes that duration. These findings establish the baseline and the range of path sampling that the collisional comparison must resolve; they do not establish a diffusion error.

## Connection to the residence hypothesis

The reference paths establish the event-location and occupation calculations needed by both stochastic descriptions. Their encounter budgets indicate where local-law evaluation and event sampling must be resolved. The next comparison must apply the matched drift/diffusion or full increments along evolving paths, include the same treatment of lighter species and excluded channels, and quantify time-step and ensemble uncertainty. The first-response prediction is then tested against the independently computed full encounter result.

The collisionless control does not test that response prediction, calibrate friction, or establish when diffusion is adequate. It gives a verified dynamical baseline against which those effects can be measured.

[Recorded comparison](../results/collisionless-residence/comparison.csv) · [Reproduction](../analysis/residence/README.md) · [Hypothesis](transport-hypothesis.md)
