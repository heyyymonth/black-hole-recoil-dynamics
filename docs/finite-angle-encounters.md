# Finite-angle encounters and matched local moments

The residence hypothesis compares finite gravitational deflections with a diffusion description derived from the same encounters. This calculation establishes that common local law in the six adopted cluster backgrounds. It evaluates encounter rates, velocity-change moments and samples of individual elastic increments. Remnant trajectories, residence distributions and the proposed residence correction have not been evaluated.

## Relation to the research question

The previous local coefficients quantify weak gravitational scattering per unit constant Coulomb logarithm. The matched comparison needs an explicit finite encounter domain and an event sampler. These are now specified together: the event rate, sampled deflections, drift and diffusion all use the same species distributions and impact-parameter bounds.

This makes a future difference between the two trajectory descriptions interpretable. Both can share the same local first two encounter moments while differing in the remaining information about individual velocity changes. The present results establish that local construction; they do not establish how much it changes the time spent near potential companions.

## Adopted local environment

The backgrounds are the verified reference reconstructions of 47 Tuc, NGC 1851, NGC 3201, Omega Centauri, M4 and the NGC 6397 standard fit. Their joint parameter rows, species masses and anisotropic number distributions are retained. The diagnostic remnant mass is 60 solar masses, and the local calculation selects BH species only. The bath is frozen and encounters are sampled independently; persistent object identities, reciprocal population response and binaries are separate parts of the proposed trajectory study.

Let `r_d` be the radius at which total BH mass density is half its central value. Define the local density scale and upper impact parameter by

```math
L(r)=\min\left[r_d,\frac{\rho_{\rm BH}(r)}{|d\rho_{\rm BH}/dr|}\right],
\qquad b_{\max}=\eta_b L(r),
\qquad \eta_b\in\{0.05,0.1,0.2\}.
```

At zero density derivative, `L=r_d`. The density radius is a defined calculation scale, distinct from an observed core radius or a BH half-mass radius. The three upper cutoffs are sensitivity controls. Their small fractions of the density scale do not by themselves establish that every encounter is local or impulsive.

At each of `r/r_d = 0, 1, 2`, the calculation evaluates zero remnant speed and speed equal to the number-weighted BH one-dimensional dispersion, with radial and tangential moving orientations. That dispersion is defined from one third of the full mean squared speed; the anisotropic distribution remains in the integrals. Combined with the three cutoffs, this gives 27 states per host and 162 across the comparison.

## Elastic domain and the excluded flux

For a field BH of mass `m_j`, write `M_pair=m_r+m_j`, the symmetric mass ratio as `eta_j=m_r m_j/M_pair^2`, and the relative speed as `w`. The leading parabolic radiation-loss estimate gives the capture scale

```math
r_{\rm cap}(w)=\frac{G M_{\rm pair}}{c^2}
\left[\frac{85\pi}{6\sqrt{2}}\eta_j
\left(\frac{c}{w}\right)^2\right]^{2/7}.
```

This scale follows the two-body radiation-loss calculation used by [O'Leary, Kocsis and Loeb (2009)](https://doi.org/10.1111/j.1365-2966.2009.14653.x). Their nuclear-cluster population predictions are not transferred to these hosts. The registered elastic control adopts

```math
r_{p,\min}=\max\left[100\frac{G M_{\rm pair}}{c^2},10r_{\rm cap}(w)\right],
\qquad b_{\min}^2=r_{p,\min}^2+
\frac{2G M_{\rm pair}r_{p,\min}}{w^2}.
```

The factors 100 and 10 are declared separation margins, not measured thresholds or a calibrated relativistic scattering model. The same velocity-dependent lower bound is retained inside every rate and moment integral. When the floor requires `b_min >= b_max`, that velocity/species combination has no allowed elastic annulus.

Incoming flux below the floor is recorded separately from the allowed elastic rate. It is an excluded-channel flux, not a capture or merger rate: the floor deliberately also excludes encounters outside the leading capture region. A future trajectory calculation must give this channel a common event treatment or stopping rule in both descriptions. No such events or trajectories have been simulated here.

## Rate, deflection and matched moments

For the number distribution `F_j`, the differential incoming rate is

```math
d\Gamma_j=F_j(\mathbf r,\mathbf u)w\,d^3u\,b\,db\,d\varphi,
\qquad \mathbf w=\mathbf v-\mathbf u.
```

Exact elastic Newtonian two-body scattering rotates the relative velocity by the deflection angle and changes the remnant velocity by

```math
\theta_j=2\arctan\frac{G(m_r+m_j)}{b w^2},
\qquad \boldsymbol\xi_j=
\frac{m_j}{m_r+m_j}(\mathbf w'-\mathbf w),
\qquad |\mathbf w'|=w.
```

Integrating the same allowed measure supplies

```math
\Gamma=\sum_j\int d\Gamma_j,\qquad
\mathbf A=\sum_j\int\boldsymbol\xi_j\,d\Gamma_j,\qquad
\mathsf D=\sum_j\int\boldsymbol\xi_j\boldsymbol\xi_j^{\mathsf T}\,d\Gamma_j.
```

`D` is the raw second jump moment per unit time. It is the tensor multiplying one half of the velocity Hessian in the diffusion generator; it is not the covariance of the conditional single-event sample. The implementation integrates impact parameter and azimuth analytically, and integrates field velocity with positive quadrature weights over the full local distribution-function support. Rates are reported per Myr, drift in `(km/s)/Myr`, and diffusion in `(km/s)^2/Myr`.

The sampler chooses a species and field-velocity node with probability proportional to its allowed incoming flux. It then samples `b^2` uniformly over the annulus and samples a uniform azimuth. Conditional increment means and raw second moments therefore match `A/Gamma` and `D/Gamma` for the same numerical bath. Velocity-quadrature refinement separately checks the approximation to the continuous distribution.

## Numerical evidence and interpretation

The registered sample contains 100,000 conditional increments per host, or 600,000 in total, at `r=r_d`, radial speed equal to the local BH dispersion, and `b_max=0.1L`. These are independent events at fixed states, not time-ordered encounter histories.

Each local record retains numerical acceptance and its refinement history. The protocol requires rate and moment refinement below 0.5%, number-density agreement below 0.2%, a positive-semidefinite diffusion tensor, the central zero-drift symmetry, and the stated separate tolerance for the deflection-scale diagnostic. The sampler compares all three mean increments and six independent raw second moments with their expected values, retaining standard errors. [The aggregate summary](../results/finite-angle-encounters/summary.json) gives the current accepted-state count; [the protocol](../analysis/encounters/protocol.json) is authoritative for the acceptance rules.

The recorded diagnostics include the rate times a local crossing time, the root-mean-square single-event increment divided by the BH dispersion, and the share of the diffusion trace supplied by relative velocities with `b90 > L`. These describe encounter frequency, increment size and the deflection length relative to the local environment. They identify assumptions that a trajectory calculation must carry forward. None alone is an established residence-error criterion or evidence that diffusion fails.

## First-response computation

An implemented numerical route avoids assuming a solved six-dimensional backward equation. For an N-step transition mixture `P_a=(1-a)P_D+aP_J`, choose a step K uniformly, replace that one diffusion transition by its jump counterpart, and retain diffusion before and after it. Then

```math
\mathbb E\left[N\left(O_{\rm replacement}-O_D\right)\right]
=\left.\frac{d}{da}\mathbb E_a O\right|_{a=0}.
```

This is a derivative of the discrete transition mixture. Exact finite-state enumeration and independent numerical differentiation test the implementation; they do not provide a physical cluster result. The estimator uses independent random streams per path and step, shared only within each baseline/replacement pair. This preserves independent paired samples even when a callback consumes a state-dependent number of random values.

Time-step and occupation-weight refinement are required before identifying this estimate with the continuous-generator response in the paper. Predicting the full jump endpoint from that first response remains an extrapolation to be tested independently, including sampling uncertainty and possible nonlinear response.

## Next dynamical calculation

The local law supplies the scattering input for the following controlled calculation:

1. Declare the launch state or joint launch law, spatial region, horizon, common lighter-species treatment, external forcing and stopping rules. Preserve the distinction between an isolated-force control boundary and physical escape.
2. Establish local-law evaluation along trajectories, including the spatial, velocity and time resolution needed by the compact density profiles. Carry excluded-channel flux and encounter-duration assumptions into the validity accounting.
3. Compute and record the diffusion residence and its first-response prediction with numerical and sampling uncertainties.
4. Independently evaluate the matched finite-encounter trajectories, preserving censored histories and their causes.
5. Compare mean residence, distribution shape, species exposure and termination probabilities across the adopted hosts. Test finite-population response and evolving backgrounds in separately controlled extensions.

The intended result is an accuracy statement for residence and exposure. The local event law and sampler are necessary inputs to that result; companion acquisition and subsequent merger probabilities additionally require their own interaction physics.

[Calculation code and reproduction](../analysis/encounters/README.md) · [Residence hypothesis](transport-hypothesis.md) · [Six-cluster backgrounds](six-cluster-comparison.md)
