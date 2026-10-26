# Step 1: merger-conditioned recoil transport

**Researcher: Venkata Hemanth Athota. Status: host and local coefficient calculations completed for one NGC 3201 realization; trajectory implementation and physical calibration pending.**

## Physical question

Given a remnant-producing merger and a consistent stellar environment, what probability distribution of remnant positions and velocities follows the recoil? Which trajectories remain within the domain where the transport approximation applies?

The output is a time-resolved stochastic trajectory ensemble, or an equivalent phase-space transport solution with its survival probability. A return time is not the complete output. [Step 2](step2-exposure.md) defines finite-region entry, residence, and environmental exposure on these histories.

## Host and coordinate conventions

Use physical units and translating, nonrotating cluster axes. Let the cluster reference position in an inertial frame be $\mathbf R(t)$ and the remnant state relative to it be $(\mathbf r,\mathbf v)$. The general equations describe conditional inputs. The [defined calculation](local-transport-calculations.md) supplies one frozen anisotropic multimass host and evaluates local coefficients, with remaining trajectory choices explicitly unassigned.

For species $j$ of individual mass $m_j$, define a **number** phase-space distribution $F_j(\mathbf r,\mathbf u,t\mid\Theta)$. Its integral over position and velocity is the expected species count. Host variables $\Theta$ include the mass components, distribution functions, boundary conditions, and any supplied evolution:

```math
n_j(\mathbf r,t)=\int F_j(\mathbf r,\mathbf u,t)\,d^3u,
\qquad
\rho_j=m_j n_j,
\qquad
\nabla^2\Phi_{\rm c}=4\pi G\sum_j\rho_j.
```

The tracked remnant's self-field is excluded. A center-of-mass approximation for an unresolved binary applies only on declared encounter scales; do not count its mass again through separate components. A continuous mass spectrum can be represented by specified bins or a mass integral with a declared normalization.

The same $F_j$ determines local streaming and the velocity-dispersion tensor:

```math
\overline{\mathbf u}_j=\frac{1}{n_j}\int \mathbf u F_j\,d^3u,
\qquad
S_{j,ik}=\frac{1}{n_j}\int
(u_i-\overline u_{j,i})(u_k-\overline u_{j,k})F_j\,d^3u.
```

These moments are defined where $n_j>0$. Rotation, velocity anisotropy, and species-dependent concentration are allowed. A density profile and one dispersion do not in general specify the full velocity distribution.

Poisson's equation requires boundary conditions and does not evolve the host. A frozen host requires a validity interval; an evolving host requires consistent $F_j(t)$ and $\Phi_{\rm c}(t)$. The remnant is a test object in that supplied background. Host backreaction, collective response, and energy feedback require additional equations.

For external potential $\Phi_{\rm ext}$, the cluster-relative external acceleration is

```math
\mathbf a_{\rm ext,rel}
=-\nabla\Phi_{\rm ext}(\mathbf R+\mathbf r,t)-\ddot{\mathbf R}(t).
```

If the reference orbit obeys $\ddot{\mathbf R}=-\nabla\Phi_{\rm ext}(\mathbf R,t)$, this equals the differential tidal acceleration. No Coriolis terms occur in nonrotating axes; a rotating frame would require its own inertial terms.

The [outer-boundary calculation](outer-boundary.md) evaluates this differential force in an explicit Galactic benchmark and separates an isolated-control accuracy limit from a circular tidal-opening scale. The [Galactic reference-orbit stage](galactic-orbit.md) then supplies time-dependent external forcing and a conservative force bound. Neither stage supplies a remnant trajectory, a tidally re-equilibrated host or an irreversible-escape rule.

## Joint merger-conditioned initial state

Let elapsed time zero denote the merger:

```math
\mathbf r(0)=\mathbf r_{\rm m},
\qquad
\mathbf v(0)=\mathbf v_{\rm COM}+\mathbf v_{\rm k},
\qquad
m_{\rm r}(0)=m_{\rm rem}.
```

The progenitor's center-of-mass velocity is distinct from its internal binary orbital velocity. Rotate the recoil vector from the binary basis into the cluster basis before addition. This uses nonrelativistic cluster kinematics; a relativistic merger calculation supplies the remnant mass and recoil.

Write the joint launch law as

```math
\pi_0(z_0\mid\Theta,\mathcal D,\mathcal M),
\qquad
z_0=(\mathbf r_{\rm m},\mathbf v_{\rm COM},\mathbf v_{\rm k},m_{\rm rem}).
```

Here $\mathcal D$ denotes observational or simulation constraints and $\mathcal M$ the merger-selection rule and relevant assembly information. Neither independence among launch variables nor independence from the host is assumed. The initial angular momentum $\mathbf r_{\rm m}\times(\mathbf v_{\rm COM}+\mathbf v_{\rm k})$ need not vanish.

Do not impose a retained-only sample without defining its escape criterion and retention interval. Instantaneous negative cluster-relative energy does not establish permanent retention in an evolving tidal environment. Define the base transport process before selecting on future survival; retained-path statistics require conditioning or reweighting the full path law. Changing the launch distribution alone while leaving subsequent stochastic evolution unchanged is insufficient. A noncentral launch need not cross the exact center.

## Transport equation

Set $\mathbf a_g=-\nabla\Phi_{\rm c}+\mathbf a_{\rm ext,rel}$. The kinetic description uses the **Itô** stochastic differential equations

```math
d\mathbf r=\mathbf v\,dt,
\qquad
d\mathbf v=(\mathbf a_g+\mathbf A)\,dt+\mathbf B\,d\mathbf W,
\qquad
\mathbf D=\mathbf B\mathbf B^{\mathsf T}.
```

The three Wiener components are independent with $\mathbb E[dW_i\,dW_k]=\delta_{ik}\,dt$. After subtracting smooth gravitational acceleration, the collision coefficients are

```math
A_i=\lim_{\Delta t\to0}
\frac{\mathbb E[\Delta v_i^{\rm coll}\mid \mathbf r,\mathbf v,t]}{\Delta t},
\qquad
D_{ik}=\lim_{\Delta t\to0}
\frac{\mathbb E[\Delta v_i^{\rm coll}\Delta v_k^{\rm coll}\mid \mathbf r,\mathbf v,t]}{\Delta t}.
```

The limits belong to a coarse-grained diffusion approximation. $\mathbf A$ has units of acceleration and $\mathbf D$ of velocity squared per time. $\mathbf A$ is the complete Itô first-moment drift; do not append an additional noise-induced drift.

Between discrete transitions, $m_{\rm r}$ is fixed on each trajectory. For a given remnant mass and host, let $p(\mathbf r,\mathbf v,t\mid m_{\rm r},\Theta)$ denote the phase-space probability density. Its interior evolution is

```math
\partial_t p
=-\partial_{r_i}(v_i p)
-\partial_{v_i}[(a_{g,i}+A_i)p]
+\frac12\partial_{v_i}\partial_{v_k}(D_{ik}p).
```

Repeated Cartesian indices are summed. If remnant masses vary across launches, condition on mass or retain it as a state label; a single mass-averaged drift is not equivalent. Initial distributions, boundary conditions, and stopping rules are necessary to define a solution.

### A candidate weak-scattering closure

For local, uncorrelated weak deflections with a leading Coulomb logarithm, define number-normalized Rosenbluth potentials at fixed position and time:

```math
h_j(\mathbf v)=\int
\frac{F_j(\mathbf r,\mathbf u,t)}{|\mathbf v-\mathbf u|}\,d^3u,
\qquad
g_j(\mathbf v)=\int
F_j(\mathbf r,\mathbf u,t)|\mathbf v-\mathbf u|\,d^3u.
```

The corresponding gravitational coefficients are

```math
A_i=4\pi G^2\sum_j m_j(m_{\rm r}+m_j)\ln\Lambda_j\,
\partial_{v_i}h_j,
\qquad
D_{ik}=4\pi G^2\sum_j m_j^2\ln\Lambda_j\,
\partial_{v_i}\partial_{v_k}g_j.
```

The number normalization is essential: a mass-normalized distribution has different explicit mass factors. The Coulomb logarithm $\ln\Lambda_j=\ln(b_{\max,j}/b_{\min,j})$ describes a declared encounter-scale range. It is held fixed during each velocity integral and differentiation. Arbitrary velocity-dependent cutoffs require a consistent derivation. No cutoff values or prescription are chosen.

The coefficients retain streaming and anisotropy through $F_j$. The sign follows from

```math
\partial_{v_i}h_j=-\int F_j\frac{w_i}{w^3}\,d^3u,
\qquad
\partial_{v_i}\partial_{v_k}g_j
=\int F_j\left(\frac{\delta_{ik}}w-\frac{w_iw_k}{w^3}\right)d^3u,
\qquad \mathbf w=\mathbf v-\mathbf u.
```

For nonnegative distributions and positive Coulomb logarithms, $\mathbf D$ is symmetric and positive semidefinite. The formulas do not impose a Maxwellian background. They are established kinetic ingredients, not validated predictions for a chosen cluster. [Rosenbluth, MacDonald & Judd (1957)](https://doi.org/10.1103/PhysRev.107.1) develops the underlying collision integrals; [Tep, Fouvry & Pichon (2022), equations 1–2](https://arxiv.org/html/2201.03985) gives the gravitational coefficients and an anisotropic application.

Local weak scattering excludes strong events, long-range correlations, and collective effects. Direct kicked-black-hole calculations show limits of simplified friction prescriptions and individual energy gains from scattering; these motivate testing the closure rather than assuming it is sufficient. [Webb et al. (2018)](https://arxiv.org/abs/1711.09100)

## Energy accounting and motion through the core

For fixed remnant mass between transitions, define specific cluster-relative energy

```math
\epsilon=\frac12|\mathbf v|^2+\Phi_{\rm c}(\mathbf r,t).
```

Itô's formula gives

```math
d\epsilon=
\left[\partial_t\Phi_{\rm c}
+\mathbf v\cdot\mathbf a_{\rm ext,rel}
+\mathbf v\cdot\mathbf A+\frac12\operatorname{tr}\mathbf D\right]dt
+\mathbf v^{\mathsf T}\mathbf B\,d\mathbf W.
```

The positive diffusion contribution is required. Individual trajectories need not lose energy monotonically. This energy is not a conserved tidal escape criterion. Entry into the core changes no position or velocity by prescription: motion continues.

## Validity, stopping, and probability conservation

Declare a common physical horizon $T$ and the stopping time

```math
T_*=\min(T,\tau_{\rm esc},\tau_{\rm strong},\tau_{\rm valid}).
```

The escape test must match the external field and intended physical question. A strong-encounter stop requires an actual event detector, resolved objects, or a justified hazard. An absent detector cannot be interpreted as absence of strong encounters. The validity rule must identify where the local-scattering and test-remnant approximations cease to apply.

Drift/diffusion and explicit close encounters must partition deflections consistently; they cannot apply the same interaction twice. A strong encounter that changes velocity, mass, or binary membership needs a separate transition model. That model is outside these two steps; stop and report the unresolved event.

If trajectories are removed on escape or another declared stop, use the unnormalized alive density. Its integral is

```math
\int p_{\rm alive}(\mathbf r,\mathbf v,t\mid m_{\rm r},\Theta,\mathcal D,\mathcal M)\,d^3r\,d^3v
=\Pr(t<T_*\mid m_{\rm r},\Theta,\mathcal D,\mathcal M).
```

This density is conditional on a fixed remnant mass. To average over masses, integrate these conditional densities against the initial mass marginal from the joint launch law; do not integrate a fixed-mass conditional density over mass without that weight. The interior Fokker–Planck operator must be supplemented by the corresponding absorbing boundary or justified killing mechanism. A history-dependent stopping rule requires an augmented state or a justified history-dependent treatment; it does not automatically give a closed equation in position and velocity alone. Do not renormalize surviving paths and report them as the full population. Step 2 distinguishes stopped residence from complete-horizon residence.

## Uncertainty propagation

For any trajectory observable $O$,

```math
\mathbb E[O\mid\mathcal D,\mathcal M]
=\int d\Theta\,p(\Theta\mid\mathcal D,\mathcal M)
\int dz_0\,\pi_0(z_0\mid\Theta,\mathcal D,\mathcal M)
\mathbb E_{\rm path}[O\mid\Theta,z_0,\mathcal D,\mathcal M].
```

The inner expectation covers Wiener forcing and any additional randomness in a supplied stopping-event model. Physical scattering, launch uncertainty, and host uncertainty are distinct. This conditional factorization imposes no independence. Merger selection can change host weights. Alternative model families require explicit conditioning or model weights; none are inferred here.

## Handoff and acceptance requirements

Step 1 must deliver physical times, positions, velocities, fixed remnant-mass labels, host/launch realization identifiers, model-validity flags, stopping times and causes, and either ensemble weights or the normalization of a transport-density solution. The record must preserve the source of each input, coordinate and unit conventions, random seeds if used, and numerical error controls.

Before physical interpretation, a numerical implementation must demonstrate dimensional consistency, force–potential agreement, positive-semidefinite diffusion, the declared Itô convention, probability accounting, appropriate collisionless and analytical stochastic limits, and convergence of trajectory distributions and stopping events. Exact conservation is required only in a static conservative limit. The stochastic energy identity supplies a separate consistency check when scattering is present. These are acceptance requirements, not executed tests or evidence of calibrated physics.

[Step 2](step2-exposure.md) · [Problem statement](problem.md) · [Calculation handoff](calculation-handoff.md) · [Input register](../formulation/inputs.json)
