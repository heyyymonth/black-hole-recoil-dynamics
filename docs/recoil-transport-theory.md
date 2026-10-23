# Recoil transport and companion access: proposed formulation

**Researcher: Venkata Hemanth Athota. Formulation revision: 23 October 2026.**

This document defines the extension motivated by the completed orbital, encounter, potential-shape, and target-count studies. It contains no new integrations, cluster fits, parameter estimates, or experimental results. No host family or numerical values are selected here. The existing Plummer results remain the collisionless limiting benchmark.

The question is: **given where a merger happens and the environment around it, how does the surviving black hole move, and how does that motion overlap with available binary companions?** Step 1 must provide position and velocity histories, with uncertainty and explicit termination rules. A return time alone is insufficient input to the later passage calculation.

## 1. What changes in step 1

| Ingredient | Completed orbital benchmark | Proposed input or output |
|---|---|---|
| Launch | At the center, with recoil as the entire initial velocity | Merger position, progenitor center-of-mass velocity, and recoil vector in a common frame |
| Host gravity | Static spherical Plummer potential | A potential consistent with the mass distributions of declared species; external acceleration and time dependence where supplied |
| Scatterers | No collisional effects | Species masses and local velocity distributions determine both mean velocity drift and stochastic diffusion |
| Motion | One radial deterministic trajectory | Three-dimensional trajectories or a phase-space probability distribution |
| Return | Next exact-center crossing | First entry or re-entry into a finite, declared region; motion continues afterward |
| Downstream output | Residence and density exposure | Joint position–velocity occupation, physical times, event status, and conditional target access |

These additional inputs are not automatically additional fitted scalar parameters. A velocity distribution, evolving host, or target history is a function requiring a physical model and evidence. Writing it symbolically does not establish its value or validity.

## 2. Host matter and reference frame

Use physical units throughout this extension. Let the cluster reference position in an inertial frame be $\mathbf R(t)$; use translating **nonrotating** axes. The remnant's cluster-relative position and velocity are $\mathbf r$ and $\mathbf v$.

For mass species $j$, define a **number** phase-space density $F_j(\mathbf r,\mathbf u,t\mid\Theta)$, where $\Theta$ describes the host. Thus $F_j\,d^3r\,d^3u$ is an expected object count, not a mass or a probability. Define

```math
n_j=\int F_j\,d^3u,\qquad
\rho_j=m_j n_j,\qquad
\rho_{\rm c}=\sum_j\rho_j,\qquad
\nabla^2\Phi_{\rm c}=4\pi G\rho_{\rm c}.
```

Boundary conditions are also required for the potential. The tracked remnant's self-field is excluded. An unresolved binary can enter the weak-scattering background as one center of mass only on encounter scales where that approximation applies; its mass must not also be counted through separate components in the same description.

The same species distributions specify streaming and velocity dispersion:

```math
\overline{\mathbf u}_j=\frac{1}{n_j}\int\mathbf u F_j\,d^3u,
\qquad
S_{j,ik}=\frac{1}{n_j}\int
(u_i-\overline u_{j,i})(u_k-\overline u_{j,k})F_j\,d^3u.
```

These moments are defined where $n_j>0$. They allow rotation, velocity anisotropy, distinct mass components, and mass segregation. Density and one scalar dispersion do not generally determine the full velocity distribution required below. The potential and collision coefficients must refer to the **same host and epoch**, rather than independent profiles assembled from incompatible models.

Let $\Phi_{\rm ext}$ be the external potential. The relative external acceleration is

```math
\mathbf a_{\rm ext,rel}(\mathbf r,t)
=-\nabla\Phi_{\rm ext}(\mathbf R+\mathbf r,t)-\ddot{\mathbf R}(t).
```

If the reference trajectory follows $\ddot{\mathbf R}=-\nabla\Phi_{\rm ext}(\mathbf R,t)$, this becomes the tidal difference

```math
\mathbf a_{\rm tide}
=-\nabla\Phi_{\rm ext}(\mathbf R+\mathbf r,t)
+\nabla\Phi_{\rm ext}(\mathbf R,t).
```

There are no Coriolis terms in these nonrotating axes. A rotating coordinate system would require its own inertial terms. Internal stellar rotation is already permitted through $F_j$.

This is a test-remnant description in a supplied background. Poisson's equation does not evolve that background. A frozen host needs a declared validity interval; an evolving host needs an externally supplied or self-consistently evolved $F_j(t)$. Host response, collective effects, and energy feedback are not obtained merely by allowing a time argument.

## 3. A merger-conditioned launch

Choose elapsed time $t=0$ at the remnant-producing merger. The initial conditions are

```math
\mathbf r(0)=\mathbf r_{\rm m},\qquad
\mathbf v(0)=\mathbf v_{\rm COM}+\mathbf v_{\rm k},\qquad
m_{\rm r}(0)=m_{\rm rem}.
```

The recoil vector must be transformed from the binary frame into the cluster frame before addition. This velocity addition uses the nonrelativistic cluster-motion approximation; relativistic merger and recoil calculations provide the initial remnant mass and kick. The progenitor's center-of-mass velocity is not its internal binary orbital velocity.

Let $\mathcal M$ denote the merger selection and relevant assembly information. Specify a **joint** launch law

```math
\pi_0(\mathbf r_{\rm m},\mathbf v_{\rm COM},\mathbf v_{\rm k},m_{\rm rem}
\mid\mathcal M,\Theta).
```

No independence between position, kick, remnant mass, and host is imposed. In general the specific angular momentum $\mathbf r_{\rm m}\times(\mathbf v_{\rm COM}+\mathbf v_{\rm k})$ is nonzero. The trajectory need not be radial or pass through the exact center. If the sample is conditioned on retention, define the escape test and time interval that constitute retention; initial negative cluster-relative energy is not a permanent-retention test in an evolving tidal environment.

## 4. Mean scattering and random deflections

Define $\mathbf a_g=-\nabla\Phi_{\rm c}+\mathbf a_{\rm ext,rel}$. A proposed local diffusion description is the **Itô** system

```math
d\mathbf r=\mathbf v\,dt,
\qquad
d\mathbf v=[\mathbf a_g+\mathbf A]\,dt+\mathbf B\,d\mathbf W,
\qquad
\mathbf D=\mathbf B\mathbf B^{\mathsf T}.
```

$\mathbf W$ has three independent standard Wiener components. $\mathbf A$ is the collisional mean velocity drift, including dynamical friction; $\mathbf D$ is the collisional second velocity increment per unit time. Smooth gravity is excluded from these increment definitions:

```math
A_i=\lim_{\Delta t\to0}
\frac{\mathbb E[\Delta v_i^{\rm coll}\mid\mathbf r,\mathbf v,t]}{\Delta t},
\qquad
D_{ik}=\lim_{\Delta t\to0}
\frac{\mathbb E[\Delta v_i^{\rm coll}\Delta v_k^{\rm coll}\mid\mathbf r,\mathbf v,t]}{\Delta t}.
```

The limits refer to a coarse-grained diffusion approximation, not resolving an individual gravitational encounter over an infinitesimal physical time. $\mathbf A$ has units of acceleration and $\mathbf D$ of velocity squared per time. The corresponding phase-space density $p$ obeys

```math
\frac{\partial p}{\partial t}
=-\frac{\partial(v_i p)}{\partial r_i}
-\frac{\partial[(a_{g,i}+A_i)p]}{\partial v_i}
+\frac12\frac{\partial^2(D_{ik}p)}{\partial v_i\partial v_k},
```

with repeated Cartesian indices summed. This convention fixes the factor of one half; $\mathbf D=\mathbf B\mathbf B^{\mathsf T}$, not $2\mathbf B\mathbf B^{\mathsf T}$. Once $\mathbf A$ is defined as the complete Itô first moment, a second noise-induced drift must not be added to it.

### An explicit candidate closure, with stated limits

For local, uncorrelated weak deflections, the leading-logarithm Rosenbluth construction gives both coefficients from the same field population. At fixed position and time, define

```math
h_j(\mathbf v)=\int\frac{F_j(\mathbf r,\mathbf u,t)}{|\mathbf v-\mathbf u|}\,d^3u,
\qquad
g_j(\mathbf v)=\int F_j(\mathbf r,\mathbf u,t)|\mathbf v-\mathbf u|\,d^3u.
```

For the **number normalization** used here,

```math
A_i=4\pi G^2\sum_j m_j(m_{\rm r}+m_j)\ln\Lambda_j
\frac{\partial h_j}{\partial v_i},
\qquad
D_{ik}=4\pi G^2\sum_j m_j^2\ln\Lambda_j
\frac{\partial^2g_j}{\partial v_i\partial v_k}.
```

Here $\ln\Lambda_j=\ln(b_{\max,j}/b_{\min,j})$ represents a declared range of encounter impact parameters. No cutoff values are selected. It is held fixed during each velocity integral and differentiation; arbitrary velocity-dependent cutoffs require a consistent rederivation. The sign of the drift follows from $\nabla_v h_j=-\int F_j(\mathbf v-\mathbf u)/|\mathbf v-\mathbf u|^3\,d^3u$. The diffusion tensor is symmetric and positive semidefinite for nonnegative distributions and positive Coulomb logarithms. Anisotropy and streaming are retained through $F_j$; no Maxwellian is imposed.

This is a candidate closure, not a claim of validated friction in any selected cluster. Its basis is [Rosenbluth, MacDonald & Judd (1957)](https://doi.org/10.1103/PhysRev.107.1); the gravitational coefficients and anisotropic application are given by [Tep, Fouvry & Pichon (2022), equations 1–2](https://arxiv.org/html/2201.03985). Their mass-normalized distribution must be converted to the number normalization above. Local weak scattering omits strong encounters and collective or long-range correlations. The appropriate regime, cutoffs, and accuracy remain to be established before use.

In particular, [Webb et al. (2018)](https://arxiv.org/abs/1711.09100) find departures from simplified friction predictions for kicked black holes and stochastic energy gains in their star-cluster experiments. Adding a deterministic drag term alone would not resolve all of those effects. Their results do not establish that the closure proposed here is sufficient.

### Motion continues after return

For constant remnant mass between discrete interactions, define cluster-relative specific energy

```math
\epsilon=\frac12|\mathbf v|^2+\Phi_{\rm c}(\mathbf r,t).
```

Itô's rule yields

```math
d\epsilon=
\left[\partial_t\Phi_{\rm c}
+\mathbf v\cdot\mathbf a_{\rm ext,rel}
+\mathbf v\cdot\mathbf A+\frac12\operatorname{tr}\mathbf D\right]dt
+\mathbf v^{\mathsf T}\mathbf B\,d\mathbf W.
```

The diffusion contribution cannot be dropped. An individual trajectory can gain energy; the energy above is not generally conserved and is not by itself a tidal escape criterion. No rule resets position or velocity to zero at core entry. The framework permits continued excursions and stochastic wandering without assuming that this closure correctly predicts their magnitude.

## 5. Define return, residence, and the information passed downstream

Specify a finite physical region $\mathcal C(t)$, such as $|\mathbf r|<r_c(t)$. Its radius is an operational definition, not automatically the same as an observed projected core radius or a model scale length.

For a launch inside the region, define first outward exit and subsequent inward re-entry:

```math
t_{\rm out}=\inf\{t>0:\mathbf r(t)\notin\mathcal C(t)\},
\qquad
t_{\rm re}=\inf\{t>t_{\rm out}:\mathbf r(t)\in\mathcal C(t)\}.
```

Only crossings that change inside/outside membership count; a tangential touch is not an exit. Crossing conventions must specify the boundary's own motion. A launch outside has a first-entry event; a remnant that never leaves has no re-entry event. A moving boundary overtaking the remnant is distinct from inward remnant motion. An implementation can declare separated inner and outer event surfaces to suppress numerical recrossings, but must report their effect rather than silently change the region definition.

Let $T$ be a common physical horizon. Follow the process until $T_*\leq T$, the earliest of that horizon, a declared escape, a strong encounter requiring an unavailable transition model, or loss of the chosen approximation's validity. Record which outcome stops it. Define

```math
C_T=\int_0^{T_*}\mathbf 1_{\{\mathbf r(t)\in\mathcal C(t)\}}\,dt,
\qquad
J_{j,T}=\int_0^{T_*}\rho_j[\mathbf r(t),t]\,dt.
```

$C_T/T$ is residence accumulated during the modeled interval divided by the full horizon. It is not an assertion of zero physical residence after an unresolved strong encounter or validity cutoff. Complete-horizon residence cannot be inferred for such stopped trajectories. Escape is a distinct physical outcome; an unobserved return before $T$ is censored rather than proven impossible.

A normalized density exposure may use a declared positive reference $\rho_{{\rm ref},j}(t)$:

```math
I_{j,T}=\int_0^{T_*}\frac{\rho_j[\mathbf r(t),t]}{\rho_{{\rm ref},j}(t)}\,dt.
```

Neither exposure is an encounter count. Preserve the joint position–velocity occupation measure

```math
\mu_T(\mathcal A)=\int_0^{T_*}
\mathbf 1_{\{(\mathbf r(t),\mathbf v(t))\in\mathcal A\}}\,dt
```

along with time-resolved trajectories and event histories. For evolving targets, this occupation measure without time labels is insufficient; the chronological trajectory or a time-augmented measure is required. Later access calculations need velocities and timing, not just the average radius, residence fraction, or mean return time.

Host uncertainty, launch uncertainty, and physical scattering are distinct. For an observable $O$, their proposed propagation has the form

```math
\mathbb E[O\mid\mathcal D,\mathcal M]
=\int d\Theta\,p(\Theta\mid\mathcal D,\mathcal M)
\int dz_0\,\pi_0(z_0\mid\Theta,\mathcal D,\mathcal M)
\mathbb E_{\mathbf W}[O\mid\Theta,z_0,\mathcal M].
```

$z_0$ collects the launch variables and $\mathcal D$ denotes observations or simulation constraints. This is a conditional factorization, not an assumption of independent errors. Merger selection can change the host weights. Different host families would require separate model conditioning or stated model weights; no posterior or weights have been estimated here.

## 6. Connection to binary companions

The scatterer distributions $F_j$ determine gravity, drift, and diffusion. The target distribution below instead describes **available binary-black-hole systems**. It excludes the consumed progenitor binary and the tracked remnant's own system. The same physical binary may contribute to background gravity and qualify as a target, but its encounters must be partitioned consistently between weak scattering and explicitly modeled close passages.

Let $Z_t=(\mathbf r_t,\mathbf v_t,m_{\rm r})$ and let $\mathcal H_t$ contain relevant merger, host, and encounter history. Define a conditional target number intensity $\mathcal F_B(\mathbf y,\mathbf u,\zeta,t\mid Z_t,\mathcal H_t)$, where $\zeta$ includes binary masses and internal orbital properties. Its integral over position, velocity, and $\zeta$ is the conditional expected available target count. A local kinetic passage model is

```math
\Gamma_B(t\mid\mathcal H_t)=
\int d^3u\,d\zeta\,
\mathcal F_B(\mathbf r_t,\mathbf u,\zeta,t\mid Z_t,\mathcal H_t)
\Sigma_B(w,\zeta;Z_t)\,w,
\qquad w=|\mathbf v_t-\mathbf u|.
```

$\Sigma_B$ is a cross section for a **declared passage**, not an unspecified binary-acquisition or merger cross section. Its definition must remain compatible with the meaning of $\zeta$, the weak-scattering cutoff, and the velocity measured before mutual focusing.

If $\Gamma_B$ is adopted as a valid conditional counting intensity, the expected number up to the stopping time is

```math
\mathbb E[N_{\rm pass}(T_*)\mid\mathcal M]
=\mathbb E\left[\int_0^{T_*}\Gamma_B(t\mid\mathcal H_t)\,dt
\,\middle|\,\mathcal M\right].
```

This does not imply an unconditional Poisson law, independent targets, or a merger probability. Depletion, correlations with remnant motion, target formation and destruction, and changes after encounters belong in the conditional history and transition model. The earlier availability decomposition remains

```math
\mathbb E[\Gamma_B(t)\mid\mathcal M]
=\sum_n P[N_B(t)=n\mid\mathcal M]
\mathbb E[\Gamma_B(t)\mid N_B(t)=n,\mathcal M].
```

It cannot generally be factorized into the product of an independently averaged target count and orbital kernel. In a sparse single realization, object-resolved encounter calculations or a justified ensemble intensity are alternatives. A list of two target positions is not itself a calibrated smooth density to insert into this expression. Observed luminous-star kinematics do not supply an inventory of hidden BBH targets.

Close encounters capable of altering the remnant's trajectory or identity need a separate few-body transition rule. Deflections already represented by $\mathbf A,\mathbf D$ must not also be applied as explicit kicks. Until that partition and transition rule are specified, the transport description covers motion before the first unresolved strong encounter. It does not predict the next merger.

## 7. Recovering the executed Plummer benchmark

Set drift, diffusion, and external relative acceleration to zero; choose a static Plummer potential; and restore a central launch with no progenitor center-of-mass motion:

```math
\mathbf A=\mathbf D=\mathbf a_{\rm ext,rel}=0,
\qquad
\Phi_{\rm c}=-\frac{GM}{\sqrt{r^2+b^2}},
\qquad
\mathbf r_{\rm m}=\mathbf v_{\rm COM}=0.
```

The trajectory is radial. With $\mathbf q=\mathbf r/b$, $s=t/t_d$, and $t_d=\sqrt{b^3/(GM)}$,

```math
\frac{d^2\mathbf q}{ds^2}=-\frac{\mathbf q}{(1+q^2)^{3/2}}.
```

This recovers the signed-coordinate equation in [the benchmark statement](problem.md). The density ratio becomes $(1+q^2)^{-5/2}$. Choosing $\mathcal C=\{r<cb\}$ and $T=Ht_d$ gives $C_T=t_d C(H,c)$ and $I_{\rho,T}=t_d I(H)$ when no stopping event occurs. Exact-center return remains the special radial benchmark; it is not the same observable as finite-region re-entry. Reinstating the specified stationary target bath and passage cross section recovers the earlier velocity-weighted integral.

## 8. How the existing findings constrain this formulation

| Completed study | Requirement imposed on the extension | What remains unknown |
|---|---|---|
| E001/E002: return and exposure | Preserve trajectories and residence, rather than impose permanent central rest after return | Collisional paths and wandering amplitudes |
| E003: velocity-weighted passages | Retain the joint relative-velocity distribution and a defined passage cross section | Encounter calibration in an actual evolving population |
| E004: frozen host and target comparison | Use a consistent mass profile, potential, and species description | Live host response, drift and diffusion accuracy |
| E005: late-time counts | Preserve discrete availability, history, and merger conditioning | The joint target phase space during remnant displacement |

The proposed formulation leaves the host family, mass spectrum, velocity distributions, launch law, external orbit, time horizon, escape definition, encounter cutoffs, and target histories **unassigned**. Lowered-isothermal multimass models provide one possible distribution-function family, not observed potentials or a selected solution; see [Gieles & Zocchi (2015)](https://arxiv.org/abs/1508.02120). Public observations can constrain host fits, while simulation histories may supply otherwise unobserved remnant and target information. Their compatibility must be established before they are used together.

The immediate theoretical advance is a consistent definition of the inputs and observables needed to test the existing accessibility hypothesis. The equations are established kinetic ingredients organized for this question; no new force law, validated transport model, calibrated companion probability, or empirical discovery is claimed.
