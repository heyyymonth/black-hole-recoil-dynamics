# Step 2: residence and environmental exposure

**Researcher:** Venkata Hemanth Athota

**Status:** mathematical formulation and calculation specification. The [host calculation](local-transport-calculations.md) supplies consistent species densities for one realization. [Collisionless controls](collisionless-residence.md) evaluate path accumulation, crossings and residence for twelve diagnostic launches, with independent checks. Stochastic residence distributions, moving-host cases and their additional acceptance checks remain unexecuted.

## Physical question

Given a merger remnant's motion through a stellar cluster, when does it enter a declared dense region, how much time does it spend there, and which matter populations does it sample during that time?

Gravitational retention, entry into a region, and sustained residence are different observables. A remnant can remain bound while repeatedly leaving the region of interest. Entry does not imply rest, and a return time alone does not describe environmental exposure. This step defines these observables on the same position–velocity histories supplied by [Step 1](step1-transport.md), preserving their physical timing and uncertainty.

This step measures residence and exposure. It does not assign companions, encounter rates, binary-acquisition probabilities, or subsequent merger probabilities. The scientific scope and required inputs are set out in the [problem statement](problem.md) and [calculation handoff](calculation-handoff.md).

## Inputs and the modeled interval

Step 1 supplies a trajectory ensemble, or an equivalent transport description, conditioned on a consistent host, merger selection, and launch law. Each realization retains:

- Physical time, cluster-relative position and velocity, remnant mass, reference frame, and units.
- Host and launch labels, statistical weights, and any physical scattering realization needed to reproduce the path.
- Species-resolved mass densities at the same positions and epochs as the trajectory.
- The observation horizon, escape rule, validity interval, event definitions, and reason for any early stop.

Let the common physical horizon be $`T`$. Define the modeled endpoint by

```math
T_*=\min\{T,\tau_{\rm esc},\tau_{\rm strong},\tau_{\rm valid}\}.
```

The three stopping times refer respectively to a declared escape event, a strong encounter requiring an unavailable transition model, and the end of the transport description's validity. A stopping time is operational only when its event rule is supplied. In particular, a strong-encounter stop requires a resolved event detector or another justified event model; missing encounter information does not establish that no strong encounter occurs.

The stopping reason accompanies every accumulated quantity. Core entry or re-entry is not a stopping condition. Motion and accumulation continue until $`T_*`$.

## A finite region and its crossings

Define the region by a declared level-set function:

```math
\mathcal C(t)=\{\mathbf r:g(\mathbf r,t)<0\}.
```

The function $`g`$ is negative inside, positive outside, and zero on the boundary. A sphere is the possible convention $`g(\mathbf r,t)=|\mathbf r|-r_c(t)`$; neither spherical shape nor a radius is selected here. A projected observational core radius is not automatically a three-dimensional region boundary. The conversion and its uncertainty belong to the region specification.

For a launch strictly inside the region, define first exit and subsequent re-entry by

```math
\tau_{\rm out}=\inf\{t>0:g(\mathbf r_t,t)>0\},
\qquad
\tau_{\rm re}=\inf\{t>\tau_{\rm out}:g(\mathbf r_t,t)<0\}.
```

These definitions require reaching the opposite side. An isolated tangential contact with $`g=0`$ does not count as an exit. For a continuous trajectory, the infimum identifies the crossing boundary time; an implementation must establish the sign change rather than count a boundary contact alone. A re-entry delay, when observed, is $`\tau_{\rm re}-\tau_{\rm out}`$.

For a launch strictly outside, the relevant event is first entry:

```math
\tau_{\rm ent}=\inf\{t>0:g(\mathbf r_t,t)<0\}.
```

A launch on the boundary needs a declared initial-side convention. It may be classified from the first definite side reached, with the contact recorded separately; it must not automatically be recorded as a return. A remnant that never leaves during the modeled interval has no observed re-entry event. The convention $`\inf\varnothing=+\infty`$ is mathematical bookkeeping, not evidence that entry or return can never happen after observation ends.

For a differentiable boundary, the relative crossing motion is

```math
\frac{d}{dt}g(\mathbf r_t,t)
=\partial_t g(\mathbf r_t,t)
+\mathbf v_t\cdot\nabla g(\mathbf r_t,t).
```

This distinguishes boundary evolution from remnant motion. A region expanding over a remnant changes membership but is not necessarily inward remnant motion. Both the membership event and its cause must remain identifiable. Sign-changing crossings remain the definition even when the derivative vanishes at a crossing. Separated inner and outer event surfaces may control numerical recrossings, but their definitions and effects must be reported; they do not silently redefine residence in $`\mathcal C(t)`$.

## Residence and species-resolved exposure

Accumulated residence in the declared region is

```math
C_T=\int_0^{T_*}\mathbf{1}_{\{\mathbf r_t\in\mathcal C(t)\}}\,dt.
```

The common-horizon quantity $`C_T/T`$ compares accumulated modeled residence across realizations. If $`T_*<T`$ because the description becomes unresolved, it is not a complete-horizon physical residence fraction. The alternative $`C_T/T_*`$, defined when $`T_*>0`$, is the fraction of the modeled interval spent inside; it uses a different denominator and answers a different question. Both require the stopping status.

For species $`j`$, let $`\rho_j(\mathbf r,t)`$ be the mass density supplied consistently with Step 1. Define the physical density exposure

```math
J_{j,T}=\int_0^{T_*}\rho_j(\mathbf r_t,t)\,dt.
```

When a dimensionless density weight is useful, declare a positive reference density function before evaluation:

```math
I_{j,T}=\int_0^{T_*}
\frac{\rho_j(\mathbf r_t,t)}{\rho_{{\rm ref},j}(t)}\,dt,
\qquad \rho_{{\rm ref},j}(t)>0.
```

| Quantity | Meaning | Dimensions |
|---|---|---|
| $`C_T`$ | Time spent inside the declared region | Time |
| $`J_{j,T}`$ | Species mass density integrated along the path | Mass × inverse volume × time |
| $`I_{j,T}`$ | Density relative to a declared reference, integrated along the path | Time |
| $`C_T/T`$ | Modeled residence accumulated per common horizon | Dimensionless |
| $`I_{j,T}/T`$ | Modeled normalized exposure accumulated per common horizon | Dimensionless |

The reference density is a convention, not an inferred universal central density. Its definition accompanies $`I_{j,T}`$. The normalized density can exceed unity unless the reference independently bounds the sampled density. Consequently, $`I_{j,T}/T`$ is not generally a fraction between zero and one.

Residence and density exposure retain different information. Density can vary along a path that remains inside the region, and a trajectory outside the region can still sample nonzero density. Exposure can be restricted to the region by including the same membership indicator, but that restricted quantity must be labeled separately. None of these integrals is an encounter count.

## Preserve time, position, and velocity together

Use a time-augmented occupation measure on physical time and phase space:

```math
\mu_T(\mathcal A)=\int_0^{T_*}
\mathbf{1}_{\{(t,\mathbf r_t,\mathbf v_t)\in\mathcal A\}}\,dt,
\qquad
\mathcal A\subseteq[0,T]\times\mathbb R^3\times\mathbb R^3.
```

Its total mass is $`T_*`$ and its units are time. It retains the chronological association between motion and environment. The observables above are weighted integrals of the same measure:

```math
C_T=\int\mathbf{1}_{\{g(\mathbf r,t)<0\}}\,
\mu_T(dt,d^3r,d^3v),
\qquad
J_{j,T}=\int\rho_j(\mathbf r,t)\,
\mu_T(dt,d^3r,d^3v).
```

Removing time labels can lose the relation to an evolving host. Removing velocity can lose information required for subsequent dynamical questions. Aggregated residence or exposure values therefore accompany, rather than replace, the timed trajectory and event record. A discretized representation must declare its temporal and phase-space resolution.

## Ensemble averages and survival weighting

At fixed host parameters $`\Theta`$, fixed remnant mass $`m_{\rm r}`$, and merger selection $`\mathcal M`$, let $`p_{\rm alive}`$ be the phase-space density before the modeled endpoint:

```math
p_{\rm alive}(\mathbf r,\mathbf v,t)\,d^3r\,d^3v
=P[\mathbf r_t\in d^3r,\mathbf v_t\in d^3v,t<T_*
\mid\Theta,m_{\rm r},\mathcal M].
```

The launch law and scattering process are included in this probability. The density is not normalized to one after paths stop:

```math
\int p_{\rm alive}(\mathbf r,\mathbf v,t)\,d^3r\,d^3v
=S(t)=P[t<T_*\mid\Theta,m_{\rm r},\mathcal M].
```

Here survival means that the transport description is still being followed, not that the black hole physically ceases to exist at a stop. For nonnegative integrands, the path definitions imply

```math
\mathbb E[C_T\mid\Theta,m_{\rm r},\mathcal M]
=\int_0^Tdt\int d^3r\,d^3v\,
\mathbf{1}_{\{g(\mathbf r,t)<0\}}p_{\rm alive}(\mathbf r,\mathbf v,t),
```

```math
\mathbb E[J_{j,T}\mid\Theta,m_{\rm r},\mathcal M]
=\int_0^Tdt\int d^3r\,d^3v\,
\rho_j(\mathbf r,t)p_{\rm alive}(\mathbf r,\mathbf v,t).
```

If a density conditional on surviving trajectories is used instead, it must be multiplied by $`S(t)`$. Silent renormalization changes the observable and gives disproportionate weight to trajectories remaining in the calculation. These identities specify a future consistency check between path accumulation and a transport-density calculation; neither has been evaluated here.

Host uncertainty, launch uncertainty, and stochastic gravitational scattering remain distinct. Average conditional quantities over the joint host and launch law defined in Step 1; preserve correlations among remnant mass, kick, merger position, and host. Region definitions and reference-density conventions must also be conditioned consistently or compared explicitly. No fitted posterior, model weights, or uncertainty ranges are assigned here.

## Event status and incomplete histories

Crossing events and termination causes occupy separate fields: a path may re-enter the region and later escape. Every record must distinguish the following outcomes:

| Status | Interpretation |
|---|---|
| Entry or re-entry observed | The declared sign-changing event occurs before the modeled endpoint |
| No exit observed | An initially interior trajectory does not reach the exterior during the modeled interval; no re-entry delay is assigned |
| Horizon reached without entry or re-entry | The event time is right-censored by the finite observation horizon |
| Escape under the declared rule | Escape is a separate physical outcome within that rule and can compete with return |
| Strong encounter unresolved | Subsequent motion and exposure need an additional transition description |
| Validity interval ended | The supplied host or transport approximation no longer supports continued prediction |

Termination and crossing events with the same recorded time require a declared ordering rule. A spatial boundary crossing is only an escape event if the chosen escape rule makes it one. Escape, numerical failure, and scientific model invalidity must not be pooled into a single absence-of-return category.

After an unresolved stop, complete-horizon residence and exposure remain unknown. Accumulated nonnegative integrals describe only the modeled interval; extrapolation requires a continuation model. Censoring can depend on the orbit and environment, so subsequent statistical analysis must not assume it is independent without justification.

## Symbolic handoff

The calculation specification passes the following information from Step 1 into this step and onward:

1. **Conditioned path:** $`\{t,\mathbf r_t,\mathbf v_t,m_{\rm r}\}`$, host and launch identifiers, units, weights, and numerical provenance when an implementation exists.
2. **Environment:** the same $`\rho_j(\mathbf r,t)`$ used by the host description, region function $`g`$, its boundary-motion rule, and any positive $`\rho_{{\rm ref},j}(t)`$.
3. **Termination:** $`T`$, $`T_*`$, stopping reason, escape/validity definitions, and availability of a strong-encounter event model.
4. **Derived observables:** entry and re-entry records, $`C_T`$, $`J_{j,T}`$, $`I_{j,T}`$ where defined, denominator conventions, and the time-augmented occupation measure or a traceable representation of it.
5. **Ensemble accounting:** original realization weights, survival probabilities, missing or failed trajectories, and host/launch/scattering uncertainty labels.

These outputs preserve the motion and environmental information needed to pose subsequent questions. They do not themselves determine a population of possible companions or an interaction law.

## Planned acceptance checks — unexecuted

- **Units and bounds:** verify the stated dimensions, nonnegative exposures, $`0\le C_T\le T_*\le T`$, and total occupation-measure mass $`T_*`$. Do not impose a unit upper bound on normalized exposure without a density bound.
- **Crossing semantics:** distinguish sign changes from tangential contacts, handle inside/outside/boundary launches explicitly, and separate boundary growth from remnant motion. Verify that re-entry never resets position or velocity.
- **Stopping semantics:** distinguish horizon censoring, declared escape, unresolved strong encounters, and invalid models; verify that no accumulation or invented continuation occurs beyond $`T_*`$.
- **Temporal resolution:** demonstrate convergence of accumulated residence, exposure, and event times under a declared numerical refinement procedure. Brief passages must remain resolved.
- **Equivalent representations:** compare weighted path averages with the survival-weighted transport-density identities where both calculations are defined. Retain original ensemble weights after termination.
- **Dependence and uncertainty:** preserve joint launch/host conditioning and distinguish stochastic variation from uncertain physical inputs. Report numerical sampling uncertainty separately.
- **Region and reference conventions:** verify that changed regions or normalizations alter the declared observable consistently, and record each convention with its result.

No pass status is assigned to these checks. They make the first numerical implementation reviewable before its output is interpreted as a physical prediction.

The [outer-boundary result](outer-boundary.md) supplies a concrete validity-stop example for an isolated-host control. Exposure accumulated before that stop is a stopped-path quantity. It cannot be interpreted as the remnant's complete future exposure or as evidence that it never returns.
