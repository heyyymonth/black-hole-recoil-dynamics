# E003: density exposure and velocity-weighted close passages

**Researcher:** Venkata Hemanth Athota  
**Data:** four deterministic model calculations; no observational data.  
**Protocol:** [question, controls, and decision rule](e003-protocol.md), committed as [4a729da](https://github.com/heyyymonth/black-hole-recoil-dynamics/commit/4a729daf9c392c0d7cec8b825b69db70d153af16) before execution.  
**Generating source:** [4161b85](https://github.com/heyyymonth/black-hole-recoil-dynamics/commit/4161b85d006bd2b5de76f4ee598393d34940eb02). The [manifest](../results/run_e003.json) records exact source and output hashes.

## Physical question

Can time integrated through a dense region predict a recoiling remnant's opportunity to pass near another binary, or does its velocity materially change the answer?

E001 established the ballistic orbit; E002 measured its density exposure. E003 supplies a specified target population and passage cross section, testing the next physical connection. It finds that the density proxy overestimates the normalized passage count by **53.1% and 98.6%** at kick fractions 0.8 and 0.9 in this configuration. Those cases fail the predeclared 10% approximation tolerance. The 0.3 case differs by 4.4% and passes.

This is a conditional approximation result, not a hierarchical-merger fraction. The original assembly question additionally requires an encounter outcome, a new binary, its subsequent evolution, and a population model.

## Host, target population, and common window

| Input | Declared value |
|---|---:|
| Fixed Plummer host mass M | 1,000,000 solar masses |
| Plummer scale b | 1 pc |
| Remnant mass | 60 solar masses |
| Tracer binaries | 100, each of total mass 60 solar masses |
| Binary semimajor axis | 0.1 AU |
| Close-passage radius R | 1 AU from the binary center of mass |
| Kick fractions eta | 0, 0.3, 0.8, 0.9 |
| Common horizon H | 10 dynamical times = 0.149096 Myr |
| Central escape speed | 92.746 km/s |
| Apocenter admission limit | 5 half-mass radii = 6.523830 b |

The target binaries are imposed, unsegregated Plummer tracers with isotropic equilibrium velocities. Their mass is 0.6% of the fixed host mass. This is not a model of an evolved black-hole subsystem. The 1 AU boundary is ten binary semimajor axes, supporting a monopole approximation for entry into the binary's neighborhood. Crossing this boundary does not establish a strong interaction or exchange.

The apocenters are 0, 0.455613, 2.591534, and 5.167285 b, so all four cases pass the declared screen. This screen is motivated by the explicit choice in [Islam, section V.1.3](https://arxiv.org/html/2604.04546v2); it does not truncate the potential or tracer distribution. In contrast, the E002 eta=0.99 orbit reaches 50.24 b and would fail this particular screen. Its isolated-model result remains valid within its original scope.

All comparisons use the same physical time. At H=10, eta=0.3 has completed approximately 2.86 radial periods, while eta=0.8 and 0.9 have not completed their first return. No case is given a longer observing window to complete an orbit.

## Rate calculation

Use velocity unit V = sqrt(GM/b), time unit td = b/V, x = r/b, and dimensionless time t. Let s = sqrt(1+x^2), so the local escape speed in these units is sqrt(2/s). The target number density is

$$
n_{\rm bin}(r)=n_0s^{-5},\qquad n_0=3N_{\rm bin}/(4\pi b^3).
$$

The isotropic Plummer distribution has f(E) proportional to (-E)^(7/2), as documented in the [versioned galpy reference](https://docs.galpy.org/en/v1.12.0/reference/dfplummer.html). At fixed radius, its normalized target-speed distribution is

$$
p(y)=\frac{512}{7\pi}y^2(1-y^2)^{7/2},\quad y=u/v_{\rm esc}(r),\quad 0\le y\le1.
$$

For remnant speed v and target speed u, the isotropic angular averages are

$$
\langle w\rangle_\Omega=\max(v,u)+\frac{\min(v,u)^2}{3\max(v,u)},\qquad
\langle w^{-1}\rangle_\Omega=\frac{1}{\max(v,u)}.
$$

The zero-speed limit is integrated explicitly. The single point v=u=0 does not cause a divergent velocity-integrated rate. Remaining speed integrals are evaluated using incomplete beta functions, with independent numerical speed and angle integration as a check.

Conservation of energy and angular momentum for an isolated two-body approach gives impact parameter squared R^2 + 2G(m_rem+m_bin)R/w^2. Thus

$$
\Gamma(r,v)=n_{\rm bin}(r)\pi R^2\left[\langle w\rangle+
\frac{2G(m_{\rm rem}+m_{\rm bin})}{R}\langle w^{-1}\rangle\right].
$$

This local kinetic-theory rate neglects host tides during the passage and the binary's internal multipoles. Its integral is the expected passage count along the **unperturbed** trajectory in the specified dilute bath. We call it optical depth tau, without interpreting it as a merger probability.

The two dimensionless orbit integrals are

$$
J_g=\int_0^H s^{-5}\langle w/V\rangle\,dt,\qquad
J_f=\int_0^H s^{-5}\langle V/w\rangle\,dt.
$$

With lambda = 2(m_rem+m_bin)b/(MR),

$$
\tau=n_0\pi R^2b\,[J_g+\lambda J_f].
$$

The factor b follows from V td = b; time conversion is applied once. Here lambda = 49.5035535. The orbital integration uses a turning-point-regularized energy quadrature and exact cycle accounting, with absolute and relative tolerances 1e-10.

The stationary central control sets the normalization: Q = tau(eta,H)/tau(0,H). E002 supplies the density proxy D = I(H)/H. Its signed relative error is D/Q - 1. The 10% criterion was chosen as a benchmark screen, not inferred from data or justified as a GW population-inference requirement.

## Executed results

| Kick fraction | Density proxy D | Passage ratio Q | Proxy error | Expected passages tau | 10% screen |
|---|---:|---:|---:|---:|---|
| 0, stationary | 1.000000 | 1.000000 | 0 | 1.608016e-6 | Control |
| 0.3 | 0.784719 | 0.751564 | +4.41% | 1.208527e-6 | Pass |
| 0.8 | 0.076256 | 0.049820 | +53.06% | 8.011086e-8 | Fail |
| 0.9 | 0.058113 | 0.029268 | +98.55% | 4.706416e-8 | Fail |

For example, the 0.9 orbit samples 5.81% of the control's density exposure but only 2.93% of its close-passage count. The density-only estimate is therefore approximately 1.99 times the velocity-weighted result. These are fractions of the declared stationary control, not percentages of remnants that merge.

![Density proxy and velocity-weighted passages at a common time window](../results/e003/proxy_comparison.png)

Figure: four selected cases at H=10 td. Positive error means the density proxy overestimates passage opportunities. Bars are deterministic outputs; they are not observations, Monte Carlo counts, or confidence intervals. [Vector PDF](../results/e003/proxy_comparison.pdf) · [SVG](../results/e003/proxy_comparison.svg) · [Full-precision table](../results/e003/encounters.csv).

## Why the proxy fails in two cases

| Kick fraction | Density ratio | Geometric contribution / its control | Focusing contribution / its control | Focusing share of total passages |
|---|---:|---:|---:|---:|
| 0.3 | 0.784719 | 0.831705 | 0.750975 | 99.19% |
| 0.8 | 0.076256 | 0.122561 | 0.049285 | 98.21% |
| 0.9 | 0.058113 | 0.111016 | 0.028668 | 97.23% |

The two velocity effects act in opposite directions. Motion increases the geometric flux relative to a density-only estimate, but faster passage through the dense center reduces the inverse-speed weighting of gravitational focusing. Focusing supplies over 97% of the calculated encounters in every case, so the latter effect controls the total. Slow motion near apocenter does not compensate for the low tracer density there.

The ordering of the three nonzero kicks is unchanged. E003 rejects quantitative accuracy of the density proxy in two specified cases; it does not find a reversal in ordering. Absolute counts are very small: even the stationary control has only 1.61e-6 expected passages in this window and imposed population. A large relative error here does not establish a large astrophysical merger-rate effect.

## Validation and stopping decision

The complete suite passes **62 tests**, including 35 added for E003. The new checks include direct numerical angle and speed integrals; analytical stationary and fast-motion limits; independent DOP853 trajectories with separately integrated numerical velocity moments for all four cases at H=10; period accounting; tighter tolerances; physical units; count, radius, and mass scaling; and rejection of the eta=0.99 case by the admission screen.

Independent trajectory comparisons pass at 1e-8 relative tolerance, below the predeclared 1e-6 validation requirement. These are acceptance tolerances, not an estimate of astrophysical model uncertainty. The result verifier recomputes all 19 fields in each of four rows, physical scales, and the decision summary, and checks every recorded hash. The existing 146 baseline rows and four-page E001/E002 paper still pass their provenance checks.

The declared experiment is complete: velocity weighting is necessary for quantitative passage counts in the two higher-kick examples. Repeating identical runs or increasing kick-grid density would not test friction, mass segregation, encounter outcomes, or evolving targets. Before using this ingredient to infer repeated mergers, the bath and trajectories require physical calibration and the passage-to-binary-outcome step remains to be supplied.

## Scope and limitations

The calculation assumes a fixed potential and an equilibrium, unsegregated target population. It omits dynamical friction, diffusion, collective core response, binary encounters that change the trajectory, target depletion, and host evolution. The short horizon does not prove those omissions are negligible in a real cluster. The finite optical depth does not validate the infinite-duration motion or a Poisson independence assumption.

These ingredients have substantial prior literature. [Leigh et al.](https://arxiv.org/abs/1711.10494) already connects encounters and hardening toward merger; [Islam](https://arxiv.org/abs/2604.04546v2) already treats displacement and recentering. This result tests our density approximation and is not a reproduction or falsification of their merger prescriptions. No observational formation-channel discrimination or publication-level novelty is established.

## Output conventions and reproduction

The [configuration](../experiments/e003.json) pins all parameters, illustrative unit-conversion constants, and tolerances. In [encounters.csv](../results/e003/encounters.csv):

| Columns | Meaning and units |
|---|---|
| `eta`, `kick_kms`, `apocenter_over_b` | Dimensionless kick fraction, kick in km/s, maximum radius in b |
| `horizon_over_tdyn`, `horizon_myr` | Common elapsed time in td and Myr |
| `density_ratio` | D = normalized density exposure I(H)/H |
| `geometric_integral`, `focusing_integral` | Dimensionless J_g and J_f defined above |
| `geometric_ratio`, `focusing_ratio` | Each integral divided by its own stationary value |
| `tau_geometric`, `tau_focusing`, `tau_total` | Dimensionless expected counts; total is the sum |
| `tau_stationary`, `encounter_ratio` | Control count and Q = total/control |
| `density_proxy_tau` | Approximate count D times the control count |
| `proxy_relative_error` | Signed D/Q - 1, stored as a fraction |
| `focus_fraction` | tau_focusing/tau_total |
| `exceeds_tolerance` | 1 if absolute proxy relative error exceeds 0.10; otherwise 0 |

From the repository root with Python 3.11 and the pinned dependencies installed:

```sh
python -m pytest -q
python experiments/verify_results.py
python experiments/verify_e003.py
```

To regenerate E003 after committing its source inputs:

```sh
python experiments/run_e003.py
python experiments/verify_e003.py
```

The runner refuses uncommitted source inputs, records their generating commit, and regenerates only E003. It performs no stochastic sampling. Reproduction requires no network access after dependency installation. E003 is a separately documented addition on `main`; the v0.1.1 PDF and citation metadata describe E001/E002. Cite an immutable repository commit containing these results when using E003, rather than attributing them to v0.1.1.
