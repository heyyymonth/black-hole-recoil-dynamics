# Central residence and density exposure of recoiling remnants in a Plummer potential

**Venkata Hemanth Athota**<br/>
Independent researcher<br/>
[bharadwaj.hemanth.002@gmail.com](mailto:bharadwaj.hemanth.002@gmail.com)

Technical note BHRD-TN-001 · Version 0.1.0 · 18 October 2026<br/>
Black-hole Recoil Dynamics · Not peer reviewed

## Abstract

We compute the first central-return time, finite-window central residence, and normalized density exposure of a test particle launched radially from the center of a static, isolated Plummer potential. Energy quadrature with an endpoint-regularizing transformation supplies the production calculation; direct orbit integration and analytical limits provide independent checks. At a kick speed equal to 0.99 of central escape speed, the center-to-center period is 792.552 dynamical times, the mean normalized density sampled over that period is 0.00129861, and the fraction spent inside one Plummer scale radius is 0.00193120. A prescription that suppresses exposure until first return and assumes permanent central residence afterward underestimates early exposure and eventually overestimates it. Near escape, the period diverges while exposure per bound excursion approaches a finite limit. These results provide a reproducible collisionless benchmark for simplified orbital-exposure prescriptions. They do not infer encounter rates, repeated-merger probabilities, or constraints from gravitational-wave observations.

## 1. Physical question and context

A retained merger remnant need not remain in the dense part of its host. Published work already studies recoil displacement and frictional recentering [2], while hierarchical-merger calculations include host evolution and uncertain companion pairing [3]. We isolate a narrower calculation: how accurately do assumed central residence and a first-return-then-central prescription reproduce the density sampled by an undamped orbit?

The experiment does not test the frictional-return model of [2]. Its first ballistic center crossing is distinct from dissipative recentering. Nor does it establish that an exact return time lacks orbital information: in this fixed one-parameter family, the return time and governing equations can identify the orbit. The subject is the error made by the two explicitly defined approximations.

We use the standard Plummer potential [1] as a smooth benchmark. The gravitational-wave connection is the post-recoil motion of a possible black-hole merger remnant; gravitational radiation, binary formation, and merger dynamics are not simulated here.

<!-- pagebreak -->

## 2. Model and numerical method

With host mass $M$, scale radius $b$, and gravitational constant $G$, define

$$\Phi(r)=-\frac{GM}{\sqrt{r^2+b^2}},\qquad t_d=\sqrt{\frac{b^3}{GM}},\qquad \eta=\frac{v_k}{\sqrt{2GM/b}}.$$

Let $x$ be signed displacement in units of $b$ and let time be measured in $t_d$. The initial-value problem and conserved specific energy are

$$\ddot{x}=-\frac{x}{(1+x^2)^{3/2}},\quad x(0)=0,\quad \dot{x}(0)=\sqrt{2}\eta,\quad E=\eta^2-1.$$

We restrict attention to $0\leq\eta<1$. For a nonzero kick, $P$ is the interval between successive center crossings. Signed position has period $2P$; radius and density have period $P$. The exactly stationary case has no return event and is treated separately.

### 2.1. E001: first central return

Energy conservation yields the apocenter and return time:

$$x_a=\frac{\eta\sqrt{2-\eta^2}}{1-\eta^2},\qquad P=2\int_0^{x_a}\frac{dx}{\sqrt{2[(1+x^2)^{-1/2}-(1-\eta^2)]}}.$$

Set $x=x_a\sin\theta$, $s=\sqrt{1+x_a^2\sin^2\theta}$, and $s_a=(1-\eta^2)^{-1}$. The endpoint singularity is removed:

$$\frac{dt}{d\theta}=\sqrt{\frac{s s_a(s+s_a)}{2}},\qquad P=2\int_0^{\pi/2}\frac{dt}{d\theta}\,d\theta.$$

### 2.2. E002: residence and density exposure

For horizon $H$ and operational sphere radius $cb$, define

$$C(H,c)=\int_0^H\mathbf{1}_{|x(t)|<c}\,dt,\qquad I(H)=\int_0^H(1+x(t)^2)^{-5/2}\,dt.$$

Here $C/H$ is a residence fraction, while $I/H$ is the mean host density relative to its central value. Neither is an encounter probability. The comparisons are $I_{\rm central}(H)=H$ and $I_{\rm wait}(H)=\max(H-P,0)$ for nonzero kicks; both equal $H$ for a stationary object.

The production calculation uses energy quadrature with absolute and relative tolerances of $10^{-11}$. Exposure integrals are split at radii 0.01, 0.1, 1, 10, and 100 when these lie on the orbit. This resolves narrow central passages at large apocenter. Complete periods accumulate exactly; the partial period is recovered by a bracketed inversion of outgoing travel time. Core crossing times are found by quadrature, not by a time-sampled indicator.

<!-- pagebreak -->

## 3. First-return results

E001 tabulates nine kick fractions from 0.1 to 0.999. The plotted curve extends to 0.9999. All lengths and times below use the model units defined in Section 2.

| Kick fraction | Apocenter / b | First return / dynamical time |
|---:|---:|---:|
| 0.50 | 0.881917 | 4.395561 |
| 0.80 | 2.591534 | 11.371439 |
| 0.90 | 5.167285 | 27.972769 |
| 0.95 | 10.207544 | 74.156668 |
| 0.99 | 50.241305 | 792.551744 |
| 0.999 | 500.249126 | 24856.344701 |

![First central-return time](../results/e001/return_time.png)

*Figure 1. Energy-quadrature return times with small-kick and near-escape limits. Binding alone does not imply a short excursion. The return is a crossing, not dissipative settling.*

The small-amplitude equation is harmonic, so $P\to\pi$ as $\eta\to0^+$. At exactly zero kick, the particle instead remains stationary. Near escape,

$$P\sim\frac{\pi}{\sqrt{2}}(1-\eta^2)^{-3/2}.$$

At $\eta=0.99$, the apocenter is about $50b$ and the excursion takes about $793t_d$. Such large excursions also expose the limitations of the isolated potential: tides, external structure, or host evolution could materially alter a real remnant's trajectory. The tabulated values are conditional on the stated ideal model.

<!-- pagebreak -->

## 4. Residence and exposure results

E002 uses 11 kicks from 0 to 0.999 and spheres with $c=0.1,0.3,1$. Seven fixed horizons span 1 to 1000 dynamical times. For the ten nonstationary kicks, a second grid uses $H/P=0.1,0.5,1,2,10$. Fixed and period-scaled windows are distinct comparisons.

| Kick fraction | Period P | Exposure I(P) | Mean density I(P)/P | Fraction inside r < b |
|---:|---:|---:|---:|---:|
| 0.30 | 3.497933 | 2.774830 | 0.793277 | 1.000000 |
| 0.80 | 11.371439 | 1.369228 | 0.120409 | 0.174465 |
| 0.99 | 792.551744 | 1.029215 | 0.00129861 | 0.00193120 |

At $\eta=0.99$, the time inside $r<b$ is 1.530574 dynamical times per radial period. Thus the mean sampled density is approximately 0.130% of the central value, while the residence fraction inside this sphere is 0.193%. These percentages describe one ideal orbit, not an observed population.

![Exposure comparison](../results/e002/exposure_comparison.png)

*Figure 2. Mean density sampled up to each horizon for three kick fractions. The first-return comparator suppresses all initial exposure and then assumes unit density permanently. The undamped orbit repeatedly departs from the center.*

At $H=P$, the wait prescription gives zero, whereas the orbit has accumulated $I(P)>0$. At $H=10P$ for $\eta=0.99$, the orbit gives 10.292151 dynamical times of exposure and the wait prescription gives 7132.965700. The error changes sign. More generally, as $H$ grows, the wait prescription approaches unit mean density while the orbit approaches $I(P)/P<1$ for every nonzero kick. It is therefore not uniformly conservative.

An orbit stays within a sphere when its apocenter does not exceed its radius:

$$\eta\leq\sqrt{1-\frac{1}{\sqrt{1+c^2}}}.$$

For $c=1$, the threshold is 0.541196. Equality touches the sphere boundary instantaneously and does not change an integrated residence fraction. Even an orbit entirely within the sphere samples a varying density; residence and density exposure are different quantities.

<!-- pagebreak -->

## 5. Analytical and numerical verification

For a stationary object, $I(H)=C(H,c)=H$. For nonzero bound kicks, both exposure and occupancy repeat every $P$, hence $I(nP)=nI(P)$. In the small-kick limit, $I(P)\to\pi$ and $I(P)/P\to1$.

Near escape, the marginally bound outgoing trajectory obeys $\dot{x}=\sqrt{2}(1+x^2)^{-1/4}$. The limiting exposure of a full bound excursion is

$$\lim_{\eta\to1^-}I(P)=\sqrt{2}\int_0^\infty(1+x^2)^{-9/4}dx=\sqrt{\frac{\pi}{2}}\frac{\Gamma(7/4)}{\Gamma(9/4)}.$$

Its value is 1.0166557018. Combining this with the period asymptote gives

$$\frac{I(P)}{P}\sim\frac{\Gamma(7/4)}{\sqrt{\pi}\Gamma(9/4)}(1-\eta^2)^{3/2}.$$

This is a limit from bound orbits. An exactly marginally bound particle does not return in finite time.

The recorded verification contains 27 passing tests. Independent checks integrate the signed equations of motion using DOP853 [4], compare central crossings, monitor energy, and integrate accumulated density as an additional state variable. Core residence is checked with event-segmented trajectories. Other tests cover stationary behavior, sphere thresholds, period accumulation, crossing continuity, density as the derivative of exposure, invalid inputs, and convergence near escape. The strongest near-escape numerical cases reach $\eta=0.9999$; no claim is made of uniform numerical accuracy arbitrarily close to escape.

Separately, the result verifier reproduces 146 recorded rows: nine return-time rows, ten cycle rows, 77 fixed-window rows, and 50 period-scaled rows. It also checks the recorded input and output hashes. This second step verifies stored results against their implementation; it is not an additional independent physics calculation. Numerical tolerances are solver controls, not error bars on the physical model.

## 6. Interpretation and limits

The two tested prescriptions fail to preserve continuous exposure for distinct reasons: one assumes uninterrupted central residence, and the other removes the initial outbound passage before imposing permanent central residence. Their errors are demonstrated within the collisionless model.

Density exposure does not supply a companion distribution, relative-velocity distribution, encounter cross section, binary fraction, or companion-acquisition mechanism. The potential is static and isolated, and the test particle does not experience friction, scattering, or tidal removal. Consequently, these results do not estimate a repeated-merger probability, distinguish astrophysical formation channels, or replace evolving cluster calculations [3]. Their contribution is an inspectable benchmark and approximation test; publication-level novelty is not established.

<!-- pagebreak -->

## 7. Reproducibility and availability

All calculations, tables, figures, and this note are distributed in the version 0.1.0 research release. Python 3.11 is the supported reference environment. Numerical dependencies and transitive packages are pinned in `requirements-lock.txt`. The paper has separate pinned build requirements in `paper/requirements.txt`.

`experiments/config.json` declares the grids and tolerances. `results/run_all.json` records the source commit, source and output SHA-256 hashes, runtime versions, and run time. The manifest's source commit identifies the code that generated the results; the release also contains the subsequent results commit. `output/pdf/build.json` records the manuscript, figures, and build inputs used for this paper.

From the repository root, install the numerical dependencies, run the independent tests, and verify the committed result tables. The documented runner regenerates all tables and plots. These operations need no network access after installation. The paper source and build script are included, and the PDF is distributed with the release.

The research artifacts are provided under the MIT license. Prior publications and third-party dependencies retain their respective rights. Citation is requested for scholarly reuse, separately from the license conditions. The release is a versioned GitHub snapshot; it has no DOI and the note has not been peer reviewed.

**Research release:** [github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.0](https://github.com/heyyymonth/black-hole-recoil-dynamics/releases/tag/v0.1.0)

**Suggested citation:** Athota, Venkata Hemanth (2026). *Black-hole recoil dynamics: central residence and density exposure in a Plummer potential*. Version 0.1.0. Software, numerical results, and technical note BHRD-TN-001. Machine-readable citation and BibTeX are included in the repository.

## References

[1] H. C. Plummer (1911). On the Problem of Distribution in Globular Star Clusters. *Monthly Notices of the Royal Astronomical Society*, 71, 460-470. [doi:10.1093/mnras/71.5.460](https://doi.org/10.1093/mnras/71.5.460).

[2] T. Islam (2026). Recoil kicks from binary black hole mergers in GWTC catalogs: implications for retention and hierarchical mergers. [arXiv:2604.04546v2](https://arxiv.org/abs/2604.04546v2). Selected displacement and recentering sections provide context; their model is not implemented here.

[3] S. Torniamenti et al. (2024). Hierarchical binary black hole mergers in globular clusters: mass function and evolution with redshift. *Astronomy & Astrophysics*, 688, A148. [doi:10.1051/0004-6361/202449272](https://doi.org/10.1051/0004-6361/202449272); [arXiv:2401.14837v2](https://arxiv.org/abs/2401.14837v2).

[4] SciPy developers. `scipy.integrate.solve_ivp` documentation. [Official method documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html). The experiments use SciPy 1.10.1; runtime versions are recorded in the result manifest.
