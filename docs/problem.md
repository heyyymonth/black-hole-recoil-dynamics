# Problem: retention and access to a dense environment

**Researcher: Venkata Hemanth Athota**

A black-hole merger remnant can be gravitationally bound to its host while spending much of its orbit outside the dense central region. A further merger requires additional dynamics, including interactions with available companions. The present experiments isolate the orbital part of that problem.

E001/E002 ask whether a retention flag or a prescription that treats first return as permanent central residence preserves orbital exposure over a finite window. We quantify the errors of those prescriptions in a collisionless benchmark. Within this fixed one-parameter model, an exact return time together with the governing equations can instead identify the orbit; we do not claim an information-theoretic degeneracy. E003 then tests whether density exposure accurately approximates a specified close-passage count once relative velocities enter the calculation.

## Model

The host potential is

$$\Phi(r)=-\frac{GM}{\sqrt{r^2+b^2}}.$$

Use length unit $b$, time unit $t_d=\sqrt{b^3/(GM)}$, and kick fraction $\eta=v_k/\sqrt{2GM/b}$. The signed radial coordinate satisfies

$$\ddot{x}=-\frac{x}{(1+x^2)^{3/2}},\quad x(0)=0,\quad\dot{x}(0)=\sqrt{2}\eta.$$

Dots denote derivatives with respect to dimensionless time. We restrict the experiment to bound or stationary cases, $0\le\eta<1$. The specific dimensionless energy is $E=\eta^2-1$.

For a nonzero kick, $P$ denotes the interval from the initial central passage to the next central passage. Signed position has period $2P$; radius, density and core occupancy have period $P$. At exactly zero kick the object remains at the center, and no return event occurs.

We define an operational central sphere by $|x|<c$, with declared $c=r_c/b$. The Plummer scale $b$ and these operational radii are not claimed to be measured cluster core radii.

## Quantities measured

For a dimensionless window $H$,

$$C(H,c)=\int_0^H\mathbf{1}_{|x(t)|<c}\,dt,$$

$$I(H)=\int_0^H(1+x(t)^2)^{-5/2}\,dt.$$

$C/H$ is a residence fraction. $I/H$ is the mean normalized Plummer density sampled along the orbit. The latter is a geometric exposure measure. Turning it into an encounter rate would additionally require the relevant species density, velocity distribution, cross section, and binary properties. It is not a merger probability.

## Comparisons

- **Assumed central:** $I_{\rm central}(H)=H$.
- **First return, then central:** $I_{\rm wait}(H)=\max(H-P,0)$ for a nonzero kick.
- **Continuous trajectory:** $I(H)$ from the actual collisionless orbit.

For a stationary object both comparisons reduce to $H$. The second prescription is an explicitly constructed simplification, not the return-time prescription from a published friction model.

## E003: connection to potential companions

The [close-passage experiment](e003-encounter-accessibility.md) supplies a tracer-binary density, an isotropic Plummer velocity distribution, and a two-body focusing cross section. It compares density exposure with the integrated passage rate, using the same physical window and a stationary control. Its declared outer-radius screen selects test orbits without truncating the host.

The larger physical sequence remains merger, recoil, host motion, companion acquisition, binary evolution, and another merger. E001/E002 validate the host-motion calculation. E003 tests an encounter-opportunity ingredient; it finds that density alone overestimates normalized passage counts in two of its three nonzero-kick cases. Neither a close passage nor its expected count determines companion acquisition or a second merger.

All experiments omit friction, scattering that changes the orbit, binary formation, tidal evolution, and host evolution. Their results apply to their declared fixed models. E003 adds a conditional rate calculation, not self-consistent encounter dynamics or observed-event inference.

## E004: external cluster comparison

The [snapshot comparison](e004-cluster-validation.md) audits actual BH–BH targets and replaces the analytical host with a frozen spherical mass profile from an evolved CMC model. Matching mass and half-mass radius is insufficient to match ballistic return here. The two available targets cannot calibrate a continuous bath. This strengthens the need for the [specific population and live-motion comparisons](cluster-dynamics-context.md); it does not validate collisional trajectories.

## E005: breadth of the target-availability constraint

The [late-time summary audit](e005-target-availability.md) finds small instantaneous BH–BH populations across the inspected release, extending the count evidence beyond E004. Its aggregate counts cannot supply joint companion phase space. The next precise question is whether following changing target numbers and properties, conditional on remnant-producing mergers, improves companion-access predictions during displacement. Counts, population weighting and live motion must be validated separately before inferring a repeated-merger rate.
