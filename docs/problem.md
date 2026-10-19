# Problem: retention and access to a dense environment

**Researcher: Hemanth AV**

A black-hole merger remnant can be gravitationally bound to its host while spending much of its orbit outside the dense central region. A further merger requires additional dynamics, including interactions with available companions. The present experiments isolate the orbital part of that problem.

The question is whether a binary retention label or a first central-return time adequately represents central residence over a finite observation window. Both descriptors discard the path between the kick and later encounters. We measure that lost information directly in a collisionless benchmark.

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

The experiments omit friction, scattering, binary formation, tidal removal, and host evolution. Their results determine the behavior of these comparisons in this specified model only.
