# Mass-resolved host reconstruction

The calculation asks whether a published joint posterior sample can supply mutually consistent gravity, density and velocity distributions for recoil transport. NGC 3201 is the first controlled case; it is not a representative sample of Galactic clusters.

## Inputs and interpretation

The versioned [protocol](../analysis/host_reconstruction/protocol.json) pins the original chain and observation files by SHA-256. The sample with the largest nested-sampling posterior quadrature weight is selected without independently combining marginal parameter estimates. This row is neither a maximum-likelihood fit nor a posterior mode, and one row does not propagate host uncertainty.

Historical GCfit, LIMEPY and stellar-evolution code revisions are pinned. They are compatible reconstruction choices, not proven original-run versions: the released chain does not preserve those exact hashes. The observation bundle supplies age, metallicity and the stellar-evolution escape-speed parameter. That parameter is distinct from the present reconstructed potential. The source constructor converts mass from million solar masses, the logarithmic anisotropy radius, and the historical percentage retention parameter.

## Consistency checks

The same mass bins supply the gravitational field and local distribution functions. The source distribution function is mass normalized; dividing each bin by its stellar mass produces the number distribution used by the transport equations. Tiny artificial tracer bins are retained in the source gravity for comparability, quantified, and excluded from scattering.

The relative potential is measured from the truncation surface. With zero physical potential at infinity, the physical potential is minus the relative potential minus GM divided by the truncation radius. The stellar cutoff speed and isolated escape speed are reported separately.

Independent quadrature checks total mass, the potential integral, the force, distribution-function density and velocity moments, projected tracer moments, and the virial relation. Numerical tolerances and refinement orders are declared in the protocol before execution. Published projected model envelopes provide a descriptive compatibility comparison. They are correlated posterior summaries, not independent measurements or evidence of a new fit.

## Scope

The calculation constructs one frozen, spherical, nonrotating candidate environment. It makes no claim to reconstruct a merger epoch, fit a Galactic tide, or calibrate a collisional closure. Numerical agreement does not establish physical fidelity. Trajectories require separately specified remnant masses, launches, collision cutoffs, validity bounds and stopping rules.
