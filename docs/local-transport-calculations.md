# Reconstructed host and local scattering

The calculation connects the observational inputs to the first step of the recoil problem: one host now supplies mutually consistent gravity, stellar mass densities, and local velocity distributions. The same density fields are available to the residence and exposure definitions in Step 2. No remnant trajectory or residence distribution has been evaluated.

## Defined experiment

NGC 3201 is the first diagnostic case. We use the actual joint posterior row with the largest nested-sampling posterior quadrature weight, preserving all thirteen parameters and the observation bundle's age, metallicity and stellar-evolution conventions. This deterministic rule does not identify the best fit, posterior mode, typical host or merger-selected environment. The selected row carries normalized weight 0.001764.

The [host protocol](../analysis/host_reconstruction/protocol.json) pins original input hashes and compatible historical model implementations. Their exact identity with the original fit run is unknown. A numerical reconstruction and a new observational fit are different operations.

The frozen spherical reconstruction contains 16 mass species and has total mass 178,056 solar masses, three-dimensional half-mass radius 10.325 pc, and a formal truncation radius of 981 pc. This extended truncation radius is a property of the selected distribution-function solution; it is not adopted as a Galactic tidal boundary. The central field-star velocity cutoff is 17.410 km/s, while isolated escape to infinity is 17.455 km/s. The source stellar-evolution escape-speed parameter, 34.807 km/s, describes a separate model ingredient.

## Checks and calculations

Independent integrations check mass, half-mass radius, force, relative potential, virial balance, and every species' number-to-mass normalization and velocity moments. Projected checks use the luminous tracer represented in the released kinematic profiles. Density and pressure interpolation preserves positivity and steep species tails.

The numerical host meets the declared 0.1% internal tolerances and 0.2% projected-tracer tolerance. Agreement with published projected posterior envelopes is a descriptive compatibility check. The selected realization is inside the 95.45% envelope at all 32 inspected radii for line-of-sight and radial proper motion, at 27 for tangential proper motion, and at 25 for total density. These are correlated model summaries, not independent tests or a fit probability.

The [collision protocol](../analysis/transport/protocol.json) evaluates twenty local states for a controlled 60-solar-mass test remnant. Radii are 0, 0.1, 1 and 3 half-mass radii. Velocities are zero, or radial/tangential at 0.25 and 0.75 of the local stellar cutoff. The full anisotropic number distribution supplies both drift and diffusion.

Results are per unit constant Coulomb logarithm, without selecting a physical impact-parameter range. Quadrature orders 40 and 64 agree to 5.50e-9 for drift and 1.31e-7 for diffusion in relative norms; zero-speed drift uses an absolute symmetry residual. Diffusion is positive semidefinite, and its independent trace identity holds to 9.0e-15. Analytic uniform-ball limits and rotated anisotropic distributions provide separate checks.

## Physical finding

The selected host contains an expected 4.28 black holes, totaling 61.2 solar masses. That species contributes about 27.7% of central mass density but 95.8% of the formal central diffusion trace at zero test velocity. By 0.1 half-mass radii its diffusion fraction is 0.000846. These are conditional model calculations, not observed counts or measured relaxation rates.

The result exposes why total density is an insufficient scattering description. Heavy species receive different mass weighting, and their velocity distributions also enter the collision integrals. It also identifies a concrete validity question: a few concentrated massive objects cannot automatically be treated as numerous independent weak scatterers.

## Connection to motion and exposure

The density used in Step 2 must come from this same host when its coefficients drive Step 1. Before a trajectory can be interpreted, the sparse remnant contribution needs a resolved or calibrated discrete-encounter treatment, with a consistent partition from weak scattering. The [boundary calculation](outer-boundary.md) proposes an isolated-control validity stop using an explicit Galactic-force tolerance and computes a separate circular opening scale. A physical outer host, eccentric orbit, escape rule, merger-conditioned launch and horizon still require specification. The earlier state at 3 half-mass radii lies outside the proposed 1% isolated-force domain.

The next physical comparison is whether that encounter treatment changes the occupation of dense regions relative to the smooth approximation under otherwise fixed assumptions. Convergence of the local integrals does not answer this trajectory question, and the current calculation supplies no merger probability or companion-acquisition rate.

[Reproduction instructions](../analysis/host_reconstruction/README.md) · [Recorded host](../results/host-reconstruction/summary.json) · [Local coefficient record](../results/local-transport-kernels/summary.json) · [Claim-to-evidence record](../paper/CLAIMS.md)
