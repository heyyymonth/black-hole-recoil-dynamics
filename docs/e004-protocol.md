# E004 protocol: comparison with an evolved CMC cluster snapshot

**Researcher:** Venkata Hemanth Athota. Recorded before calculating E004 population statistics or trajectory comparisons, 19 October 2026.

## Question

Are E003's unsegregated Plummer binary population and Plummer host-potential assumptions supported by an independently evolved cluster model? What part of trajectory validation can a CMC snapshot actually support?

## External reference and selection

Use the single full evolved snapshot distributed with the author-maintained [cmctoolkit example](https://github.com/NicholasRui/cmctoolkit/tree/a4876e7706f7b23a14c63a5be7df63a07f2dccd9/N8-RV1-RG8-Z0.1), rather than selecting a cluster after viewing discrepancies. The [configuration](../experiments/e004.json) pins its revision, file names, and SHA-256 checksums. This is one CMC simulation state, not observations or a representative cluster sample. Raw upstream files are cached locally, not redistributed or relicensed as MIT. Unit files are parsed as text and never executed.

CMC's Monte Carlo treatment of relaxation and orbital-phase sampling is described by [Rodriguez et al.](https://arxiv.org/abs/2106.02643). Snapshot radii are not phase-resolved trajectories. This dataset can test population assumptions and a frozen spherical potential; it cannot independently validate recoil damping, diffusion, or return-to-core trajectories in a live cluster.

## Declared comparisons

1. Check column meanings, physical units, finite positive masses/radii, binary component-mass accounting, and snapshot age. Mass columns explicitly labeled solar masses are already physical. Velocity uses the N-body time conversion, not the relaxation-time conversion used for age.
2. Record separate inventories for all systems, all binaries, single black holes, black-hole-bearing binaries, and binary black holes (both stellar-type codes equal 14). Do not use all binaries or all black holes as substitutes for binary black holes.
3. Match a comparison Plummer model to the snapshot's total system mass and three-dimensional half-mass radius. Compare binary-black-hole radii with both the empirical host mass distribution and the matched Plummer distribution. Summarize masses, separations, and local speed fractions. Report descriptive CDF distances without treating correlated particles as independent observational measurements.
4. If fewer than 30 binary black holes are available, declare continuous companion-distribution calibration unsupported by this snapshot. Keep the small sample visible; do not fit a flexible distribution or inflate it with repeated draws. The threshold is an operational screening rule, not a statistical sufficiency theorem.
5. Construct the frozen spherical-shell potential from system masses and radii. Compare central kicks at 0.3, 0.8, and 0.9 times the matched Plummer escape speed. Use identical physical kicks in the two potentials, recording either bound returns or escape. Do not compare different kicks by renormalizing each potential separately.
6. Compare apocenter and first ballistic return. A return-time difference greater than 10%, measured relative to the shell result, flags a potential-shape discrepancy for this benchmark. This is not a bound on live-cluster model error.

Compress the radial mass profile into 256, 512, 1024, and 2048 bins of approximately equal mass. Each shell's radius is its mass-weighted harmonic radius, preserving total mass and central potential. Require the final two resolutions to agree to 1% in return time and apocenter. These resolutions test numerical convergence, not a larger astrophysical sample. Validate shell physics against direct sums, a one-shell analytical solution, and smooth-Plummer recovery.

## Interpretation and stop rule

Complete the stated comparisons and report support, discrepancy, or insufficient information for each assumption separately. Do not call a frozen-potential comparison a validation of friction. The next experiment must address a diagnosed missing ingredient, with its required data and independent comparison specified. A denser kick grid is not a substitute for that ingredient.
