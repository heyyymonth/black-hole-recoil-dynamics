# Cluster dynamics: from E004 to the E005 findings

**Researcher:** Venkata Hemanth Athota. Updated 23 October 2026. E005's target-count audit is executed; joint companion calibration and live-recoil comparisons remain unexecuted.

The immediate question is: **How accurately can a reduced model predict a retained remnant's access to actual binary-black-hole targets while its orbit evolves after recoil?** This isolates one step toward repeated mergers. A close passage still requires an exchange or capture, subsequent binary evolution, and eventual merger before it contributes to hierarchical growth.

## Theory before further calculations

The [proposed recoil-transport formulation](recoil-transport-theory.md) now specifies the missing inputs: consistent mass and velocity distributions, merger position and progenitor motion, recoil direction, drift and diffusion, external acceleration, and finite-region entry and residence. It defines the conditional target kernel that connects those paths to the existing availability hypothesis. No host family, cluster, parameter values, or new numerical experiment is adopted by that formulation.

The background scatterers that drive friction are distinct from the BBH systems counted as potential companions. Observational host constraints and merger-linked target histories address different parts of the problem. They must describe compatible environments and times before they can be combined. The experiment descriptions below remain proposed validation requirements, not instructions to execute calculations in this revision.

## What existing research already supplies

| Established work | Consequence for this project |
|---|---|
| [Breen & Heggie (2013)](https://arxiv.org/abs/1304.3401) studies mass-segregated black-hole subsystems coupled to the energy demands of their host clusters. | An unsegregated black-hole bath is a benchmark assumption; concentration of black holes is not a new discovery. |
| [Kremer et al. (2020)](https://arxiv.org/abs/1911.00018) provides a catalog of 148 CMC cluster models spanning structural properties and evolutionary histories. | Calibration must vary host structure and age. E004 has analyzed one snapshot, not this catalog. |
| [Rodriguez et al. (2022), sections II.2 and II.5](https://arxiv.org/html/2106.02643) describes encounter sampling and randomly sampled orbital phases in CMC. | Snapshots support population and potential comparisons; connecting successive snapshot positions does not recover a physical recoil trajectory. |
| [Atallah et al. (2023), sections 2.4 and 2.8](https://arxiv.org/html/2211.09670v4) combines few-body interactions, gravitational focusing, and a friction delay in a fixed nuclear-cluster background. | Those ingredients already exist. A useful contribution would test the accuracy of encounter accessibility during displacement; nuclear-cluster parameters cannot simply be transferred to globular clusters. |
| [Islam (2026), section V.1.4](https://arxiv.org/html/2604.04546v2) combines retention and recentering in a phenomenological repeated-merger prescription. | Retention alone is already recognized as insufficient. Our ballistic calculations neither reproduce nor falsify this prescription. |

## What our calculations establish

[E003](e003-encounter-accessibility.md) shows that velocity weighting changes passage counts substantially in its imposed bath. [E004](e004-cluster-validation.md) finds only **two binary black holes** in the selected evolved snapshot, insufficient to calibrate a continuous companion distribution. Their central concentration motivates replacing the unsegregated assumption, without establishing a universal distribution.

For identical physical kicks, matching total mass and half-mass radius does not reproduce this snapshot's central potential. The matched Plummer model is shallower and predicts ballistic return times **2.96, 2.81, and 4.45 times** the shell-potential values. These are potential-shape comparisons, not measured damping times. [Executed tables and convergence checks](../results/e004/summary.json) do not validate live recoil trajectories.

## Next experiment: calibrate the actual target population

The prerequisite [E005 audit](e005-target-availability.md) now tests whether E004's target scarcity is confined to its one snapshot. It inspects 7,537 published summaries across 126 available standard model histories. Among 119 histories with a state near 12 Gyr, the median is two BH–BH binaries, 38 states have none, and the maximum is nine. No inspected saved state reaches the count-30 screen; the overall maximum is 15. Seven histories without an age-matched state remain missing, not zero. A separate scalar census and the E004 full-snapshot cross-check reproduce the result.

This supports treating target availability explicitly, but it does not invalidate a continuous ensemble description. PARAM files have no joint positions, velocities, internal binary orbits or object identities. Equal-weight old-cluster states also differ from the environments sampled by actual remnant-producing mergers. A census mean must not replace E003's 100 prescribed targets as though this calibrated the bath.

The sharper next question is: **conditional on a merger that leaves a retained remnant, does following the changing number and properties of available BH–BH targets improve predicted companion access during displacement?** E003 establishes the need for velocity weighting; E004 identifies population and potential mismatches; E005 establishes the breadth of small target counts in this late-time release. These are complementary tests, not factors that can be multiplied into a corrected merger rate.

Obtain versioned snapshots and object histories across declared cluster models and ages, anchored to merger events and including younger epochs where required. Measure the joint distribution of binary centers' radius, radial and tangential velocity, component masses, semimajor axis, and eccentricity, retaining the binary count, formation and destruction. Keep single black holes and other binaries separate. Preserve snapshot time weights and object identities; repeated appearances of one binary are not independent systems. Target-count evolution and merger conditioning cannot be reconstructed from the current aggregate summaries.

Compare the unsegregated stationary baseline with a parsimonious model that accounts for segregation and changing target availability. First supply identical independently determined remnant paths to both, isolating the population approximation. Hold out entire cluster histories, rather than randomly splitting rows. Require adequate distinct systems and held-out predictive stability before fitting continuous distributions; passing E004's count screen alone is insufficient. Sparse states may still support an ensemble model, but not a fitted single-state distribution from two objects. Test predicted passage integrals as well as radial profiles, because correlations with speed and separation determine accessibility. Record failed predictions and sparse regimes. This tests the bath, not recoil friction.

## Next experiment: test live recoil motion and accessibility

Use a validated collisional integrator with resolved close encounters. Pair kicked and un-kicked runs from the same microstate; vary independent microstates and kick direction. Declare the core boundary, common observation window, and encounter distance before execution. Define return as inward core crossing after exit, not exact arrival at the origin.

Compare live paths with frozen-potential and central-reset prescriptions using first-return distributions, core residence, cumulative predicted passages, and recorded encounter counts. Distinguish escape from finite-window non-return. Check integration accuracy and encounter regularization; any particle-number scaling must preserve or account for relaxation times. Retain the operational 10% proxy screen alongside ensemble uncertainty, without calling it an observational requirement.

The population calibration and live-motion tests still require additional data and integrations. E005 has completed the available-summary audit, with [protocol, results and validation](e005-target-availability.md). The broader catalog's full-data portal was not retrievable with verified HTTPS during this check; the pinned summary archive was. More permutations of the original static calculation cannot supply the missing phase space or trajectories. Novelty and gravitational-wave population implications remain unestablished.
