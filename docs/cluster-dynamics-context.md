# Cluster dynamics: the question after E004

**Researcher:** Venkata Hemanth Athota. Literature checked 20 October 2026. The two experiments below are proposed, not executed.

The immediate question is: **How accurately can a reduced model predict a retained remnant's access to actual binary-black-hole targets while its orbit evolves after recoil?** This isolates one step toward repeated mergers. A close passage still requires an exchange or capture, subsequent binary evolution, and eventual merger before it contributes to hierarchical growth.

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

Obtain versioned snapshots and object histories across declared cluster models and ages. Measure the joint distribution of binary centers' radius, radial and tangential velocity, component masses, semimajor axis, and eccentricity, retaining the binary count. Keep single black holes and other binaries separate. Preserve snapshot time weights and object identities; repeated appearances of one binary are not independent systems.

Compare the unsegregated baseline with a parsimonious mass-segregated model. Hold out entire cluster histories, rather than randomly splitting rows. Require adequate distinct systems and held-out predictive stability before fitting continuous distributions; passing E004's count screen alone is insufficient. Test predicted passage integrals as well as radial profiles, because correlations with speed and separation determine accessibility. Record failed predictions and sparse regimes. This tests the bath, not recoil friction.

## Next experiment: test live recoil motion and accessibility

Use a validated collisional integrator with resolved close encounters. Pair kicked and un-kicked runs from the same microstate; vary independent microstates and kick direction. Declare the core boundary, common observation window, and encounter distance before execution. Define return as inward core crossing after exit, not exact arrival at the origin.

Compare live paths with frozen-potential and central-reset prescriptions using first-return distributions, core residence, cumulative predicted passages, and recorded encounter counts. Distinguish escape from finite-window non-return. Check integration accuracy and encounter regularization; any particle-number scaling must preserve or account for relaxation times. Retain the operational 10% proxy screen alongside ensemble uncertainty, without calling it an observational requirement.

These tests require additional simulation data and live runs. They target missing physics identified by E004; more permutations of the original static calculation cannot supply it. Novelty and gravitational-wave population implications remain unestablished.
