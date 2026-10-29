# Finite-angle encounter calculations

This calculation supplies a common local encounter law, its first two velocity-change moments, and conditional event samples for the six adopted reference hosts. It provides the local inputs for the residence comparison in [the hypothesis](../../docs/transport-hypothesis.md). It does not integrate remnant trajectories or evaluate residence, exposure, or the response predictor.

## Inputs and scope

[The protocol](protocol.json) uses the existing verified joint reference rows for NGC 0104, NGC 1851, NGC 3201, NGC 5139, NGC 6121 and NGC 6397. The host meshes remain in `results/cluster-comparison/<cluster>/reference/host.npz`. The diagnostic remnant mass is 60 solar masses. Only BH species are included in this local comparison; a complete trajectory calculation must explicitly supply the common treatment of the other species.

Each host supplies its original anisotropic, number-normalized BH distribution functions. No Gaussian velocity distribution is substituted. Define the density radius by half the central total BH mass density and the local scale by the smaller of that radius and the logarithmic density-gradient scale. The tested upper impact parameters are 0.05, 0.1 and 0.2 times this local scale. These are declared conditional domains, not fitted Coulomb cutoffs.

The lower impact parameter follows a closest-approach floor that excludes the supplied strong-radiation domain from elastic scattering. The exact prescription and primary reference are given in [the calculation description](../../docs/finite-angle-encounters.md). The excluded incoming flux is recorded separately; no excluded event is continued or assigned an outcome.

## Recorded calculation

The grid contains 27 states per host: three radii, a zero-speed state and two moving orientations at each radius, and three upper impact parameters. This gives 162 local states across the six hosts. Each state records the rate, drift vector, diffusion tensor, cutoff diagnostics and quadrature-refinement history.

At one registered state per host, 100,000 independent conditional elastic increments are drawn and compared with the moments of the same numerical velocity distribution. This checks the sampler and moment calculation together. It supplies neither a chronological sequence of encounters nor a physical remnant trajectory.

Numerical acceptance is stored per state. [The aggregate summary](../../results/finite-angle-encounters/summary.json) records the accepted count and whether the complete grid meets its checks. An unresolved state retains its result and refinement history. Source, protocol and implementation hashes prevent mixing results from different calculations.

## Reproduction

From the repository root, with the dependencies in [requirements.txt](../../requirements.txt):

```bash
python analysis/encounters/run.py
python analysis/encounters/verify.py --record
python analysis/encounters/verify.py
python -m pytest -q analysis/encounters
```

To recompute one host, pass `--cluster NGC3201` to `run.py`. An aggregate can be produced only when all six records match the current inputs and implementation. The verifier requires every recorded state to meet the registered acceptance rules before it records provenance.

| File | Purpose |
|---|---|
| [physics.py](physics.py) | Elastic velocity increments, analytic impact-parameter moments, positive velocity quadrature and conditional event sampling |
| [run.py](run.py) | Adopted encounter domain, six-host calculations and sampler comparison |
| [protocol.json](protocol.json) | Physical controls, numerical refinement and interpretation |
| [verify.py](verify.py) | Recorded-result identities, acceptance checks and provenance |
| [response.py](response.py) | Generic single-transition replacement estimator and exact finite-state response |
| [test_physics.py](test_physics.py) | Conservation, quadrature, rotation and sampler checks |
| [test_response.py](test_response.py) | Independent finite-state derivative checks and isolated random streams |

## Response estimator status

For an N-step transition mixture, the estimator replaces one uniformly selected diffusion step with a jump step and resumes diffusion afterwards. N times the paired observable difference estimates the derivative at zero mixture weight. Its finite-state verification is mathematical validation of the algorithm, not a residence result for a cluster. Time-step refinement is required to establish the continuous-generator response. Extrapolating that response to a full jump trajectory remains the hypothesis to test.

Callbacks receive one path at a time with independent random streams for each path and step; randomness is shared only within a baseline/replacement pair. State-dependent draw counts cannot couple different paths. Reported sampling errors exclude discretization and physical model uncertainty.

[Local results and physical interpretation](../../docs/finite-angle-encounters.md) · [Six-cluster backgrounds](../../docs/six-cluster-comparison.md)
