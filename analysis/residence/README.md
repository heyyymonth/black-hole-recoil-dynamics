# Collisionless remnant residence

Twelve deterministic reference trajectories use the six adopted host potentials, with equal-energy radial and tangential launches. All species contribute to gravity and all collisions are disabled in the motion. BH encounter rates and moments are evaluated along the resulting paths as diagnostics. These are not the matched collisional residence ensembles.

[Physical specification and interpretation](../../docs/collisionless-residence.md) · [Protocol and acceptance thresholds](protocol.json) · [Results](../../results/collisionless-residence/summary.json)

Run from the repository root using the pinned numerical environment:

```sh
python analysis/residence/run.py
python analysis/residence/verify.py --record
python analysis/residence/verify.py
python -m pytest -q analysis/residence
```

Use `--cluster NGC3201` to calculate one declared host, or `--aggregate` to combine six existing records. Aggregation rejects mismatched source, protocol and implementation hashes. Read-only verification does not integrate trajectories or recompute encounter moments.

| File | Purpose |
|---|---|
| [orbits.py](orbits.py) | Cartesian gravity integration, crossing roots and stopped occupation |
| [radial_reference.py](radial_reference.py) | Independent spherical energy/angular-momentum residence integral |
| [run.py](run.py) | Registered launches, orbit checks and refined encounter budgets |
| [verify.py](verify.py) | Recorded identities, acceptance and provenance |

Each host record includes both launch directions, dimensionless position/velocity tracks, entry/exit and turning times, independent radial results, and every retained local and time quadrature refinement. The reported residence and exposure retain their horizon and stopping convention. Parameter sensitivity and physical calibration remain separate from numerical checks.
