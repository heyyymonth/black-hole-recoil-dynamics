# NGC 3201 outer-boundary calculation

The [registered protocol](protocol.json) applies a specified Galactic field to the existing mass profile. It computes the accuracy domain of neglecting that field and a conditional circular tidal-opening scale. No Galactic orbit, remnant trajectory, stripping evolution, or residence distribution is simulated.

Use the Python 3.11 environment in [host reconstruction](../host_reconstruction/README.md), or install NumPy 1.26.4, SciPy 1.12.0 and pytest 8.2.2 for these numerical calculations and checks:

```sh
python analysis/boundary/run.py
python -m pytest -q analysis/boundary/test_boundary.py
python analysis/boundary/verify.py
```

The final command checks recorded hashes and invariants without rerunning physics. `--record` is only for explicitly pinning a reviewed result revision. The original host and collision manifest is unchanged. The [source register](sources.json) distinguishes catalog values, primary-literature comparisons, and chosen sensitivity settings.

The enclosed mass is interpolated through its finite central quantity M(<r)/r^3, with the analytic density limit at zero. The potential uses the existing reconstructed potential and its isolated normalization. Independent density integration, radial-mesh thinning, point-mass and Plummer limits, finite-difference Galactic forces, exact angular force bounds, and exact circular saddles check the calculation.

[Physical interpretation and equations](../../docs/outer-boundary.md) · [Results](../../results/outer-boundary/summary.json)
