# Galactic orbit and host-support calculation

[Protocol](protocol.json) · [Source inputs](../../data/orbit/README.md) · [Results](../../results/orbit/summary.json)

This stage integrates the **cluster's Galactic reference orbit**, not a kicked remnant. It applies the resulting time-dependent field to the existing host as a diagnostic. It does not trim the host, refit its distribution function, simulate stripping, or claim a tidal equilibrium.

With Python 3.11, NumPy 1.26.4, SciPy 1.12.0 and pytest 8.2.2:

```sh
python analysis/orbit/run.py
python -m pytest -q analysis/orbit/test_galaxy.py
python analysis/orbit/verify.py
```

The verifier checks recorded source extraction, hashes and invariants without integrating. `--record` explicitly pins a reviewed result revision. Results use the existing host's gravity constant and interpolation, while the Galactic potential retains its original physical normalization.

Independent checks cover gradients of the potential and force, disk symmetry, halo-cutoff continuity, analytic force bounds, coordinate orientation, distance-dependent displacement/velocity changes, energy and axial angular-momentum conservation, refined integration, and time reversal. Comparison to published pericenter and apocenter estimates is descriptive; catalog covariance, Galactic-model uncertainty and the original averaging procedure are not reconstructed.

The sufficient force radius follows a componentwise Hessian bound over a spatial ball, not a finite angular sample. The frozen rotating-frame scale omits the Euler term only for its explicitly diagnostic definition; the actual external force in the nonrotating transport frame remains the exact difference of Galactic accelerations. The response-timescale indicator includes changes in tensor orientation, while a separate amplitude indicator records changes in its norm. Neither is an inferred heating rate.
