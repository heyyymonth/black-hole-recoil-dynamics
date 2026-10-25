# Mass-resolved host and local transport calculations

This implementation reconstructs one actual NGC 3201 joint posterior row, checks its local and projected distributions, then evaluates the anisotropic drift and diffusion kernels. It does not fit the observations, evolve a remnant, or measure residence.

The [host protocol](protocol.json) and [collision protocol](../transport/protocol.json) declare the choices and acceptance thresholds. The selected row has the largest posterior quadrature weight, not the largest likelihood, and is not a typical or uncertainty-marginalized host. Original run versions were not recorded in the released chain: the pinned historical packages are compatible reconstruction choices.

## Reproduce

From the repository root, use Python 3.11 in a dedicated environment. This avoids mixing the historical scientific libraries with the data-audit or archival-experiment environments.

```sh
python3.11 -m venv .venv-host
.venv-host/bin/python -m pip install --no-deps -r analysis/host_reconstruction/requirements-lock.txt
.venv-host/bin/python analysis/cluster_data/fetch.py
.venv-host/bin/python analysis/host_reconstruction/run.py
.venv-host/bin/python -m pytest -q analysis/transport/test_coefficients.py
.venv-host/bin/python analysis/transport/run.py
```

All inputs are hashed before use. The dependency lock supplies every installed package explicitly; `--no-deps` prevents historical package metadata from fetching newer default branches. Runtime versions and scientific-module hashes are recorded in [environment.json](environment.json). Host reconstruction takes several minutes on the recorded system. The collision runner refuses a host whose checks fail.

The source stellar-evolution implementation emits a known intermediate divide-by-zero warning near an evolving mass-bin edge. The warning is retained, not suppressed from the research record. All exported physical arrays, masses, densities, moments and coefficients are checked for finiteness. This does not establish that alternative stellar-evolution prescriptions would agree.

`run.py --reuse-solution` repeats validation on the already computed, parameter-matched numerical host, without re-solving the model. A full reproduction uses the default invocation. Stored host arrays use pc, km/s and solar masses; local coefficients additionally use Myr. The source model uses G=0.004302 in these units. It is preserved consistently in both gravity and collisions, rather than silently mixed with the slightly different constant used for catalog display scales.

## Validation boundaries

Mass, force, relative potential, virial balance and all-species distribution-function moments are checked independently. Projection checks target the luminous tracer bin used in the published kinematic products. Diagnostic differences for every species are also retained: a steep, effectively empty black-hole tail makes relative errors of the source's projected arrays uninformative there. These source projected arrays are not used as a local scattering distribution.

Positive density and pressure profiles are interpolated in logarithms. Direct interpolation of tiny positive tails can produce substantial relative error or negative values; velocity integration verifies the adopted interpolation. The last positive mesh point is tapered to zero with the lowered-model asymptotic exponent.

The coefficient integration uses the full anisotropic number distribution, includes both tangential velocity components correctly, and evaluates positive-weight ray integrals centered on the remnant velocity. Analytic uniform-ball limits, rotations, trace identity, positive semidefiniteness, and quadrature refinement test its implementation. Its current host validation grid covers speeds up to 0.75 of the local stellar cutoff. It is not a validated coefficient table for arbitrary trajectories.

Results are per unit constant Coulomb logarithm. No impact-parameter range, strong-encounter treatment or Gaussian continuum applicability is inferred. The few-black-hole contribution must be examined before trajectories can be interpreted physically. A finite distribution-function truncation radius is not a measured tidal escape boundary.
