# Six-cluster comparison

Repeat the established conditional host reconstruction, local scattering and environmental calculations for 47 Tuc, NGC 1851, NGC 3201, Omega Centauri, M4 and NGC 6397. Keep separate source-quality, numerical-consistency and physical-applicability decisions. No remnant trajectory or residence calculation is performed.

[Results](../../docs/six-cluster-comparison.md) · [Data suitability](../../docs/six-cluster-data-quality.md) · [Environmental results](../../docs/six-cluster-environment.md) · [Protocol](protocol.json)

Use the pinned scientific environment in [host reconstruction](../host_reconstruction/README.md). From the repository root:

```sh
.venv-host/bin/python analysis/cluster_data/fetch.py
.venv-host/bin/python analysis/cluster_comparison/posterior_inventory.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-host/bin/python analysis/cluster_comparison/run.py --stage reference --workers 2
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-host/bin/python analysis/cluster_comparison/run.py --stage pilot --workers 2
.venv-host/bin/python analysis/cluster_comparison/environment.py
.venv-host/bin/python analysis/cluster_comparison/run.py --stage collect
.venv-host/bin/python analysis/cluster_comparison/report.py
.venv-host/bin/python analysis/cluster_comparison/environment_report.py
.venv-host/bin/python analysis/cluster_comparison/verify.py
```

`--clusters` accepts explicit source identifiers for partial reruns. The reference is one maximum-quadrature-weight joint row. Eight fixed-seed weighted joint rows per cluster are an exploratory sensitivity pilot, not a converged posterior calculation. Selection records retain the full parameter vectors, source weights and row indices. Failures are retained without replacing a draw; zero-BH models remain valid outcomes.

Shared host and Rosenbluth implementations are used directly. Reference runs retain all 20 original velocity/radius/orientation states and numerical checks. Pilot runs validate reconstructed hosts and evaluate nine zero-speed radii, without repeating the velocity grid. At zero velocity, the diffusion trace is evaluated by an independent axisymmetric two-dimensional integral, checked under 64-to-128 quadrature refinement. Its central limit and original-host agreement have independent tests.

Reference host arrays are retained for reproducible environmental calculations. The accepted pilot meshes are rebuildable scratch in `tmp/cluster-comparison/`; their hashes, parameter vectors, numerical summaries, species tables and radial coefficients remain in the result records. Keeping scratch meshes out of the repository avoids hundreds of megabytes of redundant numerical arrays. Cache reuse requires matching parameters and protocol and reruns numerical validation. `--resume` skips existing attempt records, including failures, and is intended only for an unchanged run configuration. To investigate a failed reconstruction explicitly, use `--stage retry-failed` without `--resume`; the preceding failure remains in the attempt history. A known finite-extent failure is never silently reused as a fresh calculation.

Published derived arrays contain 4,000 model realizations per cluster. Their generating-row mapping is unavailable, so they are never joined by index to raw nested rows. Their empirical intervals and absent-BH cases are reported separately from the reconstructed pilot.

The environmental inputs preserve one pinned source table and its actual per-cluster phase points. A clean checkout uses the six committed exact HTML rows and their parsed metadata, checked against the recorded hashes. If the full cached table exists, its original hash must match before extraction. No fresh upstream page is silently substituted. The common logarithmic static controls and Irrgang Model I orbit integration are conditional assumptions. Tolerance refinement never changes acceptance thresholds. A missing preferred physical model or covariance is not repaired with invented values.

The notebook reads recorded results and plots editable comparisons; it does not rerun the expensive host solver. Execute it with the data-audit environment described in [the original data package](../../data/observed-clusters/README.md).
