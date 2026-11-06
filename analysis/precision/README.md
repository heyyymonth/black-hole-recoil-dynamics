# Residence-response precision pilot

This package compares three estimators of the same complete-transition residence response: eigenbasis coupling with uniform insertion, symmetric-root coupling with uniform insertion, and symmetric-root coupling with proportional time-stratified insertion. The [protocol](protocol.json) fixes the design; the [scientific description](../../docs/precision-pilot.md) gives the identities, uncertainty accounting and interpretation.

The conditional host, physical encounter law, fine table, initial states, horizon, rewards and stopping rules are inherited from the completed [response experiment](../response/README.md). Its hash-bound implementation and results are not edited. This pilot changes only within-quartet Gaussian coupling and insertion sampling. Correct marginal Gaussian covariance does not imply identical realized paths or equal estimator variance.

## Execution and retained evidence

Use the repository's pinned numerical environment. From the repository root:

```sh
python -m pytest -q analysis/precision
python analysis/precision/run.py --resume
python analysis/precision/analyze.py
python analysis/precision/verify.py
```

The runner writes a fixed execution design before any outcomes and checks the frozen response identity. It executes four blocks of 1,024 quartets for each of three launches and three methods, for 36 jobs and 36,864 independent quartets. Each quartet retains four correlated branch trajectories. No sample-size adaptation, tail deletion, or outcome-dependent method selection occurs during execution. The protocol rotates method and launch order between blocks; every job uses an independent seed.

`--resume` accepts a fresh output directory or verifies existing complete records before continuing. It refuses inconsistent inputs, changed frozen sources, stale sample hashes or incomplete record pairs. Do not edit the protocol or bound implementation to resume incompatible results. An incomplete job requires explicit investigation; completed evidence must not be silently overwritten.

Files under `results/precision/` retain per-job NPZ arrays and JSON metadata, the execution design, and a summary produced only from the full inventory. Every quartet includes branch states, normalized residence, exposure, termination time and cause, event counts, insertion step/state, zero/positive Poisson weights and the recomputable response score. Stratified records additionally retain their fixed stratum labels. Numerical-domain missing-residence allowances and computing-time measurements are recorded separately.

## Components

- `coupling.py` loads the frozen response transport into a private module namespace. `make_transport('eigen')` retains its exact Gaussian update; `make_transport('symmetric')` changes only that private instance's covariance factor. The published imported module is never mutated.
- `run.py` checks identities, assigns disjoint streams, implements the fixed insertion design and retains all branch outcomes. It reuses the tested lossless disk sampler cache under `tmp/response-sampler-cache/fine/`.
- `statistics.py` computes the design-correct mean variance and empirical efficiency ratios. Stratified uncertainty uses within-stratum variances and fixed weights, not the pooled raw-score variance.
- `analyze.py` reads the complete retained inventory, produces per-block and pooled diagnostics, and computes the six prespecified efficiency comparisons. It executes no trajectories.
- `verify.py` uses NumPy and the standard library to check the source/input chain, job inventory, chronology, seeds, exact insertion assignments, reward and stop accounting, and reproducibility of the summary. It executes no trajectories and requires no SciPy reconstruction.

The sampler cache is disposable and reconstructed from unchanged physical inputs. Its contents are not scientific evidence and are not bundled with the retained records. It preserves numeric sampler arrays and random-number calls; it is not another physical approximation.

## Uncertainty and scope

The fixed design gives 4,096 quartets per launch and method. Four equal insertion-time strata each receive 256 rows per block. Uniform insertion samples all 1,000 scheduled steps; stratified insertion samples uniformly within each assigned interval. Dead-before-insertion paths remain zero-valued members of their original sample and stratum.

The nine response means and six independent-arm mean differences are reported with descriptive standard errors. They do not establish equivalence. The six efficiency ratios use 8,192 independent within-design bootstrap replicates with nominal Bonferroni component levels. These empirical intervals are approximate, can miss unseen rare contributions, and do not certify finite-sample family coverage.

The comparison has equal quartet counts, not equal CPU budgets. The recorded variance-times-CPU metric describes efficiency on the measured machine/cache/load. Bootstrap cost-adjusted intervals condition on the observed timing; timing variability is excluded. Wall time, CPU time and cache activity remain available for interpretation.

A precision gain would support choosing an estimator for a further independent response test. It would not by itself establish endpoint prediction accuracy, convergence of the numerical process, or calibration to an observed cluster. The previously completed response benchmark and its independent validation remain distinct evidence.

## Retained findings and post-pilot diagnostic

The completed comparison supports symmetric-root coupling as a much more precise reference estimator at the measured cost. All three E/S intervals lie above one; all three S/T intervals include one. Uniform insertion is retained as the straightforward reference for a further independent test, without claiming that it has lower population variance than proportional stratification. [Numerical findings and limitations](../../docs/precision-pilot.md#completed-findings).

A separate `tails.py` diagnostic was specified after the primary comparison. It decomposes the design-correct variance estimate into nonnegative centered-square contributions and records their concentration in `results/precision/tail-diagnostic.json`. It changes no pilot scores, sample counts, confidence family or summary. The reported concentration concerns the observed sample, not population-tail probabilities. Reproduce or check it with:

```sh
python analysis/precision/tails.py
python analysis/precision/tails.py --verify
```
