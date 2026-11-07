# Independent symmetric-root residence-response test

The [protocol](protocol.json) fixes an independent test of the finite-step residence response with symmetric-root coupling. [Scientific description](../../docs/independent-response-test.md) gives the hypothesis, confidence allocation, retained budgets and limits of the numerical continuation bounds.

Use the repository's pinned Python 3.11 numerical environment. From the repository root:

```sh
python -m pytest -q analysis/retest
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python analysis/retest/run.py --stage all --resume
python analysis/retest/analyze.py
python analysis/retest/verify.py
```

The runner has separate `prediction` and `validation` stages. Validation requires every prediction and an immutable prediction register. Existing complete records are accepted only under matching source, input, design, job, dependency and array hashes. Incomplete record pairs stop execution for investigation; they are not silently overwritten.

The two workers have separate sampler caches under `tmp/retest-sampler-cache/`. They may read hard-linked copies of immutable entries from the earlier exact cache, but never share a writable cache directory. Missing entries are constructed using the unchanged physical sampler. Scratch caches can be discarded; the NPZ results, JSON records, execution design and prediction register form the evidence.

`run.py` reuses the frozen complete-transition solver through the tested private symmetric-root module. `analyze.py` recomputes all means, residuals, refinement comparisons, tail diagnostics and continuation bounds without executing trajectories. `statistics.py` contains independent-batch variance formulas and the binomial occurrence bound. `verify.py` reconstructs seeds and job inventory independently, recomputes reward/score/stop accounting, checks all numerical controls against hash-bound reference records, and enforces the stage gate while allowing independent workers to overlap.

The continuation quantile implementation requires only NumPy and the standard library; independent tests compare it with beta quantiles from SciPy. Manuscript verification requires no host reconstruction or dynamics. No previous prediction, precision-pilot or validation samples are pooled into the current estimates. Source and sample counts are fixed before production.

## Completed evidence

All 24 fixed ensembles are recorded in `results/retest/`. The independent audit verifies 368,640 groups and their 1,228,800 correlated branches, including prediction registration before validation. See [recorded findings](../../docs/independent-response-test.md#recorded-findings) for the intervals, tail diagnostics and continuation bounds. All primary accuracy and refinement classifications remain inconclusive. The sampling improvement does not establish physical-cluster calibration.
