# Discrete residence-response experiment

This package tests the first derivative of the complete jump/diffusion transition mixture against independent finite-change trajectories. The [protocol](protocol.json) fixes the controls and tolerances; the [scientific description](../../docs/response-test.md) explains the observable, estimator and scope.

`transport.py` implements complete stopped transitions, paired mixtures and the four-branch weighted replacement estimator. `cache.py` stores the existing matched encounter measure compactly without changing its sampler or random-number sequence. `statistics.py` calculates sampling intervals and prediction residuals. Exact polynomial-moment controls distinguish a correct first derivative from an accurate endpoint extrapolation. The earlier matched calculation and its hash-bound records remain independent inputs.

Use the repository's pinned numerical environment. From the repository root:

```sh
python -m pytest -q analysis/response
python analysis/response/run.py --stage pilot --resume
python analysis/response/run.py --stage prediction --resume
python analysis/response/run.py --stage validation --resume
python analysis/response/analyze.py
python analysis/response/verify.py
python paper/response_assets.py
```

The stages are ordered. Pilot variance selects and locks production sample counts before production begins. All production predictions are then saved and hash-registered before validation starts. Resume is permitted only when protocol, implementation, source inputs and existing sample hashes agree; existing incompatible records must not be silently reused.

The retained execution also has a separately recorded storage continuation:

```sh
python analysis/response/execute_cached.py --stage all --resume
```

This command requires the existing pilot, count lock and records named in [storage-adapter.json](../../results/response/storage-adapter.json). The manifest preserves the hashes of 14 records produced before the adapter and assigns the other 22 registered jobs to the adapter. Their individual metadata records retain the adapter-manifest hash. `disk_cache.py` reuses exact sampler arrays on disk; it preserves the transition law, numerical values, random-number calls, seeds, sample allocations and statistical thresholds. Independent storage-equivalence checks are retained with the manifest. This is a storage optimization within the same experiment, not an additional physical comparison.

The cached continuation calls the frozen original runner and transport implementation. Neither those sources nor the storage-adapter sources may be edited during execution: resume verifies their recorded hashes and refuses incompatible state. The cache under `tmp/response-sampler-cache/` is disposable and reconstructible from the frozen inputs. It is not scientific evidence; the versioned per-job samples, metadata and stage registers are the evidence.

The seed vector is `SeedSequence([2026110201, stage, setting, launch, method_or_a])`, with stage codes 10, 20 and 30. Indices follow the protocol order. Pilot, production prediction and validation streams are disjoint. Pairing occurs only within a constructed response quartet or mixture/control pair, and every branch retains its own state-dependent hazard.

The stored records retain final states, residence and exposure, termination times and causes, event counts and evaluated steps. Response records additionally retain the insertion step and state, branch weights, conditional count and every branch's rewards. These permit independent recomputation of the score and numerical-domain continuation allowance. A saved scalar mean alone is insufficient evidence.

The data analysis and paper build must read recorded samples; neither should rerun trajectories. Numerical execution may be expensive. A sample cap or unresolved interval is reported as such without widening the protocol tolerances or treating an unexecuted comparison as a result.

`verify.py` uses only NumPy and the standard library. It independently checks the full stage inventory, hashes, seed separation, variance-only count selection, retained branch accounting and summary arithmetic. Local UTC timestamps corroborate the hash chain; this is not an externally timestamped preregistration.

The exploratory `coupling_diagnostic.py` examines retained local covariance matrices after the pilot, without evolving trajectories. Its `--verify` mode checks equal marginal covariance and common-noise distances from stored matrices. It does not establish the cause of the residence-score variance or change the registered estimator.
