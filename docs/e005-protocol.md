# E005 protocol: late-time availability of binary-black-hole targets

**Researcher:** Venkata Hemanth Athota. Declared 21 October 2026, before reading or calculating the catalog's BBH-count outcomes. E004's already published snapshot and counts are prior information.

## Question and connection to E004

In the available late-time cluster summaries, how widely does the small instantaneous BH–BH population found in E004 occur, and does this release contain the information needed to calibrate companion access?

E001/E002 establish how ballistic motion changes residence and density exposure. E003 establishes that access also depends on relative velocity in a specified target bath. E004 tests transfer of that bath and host potential to one evolved cluster: only two actual BH–BH targets are present, and matching host mass and half-mass radius does not match ballistic returns. E005 addresses the breadth of the **target-count constraint** before attempting a continuous distribution fit. It does not repair the host potential or test live motion.

This is a simulation-release audit, not a gravitational-wave observation or an experiment producing new cluster simulations. The population relevant to a recoiling remnant is conditional on a preceding merger. Equal weighting of available cluster models at a fixed old age does not reproduce that event-conditioned population.

## External inputs and selection fixed before analysis

Use the author-maintained [cmctoolkit release](https://github.com/NicholasRui/cmctoolkit/tree/a4876e7706f7b23a14c63a5be7df63a07f2dccd9) at commit `a4876e7706f7b23a14c63a5be7df63a07f2dccd9`, the same source as E004. Inventory metadata show 7,537 `CMC/<model>/PARAM*.txt` files in 126 standard model directories. Include all of these paths, fixed in [the input lock](../experiments/e005-inputs.json), with their Git blob hashes and sizes. Exclude the separately tuned `CMC/NGC6624/` models (311 records, six histories); they are outside this declared comparison. No model is selected by its BH abundance.

The source README describes summaries older than 10 Gyr. Verify the actual age range rather than enforcing that description through silent filtering. The published [CMC catalog](https://arxiv.org/abs/1911.00018) contains 148 simulations; this late-time release exposes only 126 standard directories. Do not assign zero counts to unavailable histories or assume why they are absent. The archive has one full object snapshot, already used in E004; PARAM records contain aggregate properties, not joint binary phase space.

Fetch the pinned source archive using verified HTTPS. Read only locked PARAM members, verify file size and Git blob identity, and record SHA-256 for each consumed input. Cache original data outside version control; do not redistribute or relicense upstream files.

## Measurements and comparisons

1. Parse age, host mass, system/object counts, all binary type-pair counts, and all stellar-type counts. Reject duplicate/missing required keys, nonfinite or negative counts, noninteger counts, and inconsistent accounting. Check `sum(N_type)=Nobj`, `sum(Nbin_i_j)=Nbin`, and `Nsys+Nbin=Nobj`. Derive single BHs from total BH objects after subtracting binary membership; reject negative residuals. Unused summary fields may contain documented `nan` values for absent types.
2. For each snapshot retain model ID, source path, age, BH–BH count `Nbin_14_14`, BHs with other companions, single BH count, all-binary count, and host mass. A BH–non-BH binary is not an E003 BH–BH target.
3. For the primary cross-model comparison, take **one snapshot nearest 12 Gyr per model**, within **0.25 Gyr** inclusive; break ties by earlier age and then path. Twelve Gyr is a declared old-cluster comparison age, not an assumed merger epoch. Record models with no eligible record as missing, never zero. Do not interpolate integer populations between snapshots.
4. Report the selected sample size, zero-BBH states, states below the carried-over count-30 screen, and the median/range of instantaneous counts. Report the corresponding BH and all-binary counts to expose population substitution. Treat each model equally for this descriptive comparison; this is not a probabilistic or astrophysically weighted sample of the Universe.
5. For each available history report its number of records, age range and minimum/median/maximum BBH count, and the number of records passing the count screen. These describe saved records: **no time-average, residence duration, or probability** is inferred from snapshot frequency. Do not pool records as independent observations or sum counts as distinct systems. IDs are absent from PARAM.
6. Cross-check the already analyzed E004 state against its upstream values and the independently parsed full-snapshot results. Locate that state's count within its own available history, separately from the fixed-age primary comparison. The E004 state is not an independent new validation sample.

The count of 30 remains an operational screening rule inherited from E004. Passing it is neither a sufficiency theorem for a high-dimensional distribution nor a demonstration that a continuous bath is locally accurate. Failure does not establish that a physical kinetic description is invalid, that a cluster has no future companions, or that hierarchical mergers are absent.

## Validation and interpretation rules

Test parsers and reducers with synthetic accounting examples, zero counts, missing age coverage, corrupt inputs, and age-selection ties. Independently check the main census using a separately written scalar extraction/counting path. Recompute every published table and verify source/input/output hashes. Record numerical outcomes only after source and this protocol are committed.

If many states are sparse, the next fit must explicitly model discrete target availability and verify local support. If abundant states exist, those histories are candidates for obtaining full object data, not automatically suitable calibration samples. Either outcome still requires joint positions, velocities, masses and binary orbits, stable object identities, merger/interaction histories, and held-out cluster tests. Selection must ultimately be conditioned on remnant-producing mergers. Existing [CMC methods](https://arxiv.org/abs/2106.02643) sample orbital phases; joining snapshots does not validate a recoil trajectory.

Do not substitute the census mean for E003's 100 targets, multiply E003 and E004 discrepancy factors, calculate encounter rates from aggregate counts, or claim live friction/recoil calibration. The scientific advance of E005 is to determine what the available target inventory can and cannot support, and specify the next measurement accordingly.
