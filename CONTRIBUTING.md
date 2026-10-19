# Extending the research

Each experiment should answer one physical question with declared assumptions,
reproducible inputs, numerical verification, and a short record of actual results.
The project is maintained by Venkata Hemanth Athota as independent physics research.
Contributions use the repository's MIT license.

## Add an experiment

1. Assign the next unused identifier after E001, E002, and E003. Start its record from
   [the experiment template](experiments/templates/record.md). State the question,
   model, units, applicable regime, and relevant primary literature before interpreting results.
2. Place reusable dynamics in `src/recoil_dynamics/`. Keep experiment orchestration
   under `experiments/`, with numerical inputs in a versioned JSON configuration.
   The current runner explicitly supports E001 and E002; a new experiment must be
   registered in its argument choices and dispatch, or supplied with its own runner.
   E003 uses `experiments/e003.json`, `run_e003.py`, and `verify_e003.py` separately
   to preserve the published E001/E002 inputs and outputs.
3. Add independent validation appropriate to the physics: an analytic limit,
   conservation law, alternative formulation, or convergence study. Recomputing a
   saved table with the same method is a reproducibility check, not independent validation.
4. Save results under `results/eNNN/`. Document every column's units and conventions.
   Record configuration and source hashes, dependency versions, the generating
   source commit, and output hashes. Include stochastic seeds when applicable.
5. Extend the result verifier and workflow to cover the new outputs. Keep published
   baseline configurations intact; use a separate configuration for a changed experiment.
6. Report measured results and limitations in `docs/eNNN-<topic>.md`, and add the
   record to the README. Distinguish an approximation test from an observational
   constraint or a claim of novel physics.

## Verify a contribution

```sh
python -m pytest -q
python experiments/verify_results.py
python experiments/verify_e003.py
```

For a change to a hashed simulation input, first commit the completed source changes,
then run `python experiments/run.py`, rerun both checks, and commit the generated
results. This makes the manifest point to the generating source commit. Never edit
recorded hashes or numerical results by hand to make verification pass. For E003,
use its separate runner and verifier; its runner rejects uncommitted source inputs.

Inspect regenerated figures and compare numerical changes with their physical cause.
If the technical note changes, rebuild its PDF using [the paper instructions](paper/README.md).
Do not upload raw credentials, private datasets, or third-party articles.

## Preserve a citable result

Use a new version for a changed published result. Update package and citation
metadata together, record the change in `CHANGELOG.md`, and tag the verified commit.
Keep old release tags unchanged. A release identifies a specific snapshot; cite that
version rather than the moving `main` branch. A GitHub release is not a DOI archive.
Any future archival DOI must be added only after it has actually been issued.
