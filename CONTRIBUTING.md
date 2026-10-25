# Contributing to the research formulation

The active scope is [Step 1: recoil transport](docs/e001-ballistic-return.md) and [Step 2: residence and environmental exposure](docs/e002-core-exposure.md). Contributions should clarify, implement, or validate these quantities with explicit physical assumptions.

## Equations and models

State the question, mathematical convention, normalization, units, and validity domain. Distinguish definitions, analytical consequences, closures, numerical results, and empirical calibration. Preserve joint host/merger/launch conditioning.

Gravity and scattering must describe consistent populations. Drag, diffusion, cutoffs, boundaries, and stopping hazards require justification. Missing inputs are not zero-valued defaults.

## Preparing calculations

Use the [handoff](docs/calculation-handoff.md) and [input register](formulation/inputs.json). The general trajectory inputs remain unassigned and trajectory execution disabled. Defined host and local-coefficient stages have separate versioned protocols linked from the register. Before adding a runner, supply versioned host and launch specifications, units, closure, validity criteria, event rules, error controls, and independent validation requirements.

Plan source, configuration, seed, dependency, and output provenance before execution. Report outcomes only after appropriate checks. Numerical convergence is not physical calibration; correlated paths are not independent observations.

Retained numerical modules have their own scopes and immutable citations in [computational reference records](docs/computational-records.md). Their existence does not implement or validate the active diffusion model. Do not modify recorded hashes or result tables to make a new calculation appear reproduced.

## Documentation checks

```sh
python scripts/verify_foundation.py
python paper/verify.py
```

These checks perform no integrations, experiments, or paper compilation. Use fenced `math` display blocks and protected inline notation where needed. Check that GitHub preserves expressions, and keep links and symbols consistent across both steps.

## Manuscript and citation

Title, abstract, scope, and compilation status are defined in paper/working.json; keep them consistent with the manuscript and README. Follow [paper instructions](paper/README.md). The paper builder compiles only LaTeX and does not execute a physical model. Cite the PDF revision actually used.

Cite the precise formulation or calculated record. Keep tags immutable. A GitHub tag is not a DOI or evidence of peer review. Preserve author credit and the MIT license; external sources and datasets retain their rights.
