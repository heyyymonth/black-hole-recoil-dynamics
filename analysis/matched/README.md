# Matched finite-encounter and diffusion residence

This conditional Omega Centauri calculation uses the same positive BH encounter law for finite jumps and its first two velocity-change moments. It varies the encounter representation while retaining the host gravity, remnant mass, launches, cutoffs, region, horizon and stopping conventions. Lighter-species collisions and host response are disabled. The [results discussion](../../docs/matched-residence.md) separates numerical findings from continuum accuracy and physical calibration.

The [protocol](protocol.json) was recorded before running the ensembles. It fixes two grids, two steps, two launch directions, two arms and 4096 paths per setting. Diagnostic zero-event counts are summarized after execution; they were not a preregistered hypothesis test.

## Reproduce

Use the repository Python 3.11 environment and pinned `requirements.txt`. The existing host is checked separately by the six-cluster verifier. From the repository root:

```sh
python analysis/matched/build_kernel.py --grid coarse
python analysis/matched/build_kernel.py --grid fine
python analysis/matched/validate_sampling.py
python analysis/matched/controls.py
python analysis/matched/run.py --grid coarse
python analysis/matched/run.py --grid fine
python analysis/matched/analyze.py
python analysis/matched/verify.py
python analysis/matched/controls.py --verify
python paper/matched_assets.py
```

`--resume` on the table builder and ensemble runner accepts only matching source, protocol and output hashes. Tables retain every quadrature refinement; unaccepted local quadrature prevents trajectory execution. The finite-angle CDF cache contains up to 128 vertices and can use roughly 1.5 GB. Finite-event runs take longer than Gaussian runs because new local samplers are constructed on demand. These are independently resampled encounters, not a simulation of persistent stars.

`verify.py` and `controls.py --verify` require only NumPy and the standard library. They inspect hashes, recorded numerical criteria, path accounting and summaries without integrating trajectories. The source archive includes retained per-path final states, residence, exposure, stopping times and event counts; it does not retain complete stochastic histories. Sampler diagnostics retain summary moments and standard errors, not every sampled increment.

## Numerical law

`kernel.py` mixes full positive vertex rate measures. A sampled event chooses a vertex in proportion to interpolation weight times allowed event rate. Both arms use that mixture's drift and raw second moment. Local symmetry averaging and frame transformations are shared by the sampler and moments.

`transport.py` freezes the law at each step start, draws a common-form excluded-event clock, and applies a Gaussian or compound-Poisson increment between two gravity half-flows. Multiple Poisson events are summed at the midpoint. It integrates exact sphere occupation along each drift chord, samples density at drift midpoints, and resolves validity stops at the declared split locations. The horizon denominator is retained for stopped paths. Grid/time refinement is needed; matching frozen substep moments alone does not establish residence convergence.

`statistics.py` reports independent-ensemble normal mean intervals, Wilson stop intervals and union-DKW CDF-distance intervals. The mean/Wilson intervals have approximate coverage. All reported intervals are marginal, not simultaneous over the whole pilot. A numerical exit supplies a remaining-time bound on missing residence; it is not physical escape.
