# E005: how many BH–BH targets are actually available?

**Researcher:** Venkata Hemanth Athota. Executed 21 October 2026. This is an audit of released CMC simulation summaries, not new cluster simulations or astronomical observations.

## Question and answer

**Is E004's small target population confined to that one snapshot, and can the available release support calibration of a continuous companion distribution?**

Small instantaneous BH–BH populations occur throughout the inspected late-time release. Among 119 model histories with a saved state near 12 Gyr, the median is **two BH–BH binaries**, 38 states contain none, and the maximum is nine. Across all 7,537 inspected states the maximum is 15. The release therefore does not pass the inherited count-30 screen in any saved state. More fundamentally, its aggregate summaries do not provide the joint binary positions, velocities, internal orbits, and identities needed for the proposed calibration.

This broadens E004's target-count evidence within the available late-time model grid. It does not show that E004 is representative of the environments in which merger remnants are produced. It also does not prove that a smooth ensemble encounter model is invalid: sparse individual systems can still have useful ensemble descriptions.

[Pre-execution protocol](e005-protocol.md) · [Configuration](../experiments/e005.json) · [Summary](../results/e005/summary.json) · [All states](../results/e005/snapshots.csv) · [Per-history summaries](../results/e005/models.csv) · [Age-selected states](../results/e005/selected.csv) · [Provenance](../results/run_e005.json)

[Cite E005](../results/e005/CITATION.bib) using [immutable results snapshot 80d82f8](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/80d82f85d1c84fa7ef42fc40f71a1cbc5439c172). This experiment record is separate from the current E001–E004 paper and has no DOI. Cite the original simulation and toolkit research when reusing their data.

## Why this is the next experiment

| Earlier result | Question it leaves | E005's contribution |
|---|---|---|
| E001/E002: ballistic motion changes residence and density exposure. | Exposure to which actual companions? | Keep the target population explicitly defined as BH–BH binaries. |
| E003: relative velocity changes passage predictions in an imposed bath. | Can that bath's abundance and joint properties be calibrated? | Check abundance and data availability before fitting a velocity-dependent target model. |
| E004: the selected snapshot contains only two BH–BH targets. | Is that scarcity limited to the selected state? | Examine every standard-model summary in the same pinned release, with an age-matched comparison. |
| E004: matched Plummer and frozen-shell return times differ substantially. | Does the trajectory prescription describe a live cluster? | Leave this as a separate unresolved test; counts cannot validate motion. |

More kicks in E003 would reduce interpolation gaps in the same imposed bath. More independent cluster histories address whether that bath has empirical support. These are different sources of information. Here the natural sampling unit is the model history, not each saved row or each repeated appearance of a binary.

## Data, selection and accounting

Use [Nicholas Rui's cmctoolkit](https://github.com/NicholasRui/cmctoolkit/tree/a4876e7706f7b23a14c63a5be7df63a07f2dccd9), pinned to the same revision as E004. Inventory selection was fixed before reading the new count outcomes in protocol commit `72916d1`; generating code was committed as `55669f8`. The [input lock](../experiments/e005-inputs.json) contains 7,537 file paths and exact Git blob identities from 126 standard model directories. Six separately tuned NGC 6624 models are excluded. The published [CMC catalog](https://arxiv.org/abs/1911.00018) contains 148 simulations; this particular late-time summary release does not expose all of them. Unavailable models are not assigned zero populations.

For the primary comparison, select the saved state closest to 12 Gyr within 0.25 Gyr, with ties resolved by earlier age and then path. Seven histories have no eligible state. Their names remain in the summary and their selected fields are empty in the model table. Actual selected ages span 11.8013–12.2086 Gyr. Every age and signed offset is retained. Equal model weights give a descriptive comparison, not a cosmic population estimate or an assumed merger-time distribution.

All 7,537 records pass count-accounting checks: stellar-type counts sum to the number of objects; binary-type pair counts sum to all binaries; systems plus binaries equal objects; and BH membership leaves a nonnegative single-BH count. In particular, a BH–BH binary contains two BH objects but is one encounter target. BHs with other companions remain a separate category.

## Results

| Quantity | Executed result | Meaning |
|---|---:|---|
| Available standard model histories | 126 | Scope of this release, not all catalog simulations |
| Saved summary states | 7,537 | Correlated states, not 7,537 independent clusters |
| Saved-state ages | 10.0007–14.0022 Gyr | Simulation output ages, not a cosmological-age selection |
| Histories with a state near 12 Gyr | 119 | One selected state per history |
| Selected BH–BH counts: minimum / median / maximum | 0 / 2 / 9 | Instantaneous target systems |
| Selected states with zero BH–BH binaries | 38 / 119 | Measured zeros, separate from seven missing ages |
| Selected states with at most two BH–BH binaries | 71 / 119 | Small counts are not confined to E004 |
| Selected states meeting the count-30 screen | 0 / 119 | No single-state fit passes the operational screen |
| Maximum BH–BH count in any inspected saved state | 15 | No saved state anywhere in this audit passes that screen |
| Selected total BH-object counts: minimum / median / maximum | 0 / 38 / 1,878 | Total BHs cannot be substituted for BH–BH targets |
| Selected all-binary counts: minimum / median / maximum | 2 / 25,867 / 142,232 | Most binaries are outside the declared target class |

![Late-time target counts and comparison with all BH objects](../results/e005/target_counts.png)

Figure: each selected state contributes once. The count-30 line is an operational calibration screen, not a physical boundary. The right panel compares different populations in the same state; it does not fit an encounter-rate relation. Overlapping points are possible. [Vector figure](../results/e005/target_counts.pdf).

E004's own model has 26 available saved states, containing between one and three BH–BH binaries, with a saved-state median of two. Its near-12-Gyr selection has three targets at 11.9449 Gyr; the original E004 state has two at 12.7323 Gyr. Thus the original finding is not confined to one saved state of that model. These 26 entries do not identify 26 independent populations or distinct binaries, and they do not establish how long a particular binary survived.

The E004 state also provides a cross-check: independently derived full-snapshot values agree with the catalog parser for age, host mass, systems, all binaries, single BHs, BH-bearing binaries and BH–BH binaries. This is seven checks on the same simulation state, not seven new clusters.

## What changes in our interpretation

E003 remains a valid controlled test of density-only versus velocity-weighted passages. Its 100 identical targets are not calibrated to these late-time states. Replacing 100 with the census median would still omit the correlations between target count, concentration, velocities, component masses, semimajor axes and the host potential. The E003 and E004 discrepancy factors therefore cannot be multiplied into a corrected astrophysical rate.

The next question becomes:

> **Conditional on a merger that leaves a retained remnant, does a model that follows the number and properties of available BH–BH targets predict companion access better than a smooth, stationary bath over the remnant's displacement interval?**

For example, let \(\mathcal M\) denote that preceding merger, its host state, and remnant properties; let \(N_{\rm BBH}(t)\) be the number of available BH–BH targets. The law of total expectation separates two quantities that must be measured:

\[
\mathbb E[\Gamma(t)\mid\mathcal M]
=\sum_n P[N_{\rm BBH}(t)=n\mid\mathcal M]\,
\mathbb E[\Gamma(t)\mid N_{\rm BBH}(t)=n,\mathcal M].
\]

Here \(\Gamma\) is a close-passage rate under a declared encounter definition, not a merger rate. The first factor describes target availability; the second needs phase space, binary properties and remnant motion. This identity states the measurement problem; E005 estimates neither factor for actual merger remnants. A state with no targets has no immediate passages to that target class, but later binary formation can change that. A mean passage count alone does not determine the probability of one or more encounters without a justified counting model.

Existing research already supplies relevant physical mechanisms. [Breen & Heggie](https://arxiv.org/abs/1304.3401) couples the segregated BH subsystem to the host's energy balance. [CMC methods](https://arxiv.org/abs/2106.02643) include binary formation and strong interactions. [Atallah et al.](https://arxiv.org/html/2211.09670v4) combines encounters and recentering in a fixed nuclear-cluster framework. Our present contribution is a reproducible audit of assumptions for a reduced accessibility model, not discovery of BH segregation or a new hierarchical-merger mechanism. Nuclear-cluster results cannot automatically be transferred to this old-globular-cluster sample.

## The next comparison and what would count as progress

The next empirical study needs full snapshots and encounter/merger histories for a declared set of hosts and epochs, selected without viewing an accessibility discrepancy. Anchor the target inventory to remnant-producing merger events, including younger epochs where needed, and distinguish retained and escaping remnants. Track binary component identities, formation, exchange and destruction, together with joint radius, velocity, masses, semimajor axis and eccentricity. Preserve actual zero-target states and time coverage. Hold out whole model histories when evaluating population predictions.

Compare a smooth stationary bath with a model that accounts for changing target counts and joint properties. First compare passage predictions using the **same** independently supplied remnant motion and encounter definition, isolating the population approximation. Then compare frozen and live motion with the same population treatment. Use common windows, an explicitly defined core boundary, uncertainty from independent histories, and recorded encounter counts. Do not fit and validate on repeated appearances of the same binary. Passing the old count screen alone is insufficient; sparse states may require an ensemble treatment rather than a single-state continuous fit.

This study is **not executed here**. The inspected release contains only the E004 full snapshot. The main catalog portal could not be retrieved with verified HTTPS during this check; the reproducible census uses the accessible pinned GitHub archive. CMC's [output documentation](https://clustermontecarlo.github.io/CMC-COSMIC/output/index.html) distinguishes snapshots from interaction and merger logs. Orbit-averaged snapshot positions cannot be joined into a live recoil path. Resolving that trajectory question still requires suitable collisional integrations or independently resolved trajectories.

## Reproduction and limits

```sh
python experiments/fetch_e005.py
python experiments/verify_e005.py
```

The initial fetch downloads an approximately 79 MiB archive and reads only the 7,537 locked PARAM files. Exact Git blob checks precede analysis; the manifest also records every input's SHA-256. Subsequent verification works from the local cache. To regenerate, commit source inputs first, then run `python experiments/run_e005.py`. Original upstream files remain outside version control and retain upstream rights.

Validation includes 67 synthetic parser/selection tests and five archive-integrity tests, seven E004 cross-checks, full table recomputation, and a separately implemented scalar extraction and age-selection census. The complete repository suite now has **180 passing tests**. No inferred confidence interval treats the model grid as a random astrophysical sample. No live recoil trajectory, friction law, companion-acquisition efficiency, merger rate or observational discriminator is validated by E005.
