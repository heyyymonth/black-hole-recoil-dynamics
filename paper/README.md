# Research manuscript

**Recoil transport and environmental exposure of black-hole merger remnants**

Venkata Hemanth Athota · Independent researcher<br/>
bharadwaj.hemanth.002@gmail.com

[PDF](../output/pdf/central-residence.pdf) · [LaTeX source](manuscript.tex) · [Portable source](../output/pdf/manuscript-source.zip) · [Citation](CITATION.bib)

The manuscript develops two steps: merger-conditioned recoil transport, and finite-time residence and environmental exposure. It specifies host distributions, launches, drift and diffusion, event definitions, stopping, uncertainty, and the calculation handoff. It reports no numerical transport results or calibrated cluster model.

Revision **0.4.0-manuscript.1** is an unreviewed theoretical manuscript. No submission, acceptance, or DOI is claimed. The [claim-to-evidence record](CLAIMS.md) distinguishes definitions, established theory, and unexecuted validation requirements. [working.json](working.json) records the title, abstract, and scope.

## Build and verify

The article uses the two-column APS/PRD REVTeX format with its standard typography. See [FORMAT.md](FORMAT.md).

With Python 3.11, the paper dependencies, and Tectonic 0.16.9 installed:

```sh
python -m pip install -r paper/requirements.txt
python scripts/verify_foundation.py
python paper/build.py
python paper/verify.py
```

The builder typesets LaTeX and BibTeX only. It does not import numerical research modules, run experiments, select inputs, generate result plots, or evaluate transport equations. It writes the PDF, portable source ZIP, and build manifest under `output/pdf/`.

Hashes bind the article, bibliography, formulation, handoff, and build metadata. Verification checks citations and archive consistency, page dimensions, embedded fonts, and unresolved references or overfull equations. The archive includes the pinned unmodified REVTeX distribution; other standard TeX packages remain external dependencies.

Render and inspect every page after a layout change. Source hashes and successful compilation do not substitute for visual inspection or scientific validation.

Use [manuscript-v0.4.0.1](https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/manuscript-v0.4.0.1) for an exact citation. Computational reference records have separate [immutable citations](../docs/computational-records.md).
