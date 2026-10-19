# Technical note BHRD-TN-001

**Central residence and density exposure of recoiling remnants in a Plummer potential**<br/>
Venkata Hemanth Athota · Independent researcher · Version 0.1.0

[Read the PDF](../output/pdf/central-residence.pdf) · [Edit the manuscript](manuscript.md)

The paper records E001 and E002. It is an unreviewed technical note accompanying
a reproducible numerical benchmark. It contains no observed-event analysis.

## Build

Use the Python 3.11 environment described in the repository README:

```sh
python -m pip install -r paper/requirements.txt
python experiments/verify_results.py
python paper/build.py
python paper/verify.py
```

The build uses the committed manuscript and result figures, writes the six-page
PDF to `output/pdf/central-residence.pdf`, and records input/output hashes in
`output/pdf/build.json`. It checks the page count, author, and contact text.
Fonts come from the pinned Matplotlib package; equation images are generated
with its math renderer. Intermediate equation images stay in ignored `tmp/pdfs/`.
The original Markdown remains readable on GitHub. The small renderer supports
only the Markdown constructs used in this manuscript.

For a visual check, render the PDF with Poppler and inspect all six pages:

```sh
pdftoppm -r 120 -png output/pdf/central-residence.pdf tmp/pdfs/page
```

The PDF uses fixed document metadata for deterministic output in a fixed build
environment. Rebuild after changing the manuscript, results, or citation metadata.
Do not edit the generated PDF or provenance hashes directly. Numerical results
are reproduced independently of this typesetting step.
