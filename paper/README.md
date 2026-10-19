# Research manuscript

**Central residence and density exposure of recoiling remnants in a Plummer potential**

Venkata Hemanth Athota · Independent researcher · Version 0.1.1

[Paper PDF](../output/pdf/central-residence.pdf) · [LaTeX source](manuscript.tex) · [Portable source archive](../output/pdf/manuscript-source.zip) · [Formatting conventions](FORMAT.md)

The manuscript records E001 and E002. The two-column physics layout follows the
REVTeX family used by the recoil-paper reference. It contains no watermark,
publication banner, or acceptance claim. Review status is recorded here, outside
the paper: this is an unreviewed manuscript and has not been submitted to a venue.

## Reference build

Use Python 3.11 and the numerical dependencies from the root README. Install
[Tectonic 0.16.9](https://github.com/tectonic-typesetting/tectonic/releases/tag/tectonic%400.16.9)
for your platform, then run from the repository root:

```sh
python -m pip install -r paper/requirements.txt
python experiments/verify_results.py
python paper/build.py
python paper/verify.py
```

The builder reads the committed CSV results into native LaTeX tables, regenerates
three vector figures, compiles the article and BibTeX references, and writes:

- `output/pdf/central-residence.pdf`
- `output/pdf/manuscript-source.zip`
- `output/pdf/build.json`

The build pins Tectonic 0.16.9 and the SHA-256 of the official REVTeX 4.2f
distribution in `toolchain.json`. It stages that class and bibliography style in
ignored `tmp/latex/source/`; the repository does not relabel those external
packages as MIT. The builder first reuses the pinned REVTeX archive included in the committed
source ZIP, verifying its checksum. If that ZIP is absent, it downloads the
distribution. Tectonic downloads any missing standard TeX dependencies. Subsequent builds can use the populated caches.
`curl` is required for the initial REVTeX download.

The build rejects unresolved citations/references, missing glyphs, overfull boxes,
and bibliography-entry warnings. The known APS bibliography-style
control diagnostic is recorded separately in the manifest. It verifies US Letter page size, author/contact text,
and the absence of publication labels. Hashes bind the source, table inputs,
vector figures, PDF, and source archive. A fixed build epoch is recorded; exact
PDF bytes can still depend on the platform and external TeX dependency versions.
The source archive compiles without Python using its included tables, figures,
and unmodified REVTeX dependency archive.
See [SOURCE_README.txt](SOURCE_README.txt).

## Edit and inspect

Edit `manuscript.tex` and `references.bib`. Edit plot styling in `figures.py`; never
change the scientific inputs merely to improve a figure. `tables/*.tex` and
`figures/*.pdf` are generated files. LaTeX handles equation, section, figure,
table, and citation numbering automatically.

Render and inspect every page after a layout change:

```sh
pdftoppm -r 150 -png output/pdf/central-residence.pdf tmp/latex/page
pdffonts output/pdf/central-residence.pdf
```

Inspect figure labels at their final size, equation placement, float order,
column boundaries, references, and embedded fonts. Logs remain in `tmp/latex/`.
The previous report-style manuscript remains available in release 0.1.0; it is
not a second editable source for the current paper.
