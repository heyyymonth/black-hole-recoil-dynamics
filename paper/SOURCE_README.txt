Recoil transport and environmental exposure of black-hole merger remnants
Venkata Hemanth Athota
bharadwaj.hemanth.002@gmail.com

Manuscript revision 0.8.0-manuscript.2.
This archive contains a two-column manuscript on recoil transport and
environmental exposure, including the six-cluster data analysis, a fixed NGC 3201 host
reconstruction, local anisotropic collision calculations and static
external-field boundary diagnostics, and Galactic reference orbits with
time-dependent host-support checks.
It includes the article, seven table fragments, seven vector figures,
bibliography, compiled bibliography, citation, license, supporting data
evidence, and the unmodified REVTeX dependency.

The numerical findings describe existing observations, unit conversions,
published joint host posteriors, sampling coverage and model sensitivity.
A checked host, twenty local collision states and four static boundary
scenarios and four Galactic reference orbits are included. No new host
fit, calibrated local scattering background, remnant trajectory, or residence
calculation is supplied.

Reference build: Tectonic 0.16.9 and REVTeX 4.2f.
Unpack the dependency from this directory:
  unzip -jo vendor/revtex-4.2f.zip 'tex/latex/revtex/*' 'bibtex/bst/revtex/*'
Compile, preserving the tables/ and figures/ directories:
  tectonic --untrusted manuscript.tex

Alternatively, use a TeX installation with REVTeX 4.2 and txfonts:
  pdflatex manuscript.tex
  bibtex manuscript
  pdflatex manuscript.tex
  pdflatex manuscript.tex

Compilation typesets the existing equations, tables and figures. It does
not retrieve observations, rerun the data analysis, generate plots, fit a
host or execute transport simulations. Standard TeX packages and fonts
remain external dependencies. Tectonic can download missing packages;
a populated cache supports offline builds. Different TeX distributions
can produce pagination differences.

The evidence/ directory contains the paper asset manifest, source manifest,
candidate host records and the derived tables and diagnostics used by the
article. The observed-cluster audit refers to paths in the full repository;
not every audited repository file is duplicated here. Original large
third-party source files are not bundled. Source URLs, exact source hashes,
attribution, and instructions for the full data audit are retained in
evidence/observed-clusters/. The standalone TeX build uses only the supplied
article, bibliography, tables and figures.

The manuscript is unreviewed. No submission, acceptance, or calibrated
transport prediction is claimed. CITATION.bib identifies the versioned source:
https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/manuscript-v0.8.0.2
No DOI has been assigned.

Project-authored source is MIT-licensed. The included REVTeX distribution,
external dependencies, and third-party source material retain their own
licenses. Preserve the source attribution in the supporting data README;
the repository license does not relicense external data.
