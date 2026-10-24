Recoil transport and environmental exposure of black-hole merger remnants
Venkata Hemanth Athota
bharadwaj.hemanth.002@gmail.com

Manuscript revision 0.4.0-manuscript.1.
This archive contains a two-column theoretical manuscript, bibliography,
compiled bibliography, citation, license, and the unmodified REVTeX dependency.
It contains no calculated trajectories, model fits, or numerical result figures.

Reference build: Tectonic 0.16.9 and REVTeX 4.2f.
Unpack the dependency:
  unzip -jo vendor/revtex-4.2f.zip 'tex/latex/revtex/*' 'bibtex/bst/revtex/*'
Compile:
  tectonic --untrusted manuscript.tex

Alternatively, use a TeX installation with REVTeX 4.2 and txfonts:
  pdflatex manuscript.tex
  bibtex manuscript
  pdflatex manuscript.tex
  pdflatex manuscript.tex

Compilation typesets equations; it does not execute a physical calculation.
Standard TeX packages and fonts remain external dependencies. Tectonic can
download missing packages; a populated cache supports offline builds.
REVTeX controls typography, spacing, and margins. Different TeX distributions
can produce pagination differences.

The manuscript is an unreviewed formulation. No submission, acceptance,
selected cluster, or calibrated transport prediction is claimed.
CITATION.bib identifies the versioned source:
https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/manuscript-v0.4.0.1
No DOI has been assigned.

Project-authored source is MIT-licensed. The included REVTeX distribution and
external dependencies retain their own licenses.
