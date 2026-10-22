Companion access after black-hole recoil: dynamical benchmarks and sparse binary populations
Venkata Hemanth Athota

This archive contains the standalone LaTeX manuscript, bibliography, compiled
bibliography, six vector figures, generated table fragments, and CITATION.bib for manuscript
revision 0.3.0-manuscript.4, covering orbital exposure, velocity-weighted passages,
a cluster-snapshot comparison, and the late-time binary-target census.
This manuscript revision is separate from the baseline software release 0.1.1.
No journal publication or acceptance is claimed by the document class.

Reference build: Tectonic 0.16.9 and REVTeX 4.2f
First unpack the unmodified official REVTeX dependency into the working directory:
  unzip -jo vendor/revtex-4.2f.zip 'tex/latex/revtex/*' 'bibtex/bst/revtex/*'
Then compile:
  tectonic --untrusted manuscript.tex

Alternative: a current TeX Live or MiKTeX with REVTeX 4.2 and txfonts
  pdflatex manuscript.tex
  bibtex manuscript
  pdflatex manuscript.tex
  pdflatex manuscript.tex

REVTeX controls the two-column page dimensions and typography. Do not change
margins, line spacing, or font sizes to meet a page limit. All files needed for
article content are included. The official REVTeX distribution is included unmodified with its own sources,
documentation and license notices. Other standard TeX packages remain external
dependencies.
Initial Tectonic compilation may download these packages; a populated cache can
be used offline. Different TeX distributions can produce pagination differences.

The included CITATION.bib identifies this manuscript revision and its versioned
source tag, manuscript-v0.3.0.4. To regenerate the scientific figures and tables,
retrieve the matching full research repository at:
https://github.com/heyyymonth/black-hole-recoil-dynamics/tree/manuscript-v0.3.0.4
No DOI has been assigned.

Project-authored source and generated artifacts are covered by the accompanying
MIT license. The included REVTeX distribution and external fonts retain their own licenses;
the project MIT license does not relicense them.
