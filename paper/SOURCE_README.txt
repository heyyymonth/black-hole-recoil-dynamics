From recoil motion to companion access: controlled benchmarks and a cluster-snapshot comparison
Venkata Hemanth Athota

This archive contains the standalone LaTeX manuscript, bibliography, compiled
bibliography, five vector figures, and generated table fragments for manuscript
revision 0.2.0-manuscript.1, covering experiments E001 through E004.
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

To regenerate the scientific figures and tables, use the full research repository:
https://github.com/heyyymonth/black-hole-recoil-dynamics
Use the paper-specific citation to retrieve the matching immutable revision.

Project-authored source and generated artifacts are covered by the accompanying
MIT license. The included REVTeX distribution and external fonts retain their own licenses;
the project MIT license does not relicense them.
