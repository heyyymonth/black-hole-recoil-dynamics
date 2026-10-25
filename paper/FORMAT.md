# Manuscript formatting conventions

The selected target is an **authored, two-column physics manuscript**, using
REVTeX 4.2f with APS/PRD options and `txfonts`. The class supplies the page,
column, heading, equation, and reference conventions. Formatting sources were
checked on 18 October 2026.

| Element | Implementation |
|---|---|
| Page and columns | REVTeX APS/PRD `reprint` defaults, US Letter, two columns; `\raggedbottom` preserves natural paragraph spacing; no custom geometry or compressed spacing |
| Body and mathematics | 10-point class typography with Times-like `txfonts`; native vector mathematical notation |
| Title and authors | Centered title; Venkata Hemanth Athota; genuine independent-research affiliation; email footnote |
| Abstract | One self-contained paragraph, no numbered citations or displayed equations |
| Headings | Class-generated Roman section numbers and lettered subsections |
| Equations | Native numbered `equation` and `align`, with automatic cross-references |
| Citations | Numeric bracketed citations in order of first occurrence; BibTeX and `apsrev4-2` |
| Figures, when present | Vector PDF, approximately 3.4-inch single-column or 7-inch full-width; 8–10-point figure text at final size, line widths at least 0.5 pt, panel labels and grayscale-distinguishable line styles |
| Tables, when present | Native LaTeX with captions defining quantities and units; no colored fills or vertical rules |
| Links and markings | Black unobtrusive links; no decorative header, watermark, arXiv stamp, invented DOI, or acceptance notice |

The article includes three numerical evidence tables and three vector figures from the observed-cluster data analysis. Their captions distinguish measured profiles, published model estimates, and conditional comparison scales.

The `reprint` option names the compact layout. It does not assert that the work
has been published. The class controls dimensions and typography; a document
should not be squeezed manually to satisfy a venue's page limit.

This is a formatting baseline, not a declaration of eligibility for every
conference or journal. Actual submission requirements also depend on article
type, scope, page limits, anonymity, data availability, and disclosures. In
particular, the PRD instructions require appropriate disclosure of substantive
AI use before submission. No submission has been made as part of this formatting
change.

ICLR uses a different, single-column, double-blind template. The second example,
[*Attention Is All You Need*](https://arxiv.org/pdf/1706.03762), is a NeurIPS-family
layout, not evidence that REVTeX satisfies ICLR. If the target venue changes,
apply its official template rather than mixing incompatible layout rules.

## Primary formatting references

- [APS REVTeX distribution and guidance](https://journals.aps.org/revtex)
- [REVTeX FAQ: preprint and reprint layouts](https://journals.aps.org/revtex/revtex-faq)
- [APS style basics: figures, notation, and references](https://journals.aps.org/authors/style-basics)
- [Physical Review D author instructions](https://journals.aps.org/prd/authors)
- [ICLR 2027 author guidelines](https://iclr.cc/Conferences/2027/AuthorGuidelines)

Project-authored files retain the MIT license. REVTeX and font packages are
external dependencies with their own licenses; see [REVTeX on CTAN](https://ctan.org/pkg/revtex)
and [txfonts on CTAN](https://ctan.org/pkg/txfonts). No third-party paper content
or figures are included in the research release.
