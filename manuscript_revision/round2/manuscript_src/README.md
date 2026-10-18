Source of `../manuscript_round2_orange.pdf`, reconstructed from the authors' PDF (CSL_notes_1.pdf) because the
current .tex was not available:

- `manuscript.tex`: original submission + H→b rename + the authors' blue changes (verbatim from their PDF)
  + the new orange changes (`\textcolor{orange}`). A word-by-word comparison with the authors' PDF shows
  identical text apart from the orange additions and the two new references.
- `refs_rebuilt.tex`: the reference list rebuilt from the authors' PDF (same numbering), plus Siegel et al. (2013)
  and Guilet et al. (2017). It replaces `\bibliography{references}`; in the real source use bibtex instead.
- `fig*.pdf`: vector crops of the figures in the authors' PDF.

Compiled with the Springer Nature `sn-jnl.cls` (not included; use your own copy):
`pdflatex manuscript` (three times).
