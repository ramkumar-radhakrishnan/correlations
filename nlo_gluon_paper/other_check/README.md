# Other_terms.tex: mu coefficients, DGLAP, finite terms, running coupling

Input: `input/Other_terms.tex` (your remaining terms, after JIMWLK and BFKL).

* `Other_terms_analysis.tex/.pdf` — the note:
  1. every coefficient of mu (table, with the coupling and fragmentation counterterms);
  2. DGLAP from Group II Row IV, step by step, and the remaining terms of that row;
  3. the other finite terms of the cross section (exact in vee, Lambda -> 0 only where the limit exists term by term);
  4. running coupling.
* `diff_check.py` -> `diff_check.txt`: your file against the checked remaining terms (`../full_check/Remaining_RC_DGLAP_finite.tex`), equation by equation.
* `checks.py` -> `checks.json`: mu table, prefactors, plus prescriptions, Fourier identity, fragmentation convolution.
* `make_note.py`: writes the note; your equations are copied verbatim (Lambda -> 0 edits asserted, red/blue).

```
python3 diff_check.py; python3 checks.py; python3 make_note.py
pdflatex Other_terms_analysis.tex; pdflatex Other_terms_analysis.tex
```
