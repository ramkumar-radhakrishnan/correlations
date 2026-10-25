# Check of the new Evolution_terms.tex (evolution terms, JIMWLK, BFKL, singlet forms)

- `Evolution_terms_check.pdf` (`.tex`): the check, equation by equation.
- `Evolution_terms_corrected.pdf` (`.tex`): your file with the fixes, removed text in red, new text in blue.
- `input/Evolution_terms.tex`: the file as sent.

## Result

- Sec. 1 (genuine), all 8 equations: agree with the checked genuine terms.
- Sec. 2, all 15 equations: identical, term by term, to the previous (checked) version.
- Sec. 3: the collected form (3.3) is JIMWLK x LO exactly, for any U; dN1 + dN2 equal it with U^{bc} = U^{cb}.
- Sec. 4: the collected form (4.2) is BFKL x LO exactly, for any U; (4.1) equals it with U^{bc} = U^{cb}.
- Singlet forms (3.4), (4.3): correct for any U.
- P1: the two-rho lists of 4rho.tex (30 terms, M1) are still missing in "From 4 rho's"; without them Secs. 1+2 differ from Sec. 3 (8-54 % pointwise). With them: exact.
- P2: Sec. 4 must carry log(vee/k+), not log(vee/Lambda); total = log(vee/Lambda) J + log(vee/k+) (B - J).
- P3: singlet forms need <rho^f(x) rho^f(y)>/(N_c^2-1) and the integral over y; K(a,b,z) to be defined.
- P4: the second "2b finite" list of Row V is the counterterm; "(with c.c.)" = c.c. already written out.

## Reproducing

```
python3 parse_v4.py        # parse all terms and trace forms
python3 check_all.py       # Tables 1-3 -> check_all.json
python3 check_singlet.py   # Table 4 -> check_singlet.json
python3 check_m1.py        # P1 -> check_m1.json, fourlists_body.tex
python3 m1_singlet.py      # P1 in the singlet projection -> m1_singlet.json
python3 make_corrected.py  # corrected copy, parsed back -> corrected.json
python3 make_check.py      # Evolution_terms_check.tex
pdflatex Evolution_terms_check.tex; pdflatex Evolution_terms_corrected.tex   # three times
```
