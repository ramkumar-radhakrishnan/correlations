# Finite terms: complete simplification, running coupling and DGLAP extracted

`Finite_terms_simplified.pdf` (`.tex`): the finite terms of `input/Finite_terms_1.tex` (identical to
`../finite/input/Finite_terms.tex`), simplified completely, exact in a = k+/vee (only Lambda -> 0).
Replaces `../finite/Finite_terms_exact.pdf`.

New compared with the earlier note:
- consistent with the corrected evolution terms (1/2 log(vee/Lambda) (J + B)): every row = C_row log(vee/k+) + R_row(a);
  C_row log(vee/k+) -> C_row L*, L* = log(sqrt(vee Lambda)/k+); sum C_row L* = log(k+/sqrt(vee Lambda)) (J - B);
- c.c. checked: Group II rows are self-conjugate (no c.c.), Group I and III rows get + c.c.;
- rows with the same kernel combined (K_IV family: G1 IV, G2 I(1), G3 I, II; K_C family: G3 III(1), V, VI),
  their colour factors collapse to compact forms;
- every three-rho row = symmetric three-rho part + two-rho part (written out, merge pattern by merge pattern);
- Group II Row I part 2 collected.

## Scripts

| File | Role |
|---|---|
| `split.py` | C_row / R_row split per row (R_row constant for a = 1e-4..1e-8); c.c. self-conjugacy test -> `split.json` |
| `decomp.py` | three-rho colour operators through the reordering engine -> `decomp.pkl`, `decomp.json` |
| `compact.py` | two-rho parts by the pair rule, combined and simplified, checked against the engine -> `compact*.pkl`, `compact*.json` |
| `final_checks.py` | closed forms vs file, G2 I(2) collected, compact colour forms, sympy tallies -> `final_checks.json` |
| `render.py`, `make_note.py`, `static.tex` | the note |

```
python3 split.py; python3 decomp.py; python3 compact.py
python3 -c "import env, compact; compact.run_simplify()"
python3 final_checks.py; python3 make_note.py; pdflatex Finite_terms_simplified.tex   # three times
```
