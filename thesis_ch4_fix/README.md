# Thesis corrections (Chapters 3–4): response to the examiner's report

`thesis_corrections_note.pdf` is the correction note, and `thesis_corrections_note.tex` is its source. The note goes through the examiner's seven points one by one, lists additional errors found, gives corrected equations, recomputes Tables 4.1–4.4 by exact diagonalisation, and ends with a checklist of changes to the thesis.

| File | Purpose |
|---|---|
| `derive.py` → `derive.out` | sympy derivation of the first-order states (quartic and symmetric well, n = 0–3) |
| `verify.py` → `verify.out` | parity, Fourier-transform and normalisation checks; fidelity of first-order states against exact eigenstates |
| `compute.py` → `compute.out`, `results.json` | Chapter 3 measures (Wigner and Husimi), exact Tables 4.1–4.4, many-state S_x + S_k |
| `ch4_curves.py` → `ch4_curves.json`, `figE_ch4_exact.pdf` | corrected Figs. 4.7–4.10 |
| `figs.py` → `figA`–`figD` | corrected Figs. 3.2 and 3.4, Chapter 3 measures, validity and many-state plot |
| `check_note.py` → `check_note.out` | independent check of the formulas as printed in the note |

Requirements: Python 3 with NumPy, SciPy, SymPy and Matplotlib. Run `compute.py` before `figs.py`, because `figs.py` reads `ps_n*.npz` and `results.json`.
