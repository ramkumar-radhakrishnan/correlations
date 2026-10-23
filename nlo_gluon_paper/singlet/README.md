# Evolution terms in colour-singlet form

`Evolution_singlet.pdf` (`.tex`): every two-rho log(vee/Lambda) evolution term with the projectile in a colour
singlet, <rho^a(x) rho^b(y)> = delta^{ab} <rho^f(x) rho^f(y)>/(N_c^2-1), written as adjoint traces.

* `trace.py` - projection of a term and conversion of its colour network to traces of U, U^dagger, T^e.
* `validate.py` - trace form vs direct contraction for every term (general SU(3) U) -> `validate.json`.
* `render_singlet.py` - collected singlet forms per source, JIMWLK x LO, BFKL x LO -> `singlet.pkl`, `singlet_checks.json`.
* `compare_kernels.py` - sum of the lists vs JIMWLK x LO per kernel structure -> `compare_kernels.json`.
* `compare_rows.py`, `compare_rows_fast.py` - your rows vs the general-U reference rows -> `compare_rows*.json`,
  `rows_collected.pkl`; `ref_rows.pkl` holds the reference pieces I(a), I(b), I', II(a), II(b), ...
* `chains.py`, `orient.py` - reading rules for the index order with a general U -> `orient.json`.
* `render_tex.py`, `make_singlet.py`, `singlet_static.tex` - the note.

Run order: `validate.py`, `render_singlet.py`, `compare_kernels.py`, `compare_rows.py`, `compare_rows_fast.py`,
`orient.py`, `make_singlet.py`, then pdflatex three times.
