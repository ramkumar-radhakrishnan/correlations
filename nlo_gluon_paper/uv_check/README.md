# Transverse divergences of Other_terms.tex

Input: `input/Other_terms.tex` (your upload; Group I Row II part 1 and Group II Row IV are commented out there).

* `terms.py` — splits every live equation into its additive terms and translates them for pointwise evaluation.
* `uvscan2.py` — UV scan: every pair of points P -> Q, angle-averaged integrand A(r) (principal value where a
  singular curve crosses the circle), exponent s in A ~ r^s; s <= -2 = divergent. Equations and whole rows.
* `surfscan.py` — singular curves (|x-z|=|x-w|, the D2/D3 curves of Group II Row I, the A_0 segment):
  integrable / principal value / divergent, per term, per equation, per row.
* `uvscan.py` — first version (shared evaluator: smooth symmetric SU(3) Wilson lines and charges).
* `uvtarget.py` — targeted pairs for the rows with p+ integrals (line by line); `intcheck.py` — Group I Row V
  d^3N_2b finite at x -> z with the p+, q+ integrals done.
* `make_report.py` -> `UV_divergences.tex/.pdf` — the result.
