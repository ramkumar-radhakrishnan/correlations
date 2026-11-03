# Finite terms: final, paper-ready

Input: `input/Other_terms.tex` (your upload: all terms after JIMWLK and BFKL, exact in vee).

* `Final_result.tex/.pdf`: check of the file (25/30 equations identical to verified ones; the new Group I Row II
  part 1 and Group II Row IV checked; edits E1-E4), poles -> DGLAP + running coupling (+ the 2K constant, Q1, V_vee),
  leading power in a = k+/vee, simplified finite terms, Gaussian projectile in colour-singlet form, final result,
  AABKL comparison, exact-a dictionary, open points, checks.
* `make_final.py`: builds Final_result.tex from `front.tex`, `back.tex` and the verified sections of
  `../limiting/Limiting_terms.tex` and `../limiting/Gaussian_singlet.tex` (labels prefixed LT:, GS:).
* `checks.py` -> `checks.json`: equation identity, 2b UV -> Row II part 1 coefficients, leftover pole,
  Group III plus prescription with exact a, closed form, prefactors, a -> 0 convergence.

```
python3 checks.py; python3 make_final.py; pdflatex Final_result.tex   # three times
```
