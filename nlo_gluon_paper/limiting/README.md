# Finite terms in the limit Lambda -> 0, vee -> infinity

Input: `input/Limiting_cases.tex` (your upload). Reference: `../neglog/corrected/Other_terms.tex` (checked exact-vee terms).

* `Limiting_terms.tex/.pdf`: check of the limits, the three edits (E1-E3), finiteness, simplification
  (c.c. = 2 Re; one xi integral and Ein closed form for G1 V(1); G1 V 2a + 2b finite combined; Group III in
  zeta = k+/(p+ + k+); closed-form p+ integral of G3 III(2); CMW constant; G2 IV = DGLAP + RC + finite),
  final result, and comparison with AABKL (arXiv:1411.2869).
* `limcheck.py` -> `limcheck.json`: every equation vs the exact-vee one at a = k+/vee = 1e-2 ... 1e-8 (difference ~ a).
* `largep.py` -> `largep.json`: Group III large-p+ behaviour (Fourier transform 1/q unsubtracted, 1/q^2 subtracted).
* `checks.py` -> `checks.json`: analytic limits, prefactors, zeta measures, phases, Ein closed forms, q+ -> xi map,
  2a colour at z = x.
* `limblocks.py`: splits the TeX files into equations.

```
python3 limcheck.py; python3 largep.py; python3 checks.py
pdflatex Limiting_terms.tex; pdflatex Limiting_terms.tex
```

## Gaussian projectile, colour-singlet form
* `Gaussian_singlet.tex/.pdf`: three-rho rows -> two-rho parts (pair rule; symmetric part = 0 for a Gaussian),
  every finite term projected on the colour singlet (adjoint traces), final result.
* `singlet3.py` -> `singlet3.pkl`, `singlet3.json`: reduction + projection (machinery of ../finite_singlet),
  including the new y' -> z subtraction of Group III; checks (traces vs direct contraction, pair rule vs engine,
  z-substituted operators).
* `srender3.py` -> `singlet3_frag.tex`: the traces used in the note (generated).
