# Negative logarithms in Other_terms.tex

Input: `input/Other_terms.tex` (your upload).
* `make_corrected.py` -> `corrected/Other_terms.tex` (clean) and `corrected/Other_terms_marked.tex` (red brackets removed,
  blue bars new): the 15 logarithms of Group I Row III whose argument changes sign become log|.| (principal value).
* `pv_check.py` -> `pv_check.json`: each corrected logarithm = principal-value p+ integral, in all regions (4e-12).
* `Negative_logs.tex/.pdf`: the note.

## modulus/
* `Why_modulus.tex/.pdf`: step-by-step derivation of the log|.| (principal-value p+ integral), with `check.py`.

## iepsilon/
* `Negative_logs_iepsilon.tex/.pdf`: the i-epsilon treatment (advisor's request). Energy denominators E_g(k)-E_g(p) get +i eps,
  (i3) -> 1/(xi Z^2 - X^2 + i0), p+ integrals -> log|.| - i pi theta; in term + c.c. the i pi multiplies (O + O^dagger),
  whose Gaussian average vanishes, so log|.| is the final answer (C-odd d^{abc} part would survive).
* `check.py` -> `check.json`: i eps integrals vs log|.| - i pi theta; reversed-order colour traces (SU(2): 0, SU(3): d/2).
