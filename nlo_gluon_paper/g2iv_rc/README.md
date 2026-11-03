# Group II Row IV: DGLAP, running coupling, finite terms

`input/G2IV.tex` (your equation) -> `make_note.py` -> `G2IV_DGLAP_RC.tex/.pdf`. The pole and mu bookkeeping is checked symbolically in `make_note.py`; the plus prescription and the fragmentation convolution in `../other_check/checks.py`.

## compact/
- `G2IV_compact_RC.tex/pdf`: check of the compact form of Group II Row IV, simplified DGLAP + running coupling + finite,
  the running coupling step by step (alpha_s(c0/|y'-y|)), and why line 2 (the y'=y subtraction) does not contribute
  to the running coupling (its mu-part multiplies delta^2(k)).
- `checks.py`: sympy checks (exact bracket with d_perp = 2 - eps, line 2 = -[log(vee/k+) - (1-k+/vee)] H_4, running coupling).
