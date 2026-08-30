# z-integration of the C1 x B2 row: UV, finite and regular terms

`z_integration_UV_finite_regular.pdf` is the note. It performs the transverse
z-integral of the bracket of Eq. (1.1) in full and separates it into

  (I)   the log(1/m^2) that comes from z -> x,
  (II)  the finite companion of that log,
  (III) the eight structures that were never divergent.

The central algebraic statement is exact and pointwise in r = z - x:

  T(r) = C_UV(xi) (s.P)/r^2 + T_fin(r),   C_UV = xi*xib + xi/xib + xib/xi = Pgg/(2Nc)

## Code

* `src/algebra.py`  - sympy: builds T(r) by contracting Eq. (1.1) directly, checks it
  against the hand-reduced Eqs. (1.4)-(1.10) (Eq. (1.8) has a factor-2 typo), against
  the nine-structure basis, and against the UV split.
* `src/masters.py`  - alpha (Feynman) representations of every master integral.
* `src/quad2d.py`   - direct 2D oscillatory quadrature used to check them.
* `src/assemble.py` - assembly of the three groups.
* `src/run_all.py`  - every number quoted in the note -> `RESULTS.txt`.

    pip install numpy scipy sympy mpmath
    python3 src/run_all.py

## PDF

    cd pdf && npm install katex && pip install playwright && python3 build_pdf.py
