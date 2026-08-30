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

## Second note: the p+ integration and Lublinsky-Mulian section 4.2

`xi_integration_to_sec42.pdf` takes the same row and does the **xi (p+) integration
first**, in the style of arXiv:1610.03453 section 4.2 / appendix H.2.

Dictionary found: with `z_paper <-> x`, `z'_paper <-> z`, `y <-> y`,

    xibar * D  =  xibar (Y')^2 + xi Y^2

is exactly the energy denominator of Eq. (4.16), and the six xi-weights of the
bracket, (xi*xib, xi/xib, xib, xib/xi, xi, 1), are exactly the products the paper's
Lambda/Theta vertex decomposition produces. The p+ integral then needs only the
paper's (C.2), (C.5), (C.6) and gives

    G(r) = l_k * G_rap  +  ln(Y^2/(Y')^2) * G_log  +  G_rat

with G_rap = V_4/(Y')^2 + V_2/Y^2 -- the (H.7) pattern -- and the z->x (UV) content
sitting in one place, coefficient  2 l_k - 11/6 = int dxi Pgg/(2Nc).

* `src/xi_integration.py` - symbolic check of the vertex decomposition
* `src/xi_masters.py`     - the six elementary xi-integrals J_1..J_6
* `src/xi_assemble.py`    - the kernel G(r), its three blocks, the UV separation
* `src/run_xi.py`         - every number in the note -> `RESULTS_xi.txt`
