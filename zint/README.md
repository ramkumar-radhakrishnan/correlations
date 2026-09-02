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

## Third note: the UV part on its own

`UV_part_only.pdf` isolates the z -> x divergence and nothing else, and answers the
factor-of-2 question about S_1 = (P.r)(s.r)/r^4:

    c_1 = d_perp xi xib = 2 xi xib      <S_1> = (s.P)/(d r^2) = (s.P)/(2 r^2)
    =>  c_1 <S_1> = xi xib (s.P)/r^2    -> the xi(1-xi) of Pgg, coefficient 1

No factor of 2 is missing. The 2 in c_1 and the 1/2 in the angular average are two
sides of the same d_perp = 2. Keeping d_perp symbolic in the vertex while averaging
in exactly 2 dimensions is the inconsistent hybrid: the correct statement is
d_perp/d = 1, so the xi(1-xi) piece gets no O(eps) correction. The remarks in the
first two notes that used d_perp/2 have been corrected.

* `src/uv_check.py` - the three independent re-derivations
* `src/run_uv.py`   - every number in the note -> `RESULTS_uv.txt`

## Fourth note: the UV-projected row, integrated

`UV_projected_row_integrated.pdf` takes the row after the UV projection,

    int d^2r e^{iq.r} (1/(xibar D)) (s.P)/(P^2 r^2) [xi xib + xi/xib + xib/xi]

and does both integrations. The whole r-dependence is one master,

    G = int d^2r e^{iq.r}/[(r^2+m^2)((r-s)^2+M^2)] = (pi/D0)[ln(D0/m^2) + Ihat]

with Ihat an explicit convergent 1-parameter Bessel integral. The 1/xibar of the
longitudinal denominator cancels against 1/D0 = xibar/s^2, leaving two WW kernels
times Pgg/(2Nc) times [log + Ihat]. The xi integral of the log is elementary; the
xi integral of Ihat is fixed at the endpoints (e^{ik.(y-x)} l_k^2/2 + c0 l_k + c1).

Result:  UV = (2 l_k - 11/6) log((y-x)^2/m^2)
         FIN = (1 + e^{ik.(y-x)})/2 * l_k^2 + c0 l_k + pi^2/6 - 67/36 + c1

* `src/uv_master.py` - the master G and Ihat, checked against 2D quadrature
* `src/run_uvint.py` - every number in the note -> `RESULTS_uvint.txt`

## The UV xi-integral at finite lambda

`src/xi_exact.py` gives the UV moment

    I(a, lam) = int_lam^{1-lam} dxi C_UV(xi) log(a/(1-xi)),   a = (y-x)^2/m^2

both exactly at finite lam = Lambda/k+ and in the lam -> 0 limit,

    I -> (2 l_k - 11/6) log a + l_k^2/2 + pi^2/6 - 67/36,   l_k = log(1/lam).

The exact form is I = A(lam) log a + B(lam) with

    A(lam) = 2 log(1/lam) + 2 log(1-lam) - 11/6 + 4 lam - lam^2 + (2/3) lam^3

and B the sum of three elementary antiderivatives (only Li_2 appears). The script
also implements the Mathematica ConditionalExpression form for comparison: the
three agree to 16 digits at every lam, so the lam -> 0 formula and the exact
finite-lam result are the same integral, differing only by terms that vanish with
the rapidity cutoff.

## Reduction of the Mathematica ConditionalExpression

`mathematica_reduction.pdf` does the algebra explicitly. With lam = Lambda/k,
l = log(1/lam), w = log(1-lam), A = log a, u = A - w, v = A + l, and

    P = 12 - 12 lam + 3 lam^2 - 2 lam^3,   Q = 1 + 12 lam - 3 lam^2 + 2 lam^3

six steps (all exact identities, checked in `src/mma_reduce.py`) give

    I = A(lam) log a + B(lam),
    A(lam) = 2l + 2w - 11/6 + 4 lam - lam^2 + (2/3) lam^3

the -11/6 being (Q-P)/6 at lam = 0; and B(lam) -> l^2/2 + pi^2/6 - 67/36, the
pi^2/6 coming from 36 Li2(1) = 6 pi^2. What the limit drops is exactly
lam [ 2 log a + l_k ] + O(lam^2 log lam).

* `src/mma_reduce.py` - the six steps as sympy identities

## The xi-integral run in a Wolfram kernel

`xi_integral_wolfram_run.pdf` is the verbatim session: `Integrate` of

    Log[a/(1-z)] (z(1-z) + z/(1-z) + (1-z)/z),  {z, Lambda/k, 1 - Lambda/k}

its raw output, the kernel's own reduction (collect in Log[a]) giving

    coefficient of Log[a] = -11/6 + 2l + 4lam - lam^2 + (2/3)lam^3 + 2w

and the lam -> 0 limit plus the first correction, both produced by the kernel:

    LIMIT      = (2 lk - 11/6) Log a + lk^2/2 + Pi^2/6 - 67/36     (True)
    CORRECTION = lam (2 Log a + lk)                                (True)

Two Mathematica traps are documented there: Series will not expand
PolyLog[2, 1-lam] (branch point at 1 -- substitute the inversion identity), and
Simplify will not split Log[a/lam] without sign assumptions, so an === 0 test can
report False on two equal expressions.

## Dimensional regularisation, MS-bar (in momentum space) -- both pieces

`dimreg_msbar.pdf`. The exact split is  1/(r^2 D) = 1/(r^2 D0) + (1/r^2)[1/D - 1/D0],
i.e.  G = T/D0 + R.

**T** (the divergent piece). 1/r^2 is the square of the WW kernel, so by the
convolution theorem it is the one-loop transverse bubble,

    int d^2r e^(iqr)/r^2 = (2pi)^2 int d^2l/(2pi)^2  l.(l-q)/(l^2 (l-q)^2)

whose divergence is the large-l (UV) region. In d = 2-2eps with the standard
measure this is pi[1/epshat + log(mu^2/q^2)], 1/epshat = 1/eps - gamma + log 4pi,
so MS-bar leaves pi log(mu^2/q^2). Keeping eps_UV and eps_IR apart shows the pole
is genuinely +1/eps_UV, the IR poles cancelling between the tadpole and the bubble.

**R** (finite, no regulator). Using 1/D - 1/D0 = (2 r.s - r^2)/(D D0),

    R = (2 s.Psi - Phi)/D0,   Phi = 2 pi e^(iq.s) K0(M|q|),
    s.Psi = pi int_0^1 da e^(i(1-a)q.s)[(1-a)s^2 |q| K1(lam)/sqrt(Del) + i(q.s)K0(lam)]

**Adding them the q^2 cancels**, and the logarithm becomes the *dipole size*:

    G|MSbar = (pi/D0)[ log(mu^2 (y-x)^2 / xibar) + 2(gamma - log 2) + Ihat ]

That 2(gamma - log 2), riding on (2 l_k - 11/6) x 2Nc = 4 Nc l_k - b, is exactly
the 2 b (gamma - log 2) of K_JSJ, Eq. (2.64) of arXiv:1610.03453 -- the term their
footnote 2 records as missing from their earlier results. It appears only because
R is kept.

Final:
    P = (1/epshat)(2 l_k - 11/6)                          -> 0 in MS-bar
    F = (2 l_k - 11/6)[log(mu^2 (y-x)^2) + 2(gamma - log2)]
        + (1 + e^{i k.(y-x)})/2 * l_k^2 + pi^2/6 - 67/36 + c0 l_k + c1

prefactor alpha_s^2 Nc/(8 pi^5 k+). The l_k^2 coefficient is regulator-independent
and matches the mass-regulator computation.

* `src/remainder.py`  - R, s.Psi, Phi; checked against 2D quadrature to 1e-9
* `src/run_dimreg.py` - all the numbers -> `RESULTS_dimreg.txt`
