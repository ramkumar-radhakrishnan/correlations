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

## The four-kernel row (two WW kernels at each of z and y')

`four_kernel_row.pdf`. New cross section: four Weizsacker-Williams kernels

    (z-x)^m/(z-x)^2 * (x'-z)^{k'}/(x'-z)^2 * (y'-x')^i/(y'-x')^2 * (y-y')^j/(y-y')^2

contracted through  [ d_{k'm} d_{ij}/(p+ + k+) - d_{jm} d_{i k'}/p+ - d_{im} d_{j k'}/k+ ].

**A. the contraction collapses.** With xi = p+/k+ the bracket gives exactly

    (1/k+)[ (K1.K2)(K3.K4)/(1+xi) - (K1.K4)(K2.K3)/xi - (K1.K3)(K2.K4) ]

three products of two scalar dot-products of kernels (verified symbolically).

**B. there is NO transverse UV divergence.** Each kernel is 1/|sep| and no pair of
points appears in two kernels, so nowhere is there a 1/sep^2; every coincidence
limit is integrable in 2D. This is the structural difference from the earlier row,
where the light-cone vertex supplied a second 1/(x-z)^2.

**C. one master does both integrations.**

    B^{mj}(Q,c) = int d^2r e^{-iQ.r} (r^m/r^2) ((c-r)^j/(c-r)^2)

  int d^2y' : Q = (1+xi)k, c = x'-z ;  int d^2z : Q = -xi k, c = x-y'.
  Checked by delta_{mj} B^{mj} = -conj(A7) to 1.2e-15.

**D. exact three-way separation.**

    B^{mj} = -(pi/2) d^{mj} log(4 e^{-2gamma}/(Q^2 c^2))      logarithm
             - pi d^{mj} + pi Q^mQ^j/Q^2 + pi c^mc^j/c^2      finite constants
             + Bcal^{mj}(Q,c)                                 remainder, O(|Q||c|)

The "UV-like" log here is a *collinear/large-separation* log, not a transverse UV
one: it is cut by the dipole size, and it is finite for any physical Q, c.

**E. the p+ integral.** For the z-integration the log goes as log(1/xi^2) and hits
the 1/xi pole of the second bracket term, giving a double rapidity log,
int_lambda^Xi dxi/xi log(1/xi^2) = l_k^2 - log^2 Xi. The -d_{im}d_{jk'}/k+ term
carries no 1/xi and gives int dxi = Xi - lambda, linear in the upper cutoff.

* `src/nlo2_setup.py`  - symbolic contraction of the bracket
* `src/nlo2_master.py` - B^{mj}(Q,c) in the alpha representation
* `src/run_nlo2.py`    - all checks -> `RESULTS_nlo2.txt`

### Addendum: where the four-kernel row actually diverges

`four_kernel_row_divergences.pdf`. Complete enumeration of the singular regions of
int d^2z d^2y', answering "does it have any transverse divergence?".

**Answer: yes, exactly one, and it is infrared, not ultraviolet.**

- r -> 0 and r -> c (all four coincidence limits): integrand ~ 1/rho with zero
  angular average of the leading piece; rho^2 <integrand> falls like rho^2 over
  three decades. int rho drho / rho finite. NO UV divergence.
- r -> infinity: (c-r)^j/(c-r)^2 -> -r^j/r^2, so the integrand -> -r^m r^j/r^4,
  angular average -delta^{mj}/(2 r^2). LOG DIVERGENT, coefficient -pi delta^{mj},
  PURE TRACE (traceless part averages to zero). Measured d/dlogR = -6.283185
  against -2pi = -6.283185.
  Algebraically: r.(c-r)/(r^2(c-r)^2) = -1/2[1/r^2 + 1/(c-r)^2 - c^2/(r^2(c-r)^2)].
- the regulator is the phase: trace B = -2pi log(2 e^-gamma/(|Q||c|)), i.e.
  rho_min = |c|, rho_max = 2 e^-gamma/|Q|. Divergent only as |Q| -> 0.
- |Q| -> 0 happens in ONE place: the z-integration, Q = -xi k, at xi -> 0. The
  y'-integration has Q = (1+xi)k >= k and is finite for every xi.
- c -> 0 (the two singular points colliding) is a region of the *other* integral
  (z -> x' resp. y' -> x); log^n |c| against rho drho converges. Harmless.
- Consequence: the log multiplies the 1/xi pole of the -delta_{jm}delta_{ik'}/p+
  term (whose z-integral is exactly the pure trace) -> double log l_k^2. The
  divergence belongs to the rapidity evolution of the Wilson-line correlators,
  not to a transverse counterterm.

* `src/ir_check.py` - all of the above -> `RESULTS_nlo2_ir.txt`

## After the p+ integration: is it right, and where is the rapidity log

`pplus_integrated_row.pdf`. With kap = k.(y'-z), xi = p+/k+, lam = Lambda/k+,
Xi = (V-k+)/k+, l_k = log(k+/Lambda):

**All three p+ integrals check out** (against quadrature, to 1e-29 relative):

    -C          : int dxi e^{-i kap xi}          = [e^{-i kap Xi} - e^{-i kap lam}]/(i kap)
    -B/xi       : -int dxi e^{-i kap xi}/xi      = -[Ei(-i kap Xi) - Ei(-i kap lam)]
    +A/(1+xi)   : int dxi e^{-i kap xi}/(1+xi)   = e^{i kap}[Ei(-i kap(1+Xi)) - Ei(-i kap(1+lam))]

**One flag: the overall 1/k+.** The first line's denominator is i k.(y'-z) with no
k+, which is the fingerprint that the k+ of dp+ = k+ dxi was already spent. So
(1/k+) int dp+ = int dxi leaves NO 1/k+. The 1/k+ out front is legitimate only if
it predates the bracket (e.g. from the gluon phase space).

**Limits.** Ei(-i kap X) -> -i pi sgn(kap) as X -> infinity (a constant, not zero).
Ei(-i kap eps) = gamma + log(i kap eps). Hence, as Lambda -> 0 and V -> infinity:

    A/(1+xi) term -> e^{i kap}[-i pi sgn(kap) - Ei(-i kap)]      no l_k
    -C term       -> i/kap + oscillatory                          no l_k
    -B/xi term    -> -l_k + gamma + log(i k.(y'-z))               THE rapidity log

So there is exactly ONE rapidity logarithm, coefficient -1, on structure B.

**What multiplies l_k is the BK kernel.** B contains
(x-z).(y'-z)/[(x-z)^2 (y'-z)^2] = K(x,y';z), the real-emission piece of the
BK/JIMWLK dipole kernel. This also explains the transverse IR log of the previous
note: the full kernel M = (x-y')^2/[(x-z)^2(y'-z)^2] = 1/(x-z)^2 + 1/(y'-z)^2 - 2K
falls as rho^-4, while -2K alone falls as rho^-2. Measured angular averages at
large rho: 1.000000 + 1.000000 - 2.000000 = 0.000000. The IR log is telling you the
VIRTUAL diagrams are still missing, not that anything is wrong.

* `src/pplus_check.py`, `src/pplus_limits.py` -> `RESULTS_pplus.txt`

## The same answer in transverse position space

`position_space_msbar.pdf`. Only one thing in the MS-bar result was ever in
momentum space: the argument of the log, pi log(mu^2/q^2) with q = xibar k. Turning
it into log(mu^2 (y-x)^2) is NOT algebra -- it needs the remainder
R = int d^2r e^{iqr}(1/r^2)[1/D - 1/D0], which is UV finite (easy to drop) but
carries an opposite-sign log q^2.

**Regulate in position space instead.** Cut at |r| > r0:

    int_{|r|>r0} d^2r e^{iqr}/r^2 = 2 pi int_r0^inf dr J0(qr)/r
                                  = pi log(4 e^-2gamma /(q^2 r0^2)) + O(q^2 r0^2)

which is the SAME expression as the mass regulator 2 pi K0(m|q|) with m = r0
(checked by quadrature). MS-bar dictionary:  r0 = 2 e^-gamma / mu.

**The master, entirely in position space** (verified against 2D quadrature of the
cut integral to 1e-8):

    G = int d^2r e^{iqr}/(r^2 D) = (pi/D0)[ log(D0/r0^2) + Ihat(xi) ],
    D0 = (y-x)^2 / xibar

No gamma, no log 2, no q^2. Substituting r0 = 2 e^-gamma/mu gives
log(D0/r0^2) = log(mu^2 (y-x)^2/xibar) + 2(gamma - log 2) exactly -- so the
2(gamma - log 2) is nothing but the MS-bar translation constant.

**Final:**

    F = (2 l_k - 11/6) log((y-x)^2/r0^2)
        + (1 + e^{i k.(y-x)})/2 * l_k^2 + pi^2/6 - 67/36 + c0 l_k + c1

with the moments int C_UV = 2 l_k - 11/6, int C_UV log(1/xibar) =
l_k^2/2 + pi^2/6 - 67/36, int C_UV Ihat = e^{ik.(y-x)} l_k^2/2 + c0 l_k + c1.

* `src/posspace.py` -> `RESULTS_posspace.txt`

## Two dipole kernels integrated over the shared emission point

`two_kernel_z_integral.pdf`.

    int d^2z  (z-w').(z-w)/[(z-w')^2 (z-w)^2] * (x-z).(y-z)/[(x-z)^2 (y-z)^2]

**It converges; there is nothing to regulate.** A dipole kernel is only O(1/rho) at
its endpoints (the numerator vanishes linearly), and no point is an endpoint of both
kernels, so at each of the four singular points only one kernel is singular and the
integrand is O(1/rho): int rho drho/rho converges. At large |z| each kernel goes as
1/z^2, so the product goes as 1/z^4 (measured rho^4 <integrand> -> 1.000000).

**Closed form**, with K(u,v) = u.v/(u^2 v^2):

    = -(pi/2) [ K(w'-x, w-y) log( (w-w')^2 (x-y)^2 / [(w'-y)^2 (w-x)^2] )
              + K(w'-y, w-x) log( (w-w')^2 (x-y)^2 / [(w'-x)^2 (w-y)^2] ) ]

The log arguments are conformal cross-ratios. One term per pairing of {w',w} against
{x,y}; the log's denominator uses the other pairing.

**Method.** In complex coordinates K(a-z, b-z) = Re 1/[(zbar-abar)(z-b)], so the
product is (1/2) Re[(a)+(b)] with each term the master

    M(a1,a2;b1,b2) = int d^2z / [(zbar-a1b)(zbar-a2b)(z-b1)(z-b2)]
      = -pi/[(a1b-a2b)(b1-b2)] log[ |a1-b1|^2 |a2-b2|^2 / (|a1-b2|^2 |a2-b1|^2) ]

built from L(a,b) = int_{|z|<R} d^2z/[(zbar-abar)(z-b)] = pi log(R^2/|a-b|^2). The
four partial-fraction terms carry signs +,-,-,+ so log R^2 cancels.

Checked: all four symmetries exact to 12 digits; against direct 2D quadrature
(partition of unity about the four singular points) to 1e-7..1e-5 on four random
configurations.

**Where a divergence does appear:** in the external points, when any two collide
(w->w', x->y degenerate one kernel; w'->x etc. put two kernels on the same point).
All six are logs against a 2D measure, integrable in the remaining int_{x,y,w,w'}.
The only divergence left in the cross section is the rapidity log log(V/Lambda).

* `src/twokernel.py`, `src/tk_fast.py` -> `RESULTS_twokernel.txt`
