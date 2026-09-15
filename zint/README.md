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

**The overall 1/k+ is correct** (an earlier draft of this README flagged it as
doubled -- that flag was wrong). The w' delta sets w' = y' + (p+/k+)(y'-z), so
w'-z = (p+ + k+)/k+ (y'-z) and the first kernel supplies k+/(p+ + k+); against the
1/k+^2 in front, and pulling (p+ + k+) out of the bracket, this reproduces exactly
(1/k+)[ d/(p+ + k+) - d/p+ - d/k+ ]. And the measure does not eat it: every bracket
term carries one inverse power of momentum, so int dp+ [bracket] = int dxi [...]
with the k+ cancelling INSIDE the bracket.

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

## The + prescription for the p+ integral

`plus_prescription.pdf`. Can the rapidity divergence be handled by adding and
subtracting the p+ = 0 term?

**Yes, but not for the reason "there is no transverse divergence".**

It works because, after the w' delta is used, the ONLY p+ dependence anywhere in the
integrand is the phase e^{-i(p+/k+)k.(y'-z)}: Wilson lines, sources and all four
kernels are p+-independent. So g(0) is just "the integrand with the phase set to 1",
and the subtraction is local in the transverse variables. Only the 1/p+ term needs
it; 1/(p+ + k+) is regular at p+ = 0 and 1/k+ has no pole.

    int_Lambda^P dp+ g(p+)/p+ = int_0^P dp+ [g(p+) - g(0)]/p+ + g(0) log(P/Lambda)

equivalently 1/xi -> [1/xi]_+ + delta(xi) log(Xi/lam). This reproduces the Ei form
exactly (checked for several kappa) -- it is an exact rearrangement, not an
approximation.

**The catch.** "No transverse divergence" was established WITH the phase. At p+ = 0
the phase is gone and g(0) IS transversely log divergent: the delta_jm delta_ik'
contraction leaves the BK real kernel (y'-z).(x-z)/[(y'-z)^2 (x-z)^2], which falls
only as 1/z^2. For p+ != 0 the phase cuts it at |z| ~ k+/(p+ |k|). That entanglement
is exactly why this row gives l_k^2 rather than l_k.

**What fixes it:** add the virtual partners so -2K -> M = (x-y')^2/[(x-z)^2(y'-z)^2],
which falls as 1/z^4. Then g(0) is finite, the + prescription isolates ONE rapidity
log times the BK kernel, and the remainder is the NLO impact factor.

**Scheme point:** the prescription isolates log((V-k+)/Lambda), not l_k =
log(k+/Lambda). The finite difference log((V-k+)/k+) is large in practice (4.595 for
V = 100 k+) and must match the convention of the evolution equation being matched to.

* `src/plusprescription.py` -> `RESULTS_plus.txt`

## Checking the + prescription split of the row

`plus_split_check.pdf`.

**Correct in the second (log) term:** the + sign (the minus in the bracket against the
overall minus), the prefactor 1/(8 pi^5) = (1/4pi^4)(1/2pi), the numerator
(y'-z).(x-z) (x'-y').(y-w) as the delta_jm delta_ik' contraction, and dropping the
xi-phase (g(0) has phase 1).

**Delete the leftover bracket** from the second term: the contraction is already in
the numerator, the 1/p+ is already the log, and there is no dp+ left for a free p+.

**The first term IS transversely divergent.** Writing the z-dependence of the three
index structures, with a = x'-y', b = y-w:

    S1 = (y'-z).a/(y'-z)^2 * (x-z).b/(x-z)^2      rho^2<S1> -> a.b/2
    S2 = (y'-z).(x-z)/[(y'-z)^2 (x-z)^2]          rho^2<S2> -> 1
    S3 = (y'-z).b/(y'-z)^2 * (x-z).a/(x-z)^2      rho^2<S3> -> a.b/2

all three fall only as 1/z^2. S1 and S3 keep the full phase e^{-i xi k.(y'-z)}, so
their z-integrals converge for every xi > 0 and the leftover log(1/xi) integrates
fine (no 1/xi pole on those terms). S2 does not: [1/p+]_+ replaces the phase by
(phase - 1), and the -1 leaves the BARE BK real kernel with nothing to cut it.
Measured rho^2 |<(phase-1) S2>| stays at ~1 out to rho = 1e4, while
rho^2 |<phase S2>| falls 1.4e-1 -> 4.4e-2 -> 6.9e-3. So the first term diverges
logarithmically at large |z| for every p+. The same divergence sits in the second
term (it IS g(0)); they cancel only in the sum.

**Second, separate problem:** the unsplit p+ integral has a finite V -> infinity
limit, but each split piece slides by log(V). Fix by subtracting only below k+,
[1/p+]_+^(k+) g = [g(p+) - theta(k+ - p+) g(0)]/p+, which isolates l_k = log(k+/Lambda)
and makes each piece V-stable.

**The real fix:** collect the virtual diagrams so -2K -> M = (x-y')^2/[(x-z)^2(y'-z)^2]
(angular average 1+1-2 = 0 at order rho^-2, so M falls as rho^-4). Only then is g(0)
transversely finite and the split legitimate. Self-check: with the virtuals in, the
answer must carry a SINGLE l_k; a surviving l_k^2 means the subtraction term is still
transversely divergent.

* `src/plus_split_check.py` -> `RESULTS_plussplit.txt`

## Does the + term's z-integral converge? (computed, not power-counted)

`plus_term_divergence.pdf`. The truncated integral F(R) = int_{|z|<R} d^2z f(z) for
each structure; a convergent integral has dF/dlogR -> 0.

    structure           dF/dlogR (xi=0.3)              dF/dlogR (xi=0.05)
    S1 x ph             -0.00920 +0.00014 -0.00891     +0.45878 -0.00084 +0.00011
    S2 x ph             -0.00235 +0.00006 -0.00022     +0.15187 -0.00081 +0.00005
    S3 x ph             -0.00889 +0.00014 -0.00891     +0.45879 -0.00085 +0.00011
    S2 x (ph-1)         -6.28585 -6.28313 -6.28352     -6.13162 -6.28400 -6.28314   <- -2pi
    S2 bare (= g(0))    +6.28350 +6.28319 +6.28319     +6.28350 +6.28319 +6.28319   <- +2pi

  (decades R = 1e2 -> 1e3 -> 1e4 -> 1e5;  -2pi = -6.28319)

S1 (on 1/(p+ + k+)) and S3 (on 1/k+) keep the full phase and converge; S2 with its
phase converges too. But [1/p+]_+ replaces the phase by (ph - 1), and the -1 leaves
the bare BK real kernel: slope -2pi, i.e. LOG DIVERGENT, for every xi. The second
term is g(0), slope +2pi. Equal and opposite: they cancel only in the sum.

**The split is exact but not usable.** Checked separately that keeping the lower
limit at Lambda in the + term makes the rearrangement exact to 1e-31 (pushing it to
0 costs O(Lambda)) -- so the expression is a correct rearrangement of a finite
quantity into two infinite pieces.

The only fix is the virtual diagrams: -2K -> M = (x-y')^2/[(x-z)^2 (y'-z)^2], which
falls as rho^-4. Then use the theta(k+ - p+) subtraction so each piece is also
V-stable and the isolated log is l_k.

* `src/zdiv_direct.py`, `src/zdiv_light2.py`, `src/exactness.py` -> `RESULTS_zdiv.txt`

## Transverse UV vs IR: which divergence blocks the + prescription

`transverse_uv_vs_ir.pdf`.

**No transverse UV divergence -- not at y' -> z, not anywhere.** At y' = z only ONE
of the four kernels is singular, so the integrand is O(1/rho) and int rho drho/rho
converges. It does better than that: the leading term integrates to zero by parity
(int_{|r|<eps} d^2r r^m/r^2 = 0), so the small-disc integral is O(eps^2). Measured
|int_{|y'-z|<eps} d^2y'| with ratios 0.0100 per decade for every structure -- with
the phase, without it, and subtracted:

    structure             eps=1e-01    eps=1e-02    eps=1e-03   ratio/decade
    T1 x ph               3.386e-04    3.386e-06    3.386e-08   0.0100
    T2 x ph               3.938e-03    3.938e-05    3.938e-07   0.0100
    T3 x ph               1.178e-02    1.178e-04    1.178e-06   0.0100
    T2 x (ph-1)           5.900e-04    5.870e-06    5.870e-08   0.0100
    T2 bare (g(0))        3.894e-03    3.894e-05    3.894e-07   0.0100

Multi-point collapses also converge: making two kernels singular at once needs two
points to collide, which costs more measure than it gains integrand (net rho^1, rho^2,
rho^3 drho/rho for 2, 3, 4 points colliding).

**But y' -> z is still the crucial region:** it is where ph -> 1, i.e. where the
integrand stops depending on p+, so it is the SOURCE of the rapidity divergence.
That is a longitudinal divergence in a transverse region, not a transverse one, and
the + prescription is built for it: ph - 1 = -i xi k.(y'-z) + O(r^2) supplies a
factor xi (killing the 1/xi) AND a factor rho (softening 1/rho to O(1)). The
subtracted disc integral is 6.7x smaller than the unsubtracted one at every eps.

**The obstruction is at the far end and it is INFRARED:** |y'-z| -> infinity, where
the phase is the only suppression and the subtraction removes it. The criterion is
"no transverse divergence", not "no transverse UV divergence".

Structural reason: in a dipole calculation large-z convergence comes from the KERNEL
(full BK M ~ rho^-4) and short distance from the COLOUR structure vanishing. Here
short distance is fine because each kernel is single, but only the real piece
-2K ~ rho^-2 is present. The colour does not help either: as |z| -> infinity,
U(z) -> 1 and both brackets become z-INDEPENDENT, not zero.

**On rewriting the log:** log[(V-k+)/Lambda] = log(V/Lambda) + log(1 - k+/V) is exact
and harmless, but it does not help -- the divergent object is the COEFFICIENT g(0),
not the logarithm. Rewriting distributes the same infinite coefficient over more
terms. It earns its keep only afterwards, once g(0) is finite.

* `src/uv_ypz.py` -> `RESULTS_uvypz.txt`

## What the exponential integrals mean (Lambda -> 0 vs y' -> z)

`exponential_integrals.pdf`. After the p+ integration the pole term gives
-[Ei(-i kap Xi) - Ei(-i kap lam)] with kap = k.(y'-z), lam = Lambda/k+.

**Why it blows up in both limits:** Ei(-i u) = gamma + log(i u) + O(u) and the
argument is the PRODUCT u = kap*lam. Four very different (kap, lam) with the same
product give an identical Ei (-13.23829489 - 1.570797j). Mathematica cannot tell the
two limits apart. The log splits as

    log(i kap lam) = log(lam) + log(i k.(y'-z)) = -l_k + log(i k.(y'-z))

**The y' -> z piece is NOT a transverse divergence.** A log is integrable in 2D:
log(y'-z)^2 is also "undefined" at y'=z, yet int d^2(y'-z) log(y'-z)^2 is finite.
Measured disc integrals of the p+-integrated integrand: 9.2229e-04, 1.0765e-05,
1.0785e-07 for eps = 1e-1, 1e-2, 1e-3, i.e. ratio 0.01 per decade, F ~ eps^2 (the
leading r^m/r^2 log r still averages to zero over angles). So the + prescription is
not obstructed there.

**The same log DOES bite at large |y'-z|.** Ei(-i kap Xi) -> -i pi sgn(kap), so the
bracket -> -l_k + gamma + log(i k.(y'-z)), growing like log|r| while the kernel falls
as 1/r^2 with nonzero angular average (rho^2 <T2> = +0.028221, constant). Hence
F ~ -c l_k log R + (c/2) log^2 R. The truncated integral confirms it: the slope does
not settle but GROWS by a constant 0.4083 per decade

    F(1e2..1e6) = -6.48980, -10.05680, -12.68358, -14.37029, -15.12141
    dF/dlogR    =           -1.5491,   -1.1408,   -0.7325,   -0.3262
    differences =                      +0.4083,   +0.4083,   +0.4063

exactly matching c log(10) with c = 2 pi x 0.028221 = 0.17732. The (c/2) log^2 R
piece is where l_k^2 comes from once R is cut at 1/(xi|k|).

**Mathematica caution:** Ei(-i u) = -E_1(i u) - i pi sgn(u) for real u (checked to
1e-10). The branch term cancels in the difference since Xi, lam > 0, but do NOT drop
Ei(-i kap Xi) when taking V -> infinity: it tends to -i pi sgn(kap), not zero.

* `src/ei_limits.py`, `src/ei_slope.py` -> `RESULTS_ei.txt`

## Assembling the row from the definitions of A and C

`row_assembly_check.pdf`. Does the quoted first equation follow from

    -A^a_i(k+,z) = -(1/(sqrt2 pi))(i g/sqrt(k+)) int_x (x-z)^i/(x-z)^2 rho^a(x)
    C^{bca}_{1jki}(k+,z;p+,x) = -(g f^{abc}/(2 pi sqrt(2k+))) sqrt(p+(k+-p+))/k+
                                 delta^(2)(z + (p+/k+) x) (x^m/x^2) [ ... ]

**Yes, exactly.** Symbolic check:

- three A's and one C give g^4 and (-i)(-1)(i)(i) = -i
- numeric: (2/(2pi)^3) x (1/(sqrt2 pi))^3 x 1/(2 sqrt2 pi) = (1/(2pi)^3)(1/(4 pi^4))
- k+ powers: C's sqrt(p+ k+) cancels the 1/sqrt(p+) and 1/sqrt(k+) of the two A's,
  leaving 1/(p+ + k+)^2
- **the delta Jacobian is what makes it 1/k+^2.** C's delta is
  delta^(2)(y'-w' + (p+/K)(w'-z)) = (K/k+)^2 delta^(2)[w' - ((p+ + k+)y' - p+ z)/k+],
  and 1/(p+ + k+)^2 x (p+ + k+)^2/k+^2 = 1/k+^2. Drop the Jacobian and you would get
  1/(p+ + k+)^2 instead; the quoted 1/k+^2 is right.
- full prefactor -i g^4/(32 pi^7 k+^2) vs quoted -i g^4/(32 pi^7 k+^2): difference 0
- index bracket: C_{1jki} -> C_{1jik'} is the relabelling j->j, k->i, i->k', giving
  [d_k'm d_ij - (K/p+) d_jm d_ik' - (K/(K-p+)) d_im d_jk'] with K = p+ + k+ and
  K - p+ = k+. Checked over all 16 index assignments: 0 mismatches.
- all four WW kernels and both colour brackets reproduce exactly

Not derivable from the definitions (they are inputs): the limits int_Lambda^(V-k+),
and the convention that A(K,-z) carries (x-z)^i/(x-z)^2.

* `src/assemble_row.py` -> `RESULTS_assemble.txt`

## The second row: same transverse structure, new colour

`second_row_check.pdf`. This row differs from the previous one ONLY in the ket
colour bracket:

    previous: U^{bc}(z) Abar^a_i(k+,-w) Abar^c_j(p+,-z) - Abar^a_i(k+,-w) A^b_j(p+,-z)
    this row: U^{bc}(z) Abar^c_j(p+,-z) Abar^a_i(k+,-w)
              - U^{ac}(w) Abar^d_j(p+,-z) U^{db}(z) A^c_i(k+,-w)

Each term still has one A(p+,-z) and one A(k+,-w) at the same momenta and positions,
and the bra bracket and C are untouched, so the whole assembly (prefactor
-(1/(2pi)^3)(i g^4 f/(4 pi^4))(1/k+^2), the delta Jacobian, the index bracket) and
the w' integration carry over verbatim. Both quoted equations are correct.

**Typo:** x' is missing from both integration measures although it appears in
(x'-y')^k'/(x'-y')^2 and rho^{e'}(x'). Same slip as in the previous row.

**Transverse divergences: identical.** Same four kernels, so no UV divergence, but
the delta_jm delta_ik' structure again leaves the BK real kernel
(y'-z).(x-z)/[(y'-z)^2(x-z)^2] ~ 1/z^2 in the z-integration, uncut once the phase is
removed at p+ = 0.

**Does the new colour save the large-|z| region? No.** At U(z) -> 1 (checked to
2.7e-15 with random adjoint Wilson lines):

    bra (unchanged)    -> d^{bd'}[U^{e'c'}(y') - U^{c'e'}(x')] rho^{e'}(x')   mag 1.9703
    ket, previous row  -> U^{ad}(y)rho^d(y)[U^{be}(x)rho^e(x) - rho^b(x)]     mag 7.3631
    ket, this row      -> U^{be}(x)rho^e(x)[U^{ad}(y) - U^{ad}(w)]rho^d(y)    mag 4.7409

all z-independent and nonzero, so the log coefficient survives in both rows.

**What the new colour DOES buy:** its limit is a difference of Wilson lines at y and
w (the previous row gave a difference in x), so it vanishes linearly as U(w) -> U(y)
-- verified with U(w) = U(y) exp(eps X). That softens exactly the coincidence where
(y-w)^i/(y-w)^2 is singular: a real-virtual structure in the (y,w) pair, the analogue
in the w integration of what the virtual diagrams must still supply for z.

**So the + prescription:** yes, writable in the same form, exact with the lower limit
kept at Lambda -- and with the same caveat that neither term is separately finite.

* `src/colour_limits.py` -> `RESULTS_colour.txt`

## The U(z) -> 1 version of the second row

`src/uz_limit.py` -> `RESULTS_uzlimit.txt`. A later draft replaced the colour
brackets by

    bra: [ U^{e'c'}(y') rho^{e'}(x') - U^{c'e'}(x') rho^{e'}(x') ]
    ket: [ U^{ce}(x) U^{ad}(y) rho^e rho^d - U^{ad}(w) U^{ce}(x) rho^e rho^d ]

with f^{c'd'a} -> f^{c'ca}.

**That is exactly the old structure with U(z) -> 1**, executed correctly: the full
contraction f^{c'd'a} Bra^{c'd'b} Ket^{ab} at U(z) = 1 equals the new
f^{c'ca} Bra^{c'} Ket^{ca} to 4.4e-16, while at a generic U(z) the two differ
(1.0476 vs 5.1460). The delta^{bd'} left by U(z) -> 1 sets d' = b, which is then
relabelled b -> c.

**It does not follow from unitarity.** Of the four U(z) products only the CROSS ones
collapse:

    bra1 x ket1 : U^{d'b}(z) U^{bc}(z)  -> U(z)^2, NOT delta
    bra1 x ket2 : U^{d'b}(z) U^{db}(z)  -> delta     (unitarity)
    bra2 x ket1 : U^{bd'}(z) U^{bc}(z)  -> delta     (unitarity)
    bra2 x ket2 : U^{bd'}(z) U^{db}(z)  -> U(z)^2, NOT delta

So U(z) -> 1 is a substitution, not an identity, and z is an integration variable:
it discards the whole z-dependence of the colour, which is exactly what the virtual
diagrams must cancel against. It is useful for ONE thing -- exhibiting the
coefficient of the large-|z| logarithm, since the colour factor then comes out of
the z-integral and leaves int d^2z (y'-z).(x-z)/[(y'-z)^2 (x-z)^2], manifestly log
divergent.

## Does the bare BK real kernel's z-integral have a UV divergence?

`src/single_kernel_z.py` -> `RESULTS_singlek.txt`.

    int_{|z|<R} d^2z  (y'-z).(x-z)/[(y'-z)^2 (x-z)^2]  =  pi log( R^2 / (y'-x)^2 )

exactly, with no extra constant (checked against quadrature at R = 20, 100, 500,
2000: 17.271015/17.268607, 27.381075/27.381004, 37.493378/37.493400,
46.203719/46.203745). It follows from K = Re 1/[(zbar - ybar')(z - x)] and
L(a,b) = int_{|z|<R} d^2z/[(zbar-abar)(z-b)] = pi log(R^2/|a-b|^2).

**No UV divergence, and the short-distance region contributes EXACTLY ZERO.**
The angular average of K on a circle of radius rho about z = y' is

    <K> = (1/rho) Re \oint (d zeta / 2 pi i) 1/(rho zeta - c),   c = x - y'

whose pole zeta = c/rho lies outside the unit circle for rho < |c|, so the average
vanishes identically there. Measured: <K> = 2e-16, 6e-17, 3e-17, -1e-16 at
rho = 0.1, 0.5, 1.0, 1.27 with |c| = 1.2806; and rho^2<K> = 1.000000 for
rho = 1.29, 2, 10. The whole integral comes from rho > |x-y'|.

The only short-distance scale in the answer is |y'-x|, an external separation --
no regulator appears. The R-dependence is 2 pi log R: the divergence is INFRARED,
coefficient 2 pi, matching the slope +6.28319 measured for this structure.

The log(1/(y'-x)^2) is the y' -> x singularity: at exactly y' = x the kernel
degenerates to 1/(x-z)^2 and the z-integral IS log divergent at z -> x, but that is
a measure-zero configuration and a log is integrable against d^2(y'-x).

## Does this row contribute to DGLAP?

`src/dglap_check.py` -> `RESULTS_dglap.txt`.  **No** -- three independent reasons.

**1. No transverse UV divergence.** DGLAP logs are collinear logs from the transverse
short-distance region. This row has none: each of the four kernels is a single
1/|sep|, no pair of points appears in two kernels, so there is no 1/rho^2 anywhere;
and for the bare BK real kernel the short-distance region contributes IDENTICALLY
ZERO (the angular average vanishes for rho < |x-y'|).

**2. No collinear endpoint in the longitudinal integral.** Poles of the three bracket
terms in xi = p+/k+:

    1/(1+xi)   pole at xi = -1   -> none in [0,1]
    -1/xi      pole at xi =  0   -> the SOFT endpoint
    -1         no pole

Only xi -> 0. There is no xi -> 1 pole, i.e. no collinear endpoint. The reason is the
momentum routing: C's third term carries k+/(k+ - p+), a genuine 1/(1-xi) collinear
pole, but C is used with first argument K+ = p+ + k+, so it becomes (p+ + k+)/k+ --
finite for every p+ > 0. The collinear pole is switched off by routing.

**3. The coefficient is not a splitting function.** C_UV(xi) = P_gg/(2N_c) =
xi(1-xi) + xi/(1-xi) + (1-xi)/xi has poles at BOTH 0 and 1; this row's divergent
structure is the bare eikonal -1/xi, with the xi(1-xi) and xi/(1-xi) pieces absent.

**Contrast: the earlier UV-projected row DOES contribute to DGLAP.** There the
light-cone vertex supplied a second 1/(x-z)^2, giving a true 1/r^2 and a true
transverse UV log whose coefficient was exactly C_UV(xi) = P_gg/(2N_c). Its moment
int C_UV dxi = 2 l_k - 11/6 (checked: 39.613198424 vs 39.613198341 at lam = 1e-9),
and (2 l_k - 11/6) x 2N_c = 4 N_c l_k - b at N_f = 0 -- the constant is the
beta-function coefficient, i.e. coupling renormalisation. That is the DGLAP /
running-coupling structure; this row has none of it.

What this row contributes to instead is BFKL/BK rapidity evolution: the 1/xi pole
times the BK real kernel is the JIMWLK real-emission term, and the transverse log it
carries is the large-distance (infrared) one, not a collinear one.
