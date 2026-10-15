"""One-rho (commutator) parts of the vacuum one-gluon terms at order g^3, Calculation_main.tex Sec. 2.2.

(a) Psi^3 ([rho,rho] part of Eq. 3grhorho)  -> Eq. onerho2
(b) Psi^4 + Psi^5 ([rho,rho] part of psif + psiff) -> Eq. 4+5onerho

The transverse integral over the loop momentum is done in dimensional regularisation, d = 2 - eps, MS-bar:
    int d^2p f  ->  (2 pi)^2 mu^eps int d^dp/(2 pi)^d f,   ln(4 pi) - gamma_E absorbed in mu_MSbar.
Feynman parameter: 1/(A B) = int_0^1 dx (x A + (1-x) B)^-2.  Odd powers of l vanish, (k.l)^2 -> k^2 l^2/d.
    int (2pi)^2 d^dl/(2pi)^d          1/(l^2+M^2)^2 = pi/M^2
    int (2pi)^2 d^dl/(2pi)^d        l^2/(l^2+M^2)^2 = pi [2/eps - ln(M^2/mu^2) - 1]
    int (2pi)^2 d^dl/(2pi)^d    (k.l)^2/(l^2+M^2)^2 = (k^2/2) pi [2/eps - ln(M^2/mu^2)]
The result is the coefficient of  rho(-k) k^i/k^2 |g_i(k)>  as a function of xi, in units of g^3 N_c/(32 pi^{7/2}),
before the xi integration; it is compared with the xi-integrands of the notes.  The exact xi integrals are done after."""
import sympy as sp
import mpmath as mp

xi, c, k2, eps, L2 = sp.symbols('xi c k2 epsilon Lmu', positive=True)   # L2 = ln(k^2/mu^2)
s1, s2 = sp.symbols('s1 s2')            # s1 = k.l, s2 = l^2
t = sp.symbols('t')                     # bookkeeping for powers of l
d = 2 - eps


def lint(num, M2):
    """(2pi)^2 int d^dl/(2pi)^d num/(l^2+M^2)^2, num a polynomial in s1, s2 (the eps -> 0 limit, poles kept)"""
    num = sp.expand(num)
    P = sp.Poly(num, s1, s2)
    res = 0
    for (a1, a2), cf in P.terms():
        if a1 % 2 == 1:
            continue
        if (a1, a2) == (0, 0):
            res += cf * sp.pi / M2
        elif (a1, a2) == (0, 1):
            res += cf * sp.pi * (2 / eps - sp.log(M2 / k2) - L2 - 1)
        elif (a1, a2) == (2, 0):
            res += cf * k2 / 2 * sp.pi * (2 / eps - sp.log(M2 / k2) - L2)
        else:
            raise ValueError((a1, a2))
    return res


def kdot_vectors(c):
    """with p = c k + l: return the scalars p.k, p^2, (k-p)^2, (k-p).k in terms of s1, s2"""
    pk = c * k2 + s1
    pp = c**2 * k2 + 2 * c * s1 + s2
    qq = (1 - c)**2 * k2 - 2 * (1 - c) * s1 + s2
    qk = k2 - pk
    return pk, pp, qq, qk


# ---------------------------------------------------------------- (a) Psi^3
def psi3_onerho_xi():
    """Eq. 3grhorho, [rho,rho] part: coefficient -g^3 N_c rho(-k) B_l/(32 pi^{9/2} xi (k-p)^2 D k^2), D = (1-xi)p^2 + xi (k-p)^2.
    Feynman: x(k-p)^2 + (1-x) D = (p - c k)^2 + c(1-c) k^2, c = x + (1-x) xi, dx = dc/(1-xi), c in (xi, 1).
    Returns the xi-integrand F_a(xi) with   onerho2 = int dxi F_a(xi) * g^3 N_c rho k^i/(32 pi^{7/2} k^2)."""
    pk, pp, qq, qk = kdot_vectors(c)
    pq = pk - pp                               # p.(k-p)
    # k.B, B_l = (2p - 2 xi k)_l p.(k-p) + (2 k.p - 2/xi p^2)(k-p)_l + ((1+xi)/(1-xi)(k-p)^2 - k^2 + p^2) p_l
    kB = (2 * pk - 2 * xi * k2) * pq + (2 * pk - 2 / xi * pp) * qk + ((1 + xi) / (1 - xi) * qq - k2 + pp) * pk
    M2 = c * (1 - c) * k2
    I = lint(kB, M2)                           # int d^2p k.B / (x A + (1-x) B)^2
    # coefficient of k^l/k^2:  k.J = k.I / ...;  J_l = k_l (k.J)/k^2  -> coefficient of k^l/k^2 is k.J
    integrand_c = sp.simplify(I / (1 - xi))    # dx = dc/(1-xi)
    Ic = sp.integrate(sp.expand(integrand_c), (c, xi, 1))
    # full: -g^3 N_c/(32 pi^{9/2}) (1/xi) (k.J)/k^2 k^l/k^2 ; in units of g^3 N_c k^l/(32 pi^{7/2} k^2): -(1/pi)(1/xi) Ic/k^2
    Fa = sp.simplify(-Ic / (sp.pi * xi * k2))
    return Fa


# ---------------------------------------------------------------- (b) Psi^4 + Psi^5
def psi45_onerho_x(xv):
    """[rho,rho] part of psif + psiff at fixed xi = xv (a number).  In psif/psiff the final gluon is (xi k+, p) and the
    loop momentum is k.  rho^c(k-p) rho^a(-k) -> (i/2) f^{cae} rho^e(-p), f^{abc} f^{cae} = +N_c delta^{be}:
      psif : -i g^3 f^{abc} (...) -> -i (i/2) N_c = +N_c/2, but int d^2k (k-p)^i/D = 0 (odd in k-p);
      psiff:  i g^3 f^{abc} (...) ->  i (i/2) N_c = -N_c/2.
    psiff: xi^{3/2}/(16 pi^{9/2} k^2 p^2 D (1-xi)) V^i,  D = xi (k-p)^2 + (1-xi) p^2.
    Feynman: x k^2 + (1-x) D = a (k - s p)^2 + M2,  a = x + (1-x) xi,  s = (1-x) xi/a,  M2 = (1-x) p^2 - (1-x)^2 xi^2 p^2/a.
    Returns (pole, Lmu, finite) of F_b(xi) sqrt(xi)(1-xi), with
      4+5onerho = int dxi F_b g^3 N_c rho(-p) p^i/(32 pi^{7/2} p^2) |g_i(xi k+, p)>,  p^2 = 1 (k2 = 1)."""
    x = sp.symbols('x', positive=True)
    X = sp.nsimplify(xv)
    a = x + (1 - x) * X
    sh = (1 - x) * X / a
    pp_ = 1
    M2 = (1 - x) * pp_ - (1 - x)**2 * X**2 * pp_ / a
    kp = sh * pp_ + s1
    kk = sh**2 * pp_ + 2 * sh * s1 + s2
    qq = kk - 2 * kp + pp_
    pV = ((2 * X - 1) * kk - pp_ - (1 + X) / (1 - X) * qq) * kp + ((X - 1) / X * kp + (1 - X**2) / X * kk) * 2 * pp_
    # with k = s p + l the denominator is (a l^2 + M2)^2 = a^2 (l^2 + M2/a)^2; s1 = p.l, s2 = l^2, (p.l)^2 -> p^2 l^2/d
    I = lint(sp.expand(pV).subs(k2, 1), M2 / a) / a**2
    I = sp.expand(I.subs(k2, 1).subs(eps, 1 / t))        # t = 1/eps
    pole = sp.integrate(sp.simplify(sp.diff(I, t)), (x, 0, 1))
    lmu = sp.integrate(sp.simplify(sp.diff(I, L2)), (x, 0, 1))
    fin_expr = I.subs({t: 0, L2: 0})
    f = sp.lambdify(x, fin_expr, 'mpmath')
    fin = mp.quad(f, [0, 0.5, 1])
    pre = -sp.Rational(1, 2) * X**sp.Rational(3, 2) / (1 - X) * 2 / sp.pi   # -N_c/2, 1/16 = 2/32, units /pi
    norm = sp.sqrt(X) * (1 - X)
    return (float(pre * pole * norm), float(pre * lmu * norm), float(pre * norm) * float(fin))


if __name__ == '__main__':
    import numpy as np, os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'three_rho'))
    import numerics as nm
    _, F = nm.structure_constants(3)
    print('f^{abc} f^{cae} = -N_c delta^{be}:', np.allclose(np.einsum('abc,cae->be', F, F), -3 * np.eye(8)))
    print('f^{bad} f^{bae} = +N_c delta^{de}:', np.allclose(np.einsum('bad,bae->de', F, F), 3 * np.eye(8)))
    Fa = psi3_onerho_xi()
    print('(a) Psi^3 one-rho xi-integrand F_a(xi) =')
    sp.pprint(sp.collect(sp.expand(Fa), [1 / eps, L2]))
    print('(b) Psi^4+5 one-rho: F_b sqrt(xi)(1-xi) = pole/eps + Lmu*ln(p^2/mu^2) + finite, vs 4+5onerho:')
    for xv in (0.2, 0.5, 0.7):
        pol, lm, fi = psi45_onerho_x(xv)
        nt = (xv - 2) * np.log(xv) + (xv + 1) * np.log(1 - xv)
        print('  xi=%.2f  pole %.12f (notes %.12f)  Lmu %.12f (notes %.12f)  finite %.12f (notes %.12f)'
              % (xv, pol, 2 * (xv - 2), lm, -(xv - 2), fi, nt))
