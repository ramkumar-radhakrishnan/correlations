"""Referee 3, collinear cross-check: one-hadron inclusive gluon coefficient functions
in d = 4 - 2 eps, with the MS-bar counterterm in two Jacobian conventions.

(A) e+e- -> gamma* -> q qbar g, observed gluon with energy fraction z (T+L summed).
    Literature (Kretzer hep-ph/0003177, App., from AEMP/FP/NW), in units (alpha_s/2pi) C_F:
       c_T^g + c_L^g = 2 (1+(1-z)^2)/z * ( ln(1-z) + 2 ln z + ln(Q^2/mu^2) ).
(B) H -> g g (scalar phi G G), NLO real H -> g g g, observed gluon, z < 1, N_f = 0.
    Literature (Almasy-Moch-Vogt 1107.2263, Eq. (A.1)), a_s = alpha_s/(4 pi), mu = M:
       c_{phi,g}^{(1)}(x)|_{x<1} = C_A [ 11/3 (1+x+x^2-1/(1-x)) + p_gg(x) (8 ln x + 4 ln(1-x)) ].

Conventions compared (CS = Catani-Seymour hep-ph/9605323 Eq. (saCinc)):
   number density ("ND"): partonic counterterm with dz/z^{2-2eps} for the fully differential
       cross section, i.e. dz/z for d sigma/dz  (this is the convention of the draft, xi^{-d});
   "2D Jacobian" ("J2"): dz/z^2 also in d dimensions (the CXY form), i.e. dz/z^{1+2eps}
       for d sigma/dz.
The real emission is normalised to the d-dimensional Born; prefactor (4 pi mu^2/Q^2)^eps/Gamma(1-eps),
counterterm (alpha_s/2pi) (4 pi mu^2/mu_F^2)^eps/Gamma(1-eps) (1/eps) P(z).  We set mu = mu_F = Q.
Run: python3 -I c2_sia_coefficients.py
"""
import sympy as sp

eps, z, v = sp.symbols('epsilon z v', positive=True)


w = sp.Symbol('w', positive=True)  # stands for (1-v)
Lz, L1z, L4pi = sp.symbols('Lz L1z L4pi')  # ln z, ln(1-z), ln(4 pi)
_beta_cache = {}


def beta_series(a, b, order=2):
    """Laurent series in eps of int_0^1 dv v^(a-eps) (1-v)^(b-eps) = B(1+a-eps, 1+b-eps), integers a,b."""
    key = (a, b)
    if key not in _beta_cache:
        B = sp.gamma(1 + a - eps) * sp.gamma(1 + b - eps) / sp.gamma(2 + a + b - 2 * eps)
        _beta_cache[key] = sp.series(B, eps, 0, order).removeO()
    return _beta_cache[key]


def vint(expr):
    """int_0^1 dv [v w]^{-eps} * expr, expr a Laurent polynomial in v, w=(1-v); returned as eps series."""
    expr = sp.expand(expr)
    res = 0
    for term in sp.Add.make_args(expr):
        a = b = 0
        for f in sp.Mul.make_args(term):
            base, ex = f.as_base_exp()
            if base == v:
                a += ex
            elif base == w:
                b += ex
        c = sp.simplify(term / (v ** a * w ** b))
        res += c * beta_series(int(a), int(b))
    return sp.expand(res)


def pref_series(extra_log):
    """(4 pi)^eps / Gamma(1-eps) * exp(eps*extra_log) as a series to O(eps)."""
    t = sp.symbols('t')
    ser = sp.series(sp.exp(t * (L4pi + extra_log)) / sp.gamma(1 - t), t, 0, 3).removeO()
    return sp.expand(ser.subs(t, eps))


def laurent(expr):
    expr = sp.expand(expr)
    pole = sp.simplify(expr.coeff(eps, -1))
    fin = sp.simplify(expr.coeff(eps, 0))
    dpole = sp.simplify(expr.coeff(eps, -2))
    return dpole, pole, fin


def finite_part(real_density, ct_density):
    dpole, pole, fin = laurent(real_density + ct_density)
    assert dpole == 0
    back = {Lz: sp.log(z), L1z: sp.log(1 - z)}
    return sp.simplify(pole.subs(back)), sp.simplify(fin.subs(back))


# ---------------------------------------------------------------- (A) e+e-
# y13 = z v, y23 = z w (w = 1-v), y12 = 1 - z ; measure z^{1-2eps} (1-z)^{-eps} [v w]^{-eps} dv
y12, y13, y23 = 1 - z, z * v, z * w
Mhat = (1 - eps) * (y13 / y23 + y23 / y13) + 2 * y12 / (y13 * y23) - 2 * eps
Iv = vint(Mhat)
# (1/sigma0) dsigma_g/dz / (alpha_s C_F/2pi) = (4pi)^eps/Gamma(1-eps) z^{1-2eps} (1-z)^{-eps} Iv   (mu = Q = 1)
real_A = sp.expand(z * pref_series(-2 * Lz - L1z) * Iv)
Pgq = (1 + (1 - z) ** 2) / z                            # P_gq / C_F
ct_ND = sp.expand(pref_series(0) / eps * 2 * Pgq)
ct_J2 = sp.expand(pref_series(-2 * Lz) / eps * 2 * Pgq)
lit_A = 2 * Pgq * (sp.log(1 - z) + 2 * sp.log(z))

print("(A) e+e- -> gX, T+L, units (alpha_s/2pi) C_F, mu = Q")
for name, ct in (("number density (draft, CS)", ct_ND), ("2D Jacobian (CXY form)", ct_J2)):
    pole, fin = finite_part(real_A, ct)
    d = sp.simplify(sp.expand_log(fin - lit_A, force=True))
    print("  %-28s pole: %s ; finite - literature = %s" % (name, sp.simplify(pole), sp.factor(d)))
# numerical spot check of the finite part
pole, finND = finite_part(real_A, ct_ND)
for zz in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(4, 5)):
    print("    z=%s  finite(ND)=%.12f  literature=%.12f" % (zz, float(finND.subs(z, zz)), float(lit_A.subs(z, zz))))

# energy sum rule test of the convention: shift of (1/2) int z [2 C_q + C_g] dz caused by J2
# delta C_q = -P_qq ln z^2 (per quark), delta C_g = -2 P_gq ln z^2 ; P in units of C_F
x = sp.Symbol('x', positive=True)
Pqq_reg = (1 + x ** 2) / (1 - x)   # ln x^2 kills the endpoint, no plus/delta needed
Pgq_x = (1 + (1 - x) ** 2) / x
shift = sp.integrate(sp.expand(-x * sp.log(x ** 2) * (Pqq_reg + Pgq_x)).rewrite(sp.log), (x, 0, 1))
print("  J2 shift of the energy sum (1/2) int z [2C_q + C_g], units (alpha_s/2pi) C_F:", sp.nsimplify(sp.simplify(shift)), "=", float(shift))
# check that the literature coefficients obey the energy sum: (1/2) int z [2 c_q + c_g] = 3/2
cqT = None  # c_q requires distributions; done in c1 (not needed here)

# ---------------------------------------------------------------- (B) H -> g g g
CA = sp.Symbol('C_A', positive=True)
MhatB = (1 + y12 ** 4 + y13 ** 4 + y23 ** 4) / (y12 * y13 * y23)
IvB = vint(MhatB)
# n(z)/Gamma0 = (alpha_s C_A/2pi) * pref * z^{1-2eps}(1-z)^{-eps} IvB ; per Born gluon divide by 2
real_B = sp.expand(sp.Rational(1, 2) * z * pref_series(-2 * Lz - L1z) * IvB)   # units (alpha_s/2pi) C_A
PggR = 2 * (z / (1 - z) + (1 - z) / z + z * (1 - z))   # P_gg+ / C_A
ctB_ND = sp.expand(pref_series(0) / eps * PggR)
ctB_J2 = sp.expand(pref_series(-2 * Lz) / eps * PggR)
pgg = 1 / z - 2 + z - z ** 2 + 1 / (1 - z)
lit_B = sp.Rational(1, 2) * (sp.Rational(11, 3) * (1 + z + z ** 2 - 1 / (1 - z)) + pgg * (8 * sp.log(z) + 4 * sp.log(1 - z)))  # a_s -> alpha_s/2pi: /2
print("(B) H -> gX, z<1, N_f=0, units (alpha_s/2pi) C_A, mu = M_H")
for name, ct in (("number density (draft, CS)", ctB_ND), ("2D Jacobian (CXY form)", ctB_J2)):
    pole, fin = finite_part(real_B, ct)
    d = sp.simplify(sp.expand_log(fin - lit_B, force=True))
    print("  %-28s pole: %s ; finite - literature = %s" % (name, sp.simplify(pole), sp.factor(d)))
pole, finB = finite_part(real_B, ctB_ND)
for zz in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(4, 5)):
    print("    z=%s  finite(ND)=%.12f  literature=%.12f" % (zz, float(finB.subs(z, zz)), float(lit_B.subs(z, zz))))
