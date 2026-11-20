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


def vint(expr):
    """int_0^1 dv [v(1-v)]^{-eps} * expr, expr a Laurent polynomial in v,(1-v) -> Beta functions."""
    expr = sp.expand(expr)
    res = 0
    for term in sp.Add.make_args(expr):
        c, (a, b) = coeff_powers(term)
        res += c * sp.gamma(1 + a - eps) * sp.gamma(1 + b - eps) / sp.gamma(2 + a + b - 2 * eps)
    return res


w = sp.Symbol('w', positive=True)  # stands for (1-v)


def coeff_powers(term):
    # term = c * v**a * w**b
    c = term
    a = 0
    b = 0
    for f in sp.Mul.make_args(term):
        if f.has(v):
            base, ex = f.as_base_exp()
            a += ex
        elif f.has(w):
            base, ex = f.as_base_exp()
            b += ex
    c = sp.simplify(term / (v ** a * w ** b))
    return c, (a, b)


def series0(expr, n=1):
    return sp.series(expr, eps, 0, n).removeO()


def finite_part(real_density, ct_density):
    tot = sp.simplify(real_density + ct_density)
    s = sp.series(tot, eps, 0, 1).removeO()
    pole = sp.simplify(s.coeff(eps, -1))
    fin = sp.simplify(s.coeff(eps, 0))
    return pole, fin


# ---------------------------------------------------------------- (A) e+e-
# y13 = z v, y23 = z w (w = 1-v), y12 = 1 - z ; measure z^{1-2eps} (1-z)^{-eps} [v w]^{-eps} dv
y12, y13, y23 = 1 - z, z * v, z * w
Mhat = (1 - eps) * (y13 / y23 + y23 / y13) + 2 * y12 / (y13 * y23) - 2 * eps
Iv = vint(sp.expand(Mhat))
pref = (4 * sp.pi) ** eps / sp.gamma(1 - eps)          # mu = Q = 1 units
real_A = pref * z ** (1 - 2 * eps) * (1 - z) ** (-eps) * Iv   # (1/sigma0) dsigma_g/dz / (alpha_s C_F/2pi)
Pgq = (1 + (1 - z) ** 2) / z                            # P_gq / C_F
ct_ND = (4 * sp.pi) ** eps / sp.gamma(1 - eps) / eps * 2 * Pgq
ct_J2 = ct_ND * z ** (-2 * eps)
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
IvB = vint(sp.expand(MhatB))
# n(z)/Gamma0 = (alpha_s C_A/2pi) * pref * z^{1-2eps}(1-z)^{-eps} IvB ; per Born gluon divide by 2
real_B = sp.Rational(1, 2) * pref * z ** (1 - 2 * eps) * (1 - z) ** (-eps) * IvB   # units (alpha_s/2pi) C_A
PggR = 2 * (z / (1 - z) + (1 - z) / z + z * (1 - z))   # P_gg+ / C_A
ctB_ND = (4 * sp.pi) ** eps / sp.gamma(1 - eps) / eps * PggR
ctB_J2 = ctB_ND * z ** (-2 * eps)
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
