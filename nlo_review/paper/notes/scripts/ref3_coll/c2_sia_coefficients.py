"""Referee 3, collinear cross-check: one-hadron inclusive gluon coefficient functions in d=4-2eps (SIA),
MS-bar counterterm with the number-density (d-dim) Jacobian vs the two-dimensional (CXY) Jacobian.
Literature: e+e- (Kretzer hep-ph/0003177 App., from AEMP/FP/NW): c_T^g+c_L^g = 2(1+(1-z)^2)/z [ln(1-z)+2ln z+ln Q^2/mu^2];
H->gX (Almasy-Moch-Vogt 1107.2263, Eq. (A.1)), a_s=alpha_s/4pi.  Standard MS-bar counterterm: Catani-Seymour
hep-ph/9605323 Eq. (6.17) "saCinc", dz/z^{2-2eps} for the fully differential cross section.
Run: python3 -I c2_sia_coefficients.py

v-integrals, int_0^1 dv [v(1-v)]^{-eps} v^a (1-v)^b = B(1+a-eps, 1+b-eps):
 (A) e+e-:  Mhat = (1-eps)(v/w + w/v) + 2(1-z)/(z^2 v w) - 2 eps
     Iv = 2(1-eps) B(2-eps,-eps) + 2(1-z)/z^2 B(-eps,-eps) - 2 eps B(1-eps,1-eps)
 (B) H->ggg: Mhat = [1+(1-z)^4 + z^4 (v^4+w^4)] / ((1-z) z^2 v w)
     Iv = [1+(1-z)^4]/((1-z) z^2) B(-eps,-eps) + 2 z^2/(1-z) B(4-eps,-eps)
 real density = z (4pi)^eps/Gamma(1-eps) z^{-2eps} (1-z)^{-eps} Iv    [x 1/2 per Born gluon in (B)]
 counterterm  = (4pi)^eps/Gamma(1-eps) (1/eps) P(z) x {1 (number density) | z^{-2eps} (2D Jacobian)}
"""
import sympy as sp

e, z = sp.symbols('epsilon z', positive=True)
Lz, L1z, L4 = sp.symbols('Lz L1z L4pi')


def ser(expr, n=2):
    return sp.expand(sp.series(expr, e, 0, n).removeO())


def B(a, b):
    return ser(sp.gamma(a) * sp.gamma(b) / sp.gamma(a + b))


pref = ser(sp.exp(e * L4) / sp.gamma(1 - e), 3)
zfac = ser(sp.exp(e * (-2 * Lz - L1z)), 3)
zJ2 = ser(sp.exp(e * (-2 * Lz)), 3)


def laurent_coeffs(expr):
    """coefficients of eps^-2, eps^-1, eps^0 of expr (a Laurent polynomial in eps with z-rational
    coefficients); robust against sympy merging eps into denominators."""
    f2 = sp.cancel(sp.together(expr * e ** 2))
    return [sp.simplify(sp.diff(f2, e, k).subs(e, 0) / sp.factorial(k)) for k in (0, 1, 2)]


def trunc(expr):
    c = laurent_coeffs(expr)
    return c[0] / e ** 2 + c[1] / e + c[2]


def finite(real, ct):
    dp, p, f = laurent_coeffs(real + ct)
    back = {Lz: sp.log(z), L1z: sp.log(1 - z)}
    return dp, sp.simplify(p.subs(back)), sp.simplify(f.subs(back))


# (A)
IvA = 2 * (1 - e) * B(2 - e, -e) + 2 * (1 - z) / z ** 2 * B(-e, -e) - 2 * e * B(1 - e, 1 - e)
realA = trunc(z * pref * zfac * IvA)
PgqA = 2 * (1 + (1 - z) ** 2) / z          # q and qbar, units C_F
litA = PgqA * (sp.log(1 - z) + 2 * sp.log(z))
print("(A) e+e- -> gX (T+L), units (alpha_s/2pi) C_F, mu = Q;  literature: 2(1+(1-z)^2)/z [ln(1-z) + 2 ln z]")
for name, fac in (("number density (draft, CS)", 1), ("2D Jacobian (CXY form)", zJ2)):
    ct = trunc(pref * fac / e * PgqA)
    dp, p, f = finite(realA, ct)
    print("   %-27s double pole %s, pole %s, finite - literature = %s" % (name, dp, p, sp.factor(sp.expand(sp.expand_log(f - litA, force=True)))))

# (B)
IvB = (1 + (1 - z) ** 4) / ((1 - z) * z ** 2) * B(-e, -e) + 2 * z ** 2 / (1 - z) * B(4 - e, -e)
realB = trunc(sp.Rational(1, 2) * z * pref * zfac * IvB)
PggB = 2 * (z / (1 - z) + (1 - z) / z + z * (1 - z))   # P_gg+ / C_A
pgg = 1 / z - 2 + z - z ** 2 + 1 / (1 - z)
litB = sp.Rational(1, 2) * (sp.Rational(11, 3) * (1 + z + z ** 2 - 1 / (1 - z)) + pgg * (8 * sp.log(z) + 4 * sp.log(1 - z)))
print("(B) H -> gX, z<1, N_f=0, units (alpha_s/2pi) C_A, mu = M_H;  literature: AMV (A.1)/2")
for name, fac in (("number density (draft, CS)", 1), ("2D Jacobian (CXY form)", zJ2)):
    ct = trunc(pref * fac / e * PggB)
    dp, p, f = finite(realB, ct)
    print("   %-27s double pole %s, pole %s, finite - literature = %s" % (name, dp, p, sp.factor(sp.expand(sp.expand_log(f - litB, force=True)))))
    if fac == 1:
        print("      finite (number density) =", sp.simplify(f))

# energy sum rule of the e+e- coefficients with the 2D Jacobian: shift of (1/2) int z [2 C_q + C_g]
x = sp.Symbol('x', positive=True)
sh = 0
for t in sp.Add.make_args(sp.apart(sp.cancel(-x * ((1 + x ** 2) / (1 - x) + (1 + (1 - x) ** 2) / x)), x)):
    sh += sp.integrate(t * 2 * sp.log(x), (x, 0, 1))
print("(A') 2D-Jacobian shift of (1/2) int z [2C_q + C_g] dz, units (alpha_s/2pi) C_F:", sp.simplify(sh), "=", float(sh))
print("     (standard coefficients give 3/2, i.e. sigma_tot = sigma_0 (1 + alpha_s/pi) for C_F = 4/3)")
# same for H -> gg (pure YM): shift = -int z P_gg ln z^2, units (alpha_s/2pi) C_A
sh2 = 0
for t in sp.Add.make_args(sp.apart(sp.cancel(-x * 2 * (x / (1 - x) + (1 - x) / x + x * (1 - x))), x)):
    sh2 += sp.integrate(t * 2 * sp.log(x), (x, 0, 1))
print("(B') 2D-Jacobian shift of int z C_g dz for H -> gX, units (alpha_s/2pi) C_A:", sp.simplify(sh2), "=", float(sh2))

# (B'') energy sum of the full AMV coefficient (N_f=0), a_s = alpha_s/(4 pi) units, per C_A:
#   c = 11/3 (1+x+x^2 - [1/(1-x)]_+) + p_gg (8 ln x - 4 H_1) + (67/9 + 8 zeta2) delta(1-x),  H_1 = -ln(1-x)
#   p_gg H_1 contains -[ln(1-x)/(1-x)]_+ ; int x [f]_+ = int (x-1) f
def Ireg(rat, lg):
    tot = 0
    for t in sp.Add.make_args(sp.apart(sp.cancel(rat), x)):
        tot += sp.integrate(t * lg, (x, 0, 1))
    return tot
pg_reg = 1 / x - 2 + x - x ** 2
mom = (sp.Rational(11, 3) * (sp.integrate(x * (1 + x + x ** 2), (x, 0, 1)) - sp.integrate(sp.cancel((x - 1) / (1 - x)), (x, 0, 1)))
       + 8 * (Ireg(x * pg_reg, sp.log(x)) + Ireg(x / (1 - x), sp.log(x)))
       + 4 * (Ireg(x * pg_reg, sp.log(1 - x)) + sp.integrate(sp.cancel((x - 1) / (1 - x)) * sp.log(1 - x), (x, 0, 1)))
       + sp.Rational(67, 9) + 8 * sp.zeta(2))
import mpmath as mp
num = (mp.mpf(11) / 3 * (mp.mpf(13) / 12 + 1)
       + 8 * mp.quad(lambda t: t * (1 / t - 2 + t - t ** 2 + 1 / (1 - t)) * mp.log(t), [0, 1])
       + 4 * (mp.quad(lambda t: (1 - 2 * t + t ** 2 - t ** 3) * mp.log(1 - t), [0, 1]) + 1)
       + mp.mpf(67) / 9 + 8 * mp.zeta(2))
print("(B'') numerical check of the same energy sum:", mp.nstr(num, 15), " 73/3 =", mp.nstr(mp.mpf(73) / 3, 15),
      " (symbolic result below may carry a spurious imaginary branch term from sympy's apart/log handling)")
print("(B'') int x c_{phi,g}^{(1)} dx (AMV, N_f=0, a_s units, per C_A) =", sp.simplify(mom),
      " [known H->gg NLO in the effective theory: (alpha_s/pi) 73/4 = a_s C_A 73/3 for N_c=3]")
print("      contribution of the ln x^2 (number-density Jacobian) term, 8 int x p_gg ln x =",
      sp.simplify(8 * (Ireg(x * pg_reg, sp.log(x)) + Ireg(x / (1 - x), sp.log(x)))))
