"""Referee 3: exact moments used in ref3_collinear.tex (sympy).  Run: python3 -I c1_moments.py

xi = momentum fraction of the observed gluon in its parent.
P_gg(xi) = 2 Nc [ xi/(1-xi)_+ + (1-xi)/xi + xi(1-xi) ] + beta0/2 delta(1-xi),  beta0 = 11 Nc/3 - 2 Nf/3 (T_R=1/2).
"""
import sympy as sp

x = sp.Symbol('xi', positive=True)
Nc, Nf = sp.symbols('N_c N_f', positive=True)
TR = sp.Rational(1, 2)
z2 = sp.zeta(2)

Preg = 2 * Nc * (x / (1 - x) + (1 - x) / x + x * (1 - x))     # real part (no plus needed when f(1)=0)
Pqg = TR * (x ** 2 + (1 - x) ** 2)                             # g -> q, per flavour


def I(f):
    return sp.nsimplify(sp.simplify(sp.integrate(sp.expand(f), (x, 0, 1))), [sp.pi ** 2])


print("[1] momentum sum rule  int xi P_gg = 0 (Nf=0):",
      sp.simplify(sp.integrate(sp.expand(x * Preg - 2 * Nc * 1 / (1 - x)), (x, 0, 1)) + 2 * Nc * sp.integrate(x * x / (1 - x) - x * x / (1 - x), (x, 0, 1))))
# int_0^1 xi * [xi/(1-xi)]_+ = int (xi^2 - 1)/(1-xi) = -3/2 ; int xi[(1-xi)/xi + xi(1-xi)] = 1/2 + 1/12
msr = 2 * Nc * (sp.integrate(sp.cancel((x ** 2 - 1) / (1 - x)), (x, 0, 1)) + sp.integrate(1 - x + x ** 2 * (1 - x), (x, 0, 1))) + sp.Rational(11, 6) * Nc
print("    2Nc(-3/2 + 1/2 + 1/12) + 11Nc/6 =", sp.simplify(msr))

M = sp.integrate(sp.expand(x * Preg * 2 * sp.log(x)), (x, 0, 1))
M = sp.simplify(M)
print("[2] int xi P_gg ln xi^2  =", M, " = 4Nc(65/144 - zeta2)?", sp.simplify(M - 4 * Nc * (sp.Rational(65, 144) - z2)) == 0,
      "; Nc=3:", float(M.subs(Nc, 3)))
Mq = sp.simplify(sp.integrate(sp.expand(x * 2 * Nf * Pqg * 2 * sp.log(x)), (x, 0, 1)))
print("    with quarks: int xi [P_gg + 2 Nf P_qg] ln xi^2 = ", sp.simplify(M + Mq), " (Nc=3,Nf=3: %.6f)" % float((M + Mq).subs({Nc: 3, Nf: 3})))

# [3] Row II.IV collinear terms, momentum-weighted (inclusive over the collinear pair)
#   DGLAP (draft form) = -C int P F,  F = xi^-2 [ G_log(k/xi) - ln xi^2 Lhat(k/xi) ]
#   map xi^-2 G(k/xi) -> weight xi, G(Q):  -C [int xi P] G_log + C [int xi P ln xi^2] Lhat
#   relative to LO = N0 Lhat, C/N0 = alpha_s/(2 pi):
const_draft = M / 2                       # in units alpha_s/pi ... keep explicit below
a_over_4pi = sp.Symbol('alpha_s_Nc_over_4pi')
c_draft = sp.simplify(M * 2 / Nc)          # (alpha_s/2pi) M = (alpha_s Nc/4pi) * (2M/Nc)
print("[3] Row II.IV DGLAP, momentum-weighted, LO-proportional constant (draft convention):")
print("    (alpha_s Nc/4pi) x", sp.simplify(c_draft), "=", sp.nsimplify(sp.expand(c_draft), [sp.pi ** 2]), "=", float(c_draft))
print("    CXY convention: 0 (int xi P_gg = 0)")
target = -(sp.Rational(67, 9) - sp.pi ** 2 / 3)
print("    target of Sec. 7.6: (alpha_s Nc/4pi) x", target, "=", float(target))
print("    target - draft value = (alpha_s Nc/4pi) x", sp.simplify(target - c_draft), "=", float(target - c_draft))
cT = sp.integrate(x * 2 * x * (1 - x), (x, 0, 1))    # tensor term: 2Nc C int dxi (1-xi)/xi ... = 2Nc C int xi(1-xi) xi^-2 G_T
print("    tensor term (eq:app_IVfin): weight 2Nc int xi^2(1-xi) =", cT, "* Nc -> (alpha_s Nc/4pi) x", sp.simplify(2 * cT), "x G_T (not LO structure)")

# [4] KLSZ double-plus integrals, N = xi^2+(1-xi)^2+xi^2(1-xi)^2
N = x ** 2 + (1 - x) ** 2 + x ** 2 * (1 - x) ** 2
dp = sp.simplify(sp.integrate(sp.expand(sp.cancel((N - (1 - x) - x) / (x * (1 - x)))), (x, 0, 1)))
print("[4] double plus int N/(xi_+(1-xi)_+) =", dp)
f = sp.cancel(N / (x * (1 - x))) - 1 / x - 1 / (1 - x)
for name, g in (("ln xi", sp.log(x)), ("ln(1-xi)", sp.log(1 - x)), ("ln xi(1-xi)", sp.log(x) + sp.log(1 - x)), ("ln xi^2", 2 * sp.log(x))):
    val = sp.simplify(sp.integrate(sp.expand(f * g), (x, 0, 1)))
    print("    double plus x %-12s =" % name, sp.nsimplify(val, [sp.pi ** 2]), "=", float(val))

# [5] Mellin-type pieces
print("[5] pieces: int xi^2 ln xi/(1-xi) =", sp.simplify(sp.integrate(x ** 2 * sp.log(x) / (1 - x), (x, 0, 1))),
      "; int (1-xi) ln xi =", sp.integrate((1 - x) * sp.log(x), (x, 0, 1)),
      "; int xi^2(1-xi) ln xi =", sp.integrate(x ** 2 * (1 - x) * sp.log(x), (x, 0, 1)))
