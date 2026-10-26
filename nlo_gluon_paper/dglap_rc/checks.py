"""Pole, beta_0 and constant bookkeeping of the remaining terms (Finite_terms_split.pdf, Secs. 7-9), in units of LO,
and the AABKL constants.

Remaining terms that contain poles or beta_0 (all others are finite):
  G1 Row II part 1 + c.c. = (as N_c/4pi) LO[(E_X+E_X')(11/3 + l1)] + (as N_c/2pi)(67/9 - pi^2/3 - l1^2 - Lv l1) LO
  G1 Row V 2b UV + c.c.   = (as N_c/2pi) LO[(E_X+E_X')(-11/6 + 2 lam)],   lam = log(1 - Lambda/k+)
  G2 Row IV, lines 1+2    = (as/2pi) int_a^1 dxi/xi^2 P^(r)_{gg,+}(xi) LO_xi[-E_Delta]   (DG = the integral with LO_xi[1],
                            DGl = with LO_xi[l_Delta])
  LO coupling counterterm = -(as beta0/4pi)(2/eps) LO
E_r = 2/eps + l_r,  l_r = log(r^2 mu^2/c0^2),  c0 = 2 e^{-gamma_E},  l1 = log(1-a),  a = k+/vee.
(The units LO[E_X+E_X'] of the rows were checked against your file in ../finite_v2, final_checks.json.)
Output: checks.json"""
import json
import sympy as sp

a, eps, Nc, al = sp.symbols('a epsilon N_c alpha_s', positive=True)
lX, lXp, lD, lam, DG, DGl = sp.symbols('ell_X ell_Xp ell_Delta lambda DG DGl')
Lv, l1 = -sp.log(a), sp.log(1 - a)
EX, EXp = 2 / eps + lX, 2 / eps + lXp
b0 = sp.Rational(11, 3) * Nc
out = {}

rowII_b0 = al / (4 * sp.pi) * b0 * (EX + EXp)                       # beta_0 part of G1 Row II part 1
rowV_b0 = al * Nc / (2 * sp.pi) * (EX + EXp) * (-sp.Rational(11, 6))  # beta_0 part of G1 Row V 2b UV
out['beta0: G1 II part 1 + G1 V 2b UV'] = str(sp.simplify(rowII_b0 + rowV_b0))

rowII = al * Nc / (4 * sp.pi) * ((EX + EXp) * (sp.Rational(11, 3) + l1) + 2 * (sp.Rational(67, 9) - sp.pi**2 / 3 - l1**2 - Lv * l1))
rowV = al * Nc / (2 * sp.pi) * (EX + EXp) * (-sp.Rational(11, 6) + 2 * lam)
rowIV = al / (2 * sp.pi) * (-(2 / eps) * DG - DGl)
ct = -al * b0 / (4 * sp.pi) * 2 / eps
tot = sp.expand(rowII + rowV + rowIV + ct)
pole = sp.expand(tot.coeff(eps, -1))
fin = sp.expand(tot - pole / eps)
# expected: pole = (as/2pi)(-2)[DG + (11 N_c/6) - N_c l1] + (as N_c/2pi) 8 lam   (lam -> 0 for Lambda -> 0)
exp_pole = al / (2 * sp.pi) * (-2) * (DG + sp.Rational(11, 6) * Nc - Nc * l1) + al * Nc / (2 * sp.pi) * 8 * lam
out['pole - [(as/2pi)(-2/eps)(P_gg x LO - N_c log(1-a) LO) + O(Lambda)]'] = str(sp.simplify(pole - exp_pole))
exp_fin = (al * Nc / (2 * sp.pi) * (sp.Rational(67, 9) - sp.pi**2 / 3)               # constant (CMW combination)
           - al * Nc / (2 * sp.pi) * (l1**2 + Lv * l1)                                # V_vee
           + al * Nc / (4 * sp.pi) * l1 * (lX + lXp)                                   # finite partner of the Q1 pole
           - al / (2 * sp.pi) * DGl                                                    # DGLAP finite term
           + al * Nc / (2 * sp.pi) * 2 * lam * (lX + lXp))                             # O(Lambda/k+)
out['finite - [CMW + V_vee + Q1 finite + DGLAP + O(Lambda)]'] = str(sp.simplify(fin - exp_fin))
out['coefficient of ell_X in the finite part minus (as N_c/4pi) log(1-a)  (no beta0 log left)'] = str(sp.simplify(sp.diff(fin, lX).subs(lam, 0) - al * Nc / (4 * sp.pi) * l1))

# AABKL: virtual integral, splitting function, scheme constant
x = sp.symbols('xi', positive=True)
f = sp.Function('f')
# int_0^1 dxi [(1-xi)/[xi]_+ + xi/[1-xi]_+ + xi(1-xi)] acting on 1:  [1/xi]_+ g = (g(xi)-g(0))/xi
v = sp.integrate(((1 - x) - 1) / x + (x - 1) / (1 - x) + x * (1 - x), (x, 0, 1))
out['AABKL virtual: int_0^1 dxi [(1-xi)/[xi]_+ + xi/[1-xi]_+ + xi(1-xi)]'] = str(v)
out['  times 2 C_A, = -(beta0/2)*2 ?  (-11/6 * 2N_c = -beta0)'] = str(sp.simplify(2 * Nc * v + b0))
PA = lambda xi: (1 - xi) / xi + xi / (1 - xi) + xi * (1 - xi)                   # AABKL, xi = unobserved fraction
Pours = lambda xi: 2 * Nc * (xi / (1 - xi) + (1 - xi) / xi + xi * (1 - xi))     # yours (regular part), xi = observed
out['P_ours(xi) - 2N_c P_AABKL(1-xi)'] = str(sp.simplify(Pours(x) - 2 * Nc * PA(1 - x)))
R = 2 * sp.exp(sp.polygamma(0, 1))
c0 = 2 * sp.exp(-sp.EulerGamma)
out['AABKL R = 2 e^{Psi(1)} vs your c0 = 2 e^{-gamma_E}'] = dict(R=float(R), c0=float(c0), R2=float(R**2), difference=str(sp.simplify(R - c0)))
for k, val in out.items():
    print('%-75s %s' % (k, val))
json.dump(out, open('checks.json', 'w'), indent=1)
