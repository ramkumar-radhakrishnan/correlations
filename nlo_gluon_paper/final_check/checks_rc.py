"""Symbolic checks for the sections on the log(vee/Lambda) convention, the scale dependence and the running coupling.
-> checks_rc.json"""
import json
import sympy as sp

xi, lam, a, L1, L2, Lv, Lk = sp.symbols('xi lambda alpha_s ell_X ell_Xp L_vL L_k', positive=True)
Nc, b0 = sp.symbols('N_c beta_0', positive=True)
out = {}

# 1. Row V longitudinal weight = the paper's integrand; double-plus constant; exact endpoint form
Phi = xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi
paper = (xi**2 + (1 - xi)**2 + xi**2 * (1 - xi)**2) / (xi * (1 - xi))
out['RowV_weight_equals_paper_integrand'] = str(sp.simplify(Phi - paper))
dplus = sp.integrate(sp.simplify(Phi - 1 / xi - 1 / (1 - xi)), (xi, 0, 1))
out['double_plus_constant'] = str(sp.nsimplify(dplus))
exact = sp.integrate(Phi, (xi, lam, 1 - lam))
target = 2 * sp.log((1 - lam) / lam) - sp.Rational(11, 6) + 4 * lam - lam**2 + sp.Rational(2, 3) * lam**3
out['endpoint_form_minus_exact_numeric'] = max(abs(float(sp.re(sp.N((exact - target).subs(lam, v), 30)))) for v in (sp.Rational(1, 10), sp.Rational(1, 1000), sp.Rational(3, 7)))

# 2. the three ways of writing the LL part with log(vee/Lambda);  L_vL = log(vee/Lambda), L_k = log(vee/k+)
J, B = sp.symbols('J B')
ll = (Lv - Lk) * J + Lk * B                  # log(k+/Lambda) J + log(vee/k+) B
A_ = Lv * J + Lk * (B - J)
B_ = Lv * B + (Lv - Lk) * (J - B)
C_ = Lv * (J + B) / 2 + (Lv - 2 * Lk) * (J - B) / 2     # log(k+^2/(vee Lambda)) = L_vL - 2 L_k
out['A_minus_LL'] = str(sp.simplify(A_ - ll))
out['B_minus_LL'] = str(sp.simplify(B_ - ll))
out['C_minus_LL'] = str(sp.simplify(C_ - ll))
out['C_at_mid_rapidity'] = str(sp.simplify(C_.subs(Lk, Lv / 2)))

# 3. running coupling: sqrt(alpha(c0^2/X^2) alpha(c0^2/X'^2)) to O(alpha^2), one-loop running
aX = a / (1 - a * b0 / (4 * sp.pi) * L1)
aXp = a / (1 - a * b0 / (4 * sp.pi) * L2)
ser = sp.series(sp.sqrt(aX * aXp), a, 0, 3).removeO()
out['sqrt_aX_aXp_O(a^2)'] = str(sp.simplify(sp.expand(ser - a)))
# Row II + Row V beta0 terms relative to LO (both orderings), and their split
rowII = a * b0 / (4 * sp.pi) * (L1 + L2)
rowV = -a * b0 / (4 * sp.pi) * (L1 + L2)
out['RowII_plus_RowV'] = str(sp.simplify(rowII + rowV))
rc = sp.expand(ser / a - 1)
out['running_part_minus_half_RowII'] = str(sp.simplify(rc - rowII / 2))

# 4. Row II relative to LO = K_JSJ/K_LO of arXiv:2308.15545, Eq. (7) (Y -> X')
g = sp.symbols('g', positive=True)
perord = (g**4 * Nc / (32 * sp.pi**4)) / (g**2 / (2 * sp.pi**2))      # Row II prefactor / LO prefactor
ratio = sp.simplify(perord * 2 / (a / (2 * sp.pi)) .subs(a, g**2 / (4 * sp.pi)))
out['RowII_both_orderings_over_LO_in_units_alpha_over_2pi_times_Nc'] = str(sp.simplify(ratio / Nc))
# K_JSJ/K_LO = [alpha^2/(4 pi^3)] / [alpha/(2 pi^2)] = alpha/(2 pi)

# 5. asymptotics of the full logarithms (A=(x-z)^2, B=(x-w)^2)
A, Bb, kp, vee, Lam = sp.symbols('A B kp vee Lam', positive=True)
f1 = sp.log(((vee - kp) * A + kp * Bb) / (Lam * A + kp * Bb)) - (sp.log(vee / kp) + sp.log(A / Bb))
out['full_log_I_minus_asymptote'] = str(sp.limit(sp.limit(f1, vee, sp.oo), Lam, 0))
f2 = sp.log((kp * Bb - vee * A) / (kp * Bb - (kp + Lam) * A)) - (sp.log(vee / kp) + sp.log(A / (A - Bb)))
out['full_log_III_minus_asymptote'] = str(sp.simplify(sp.limit(sp.limit(f2, vee, sp.oo), Lam, 0)))
json.dump(out, open('checks_rc.json', 'w'), indent=1)
for k, v in out.items():
    print(k, ':', v)
