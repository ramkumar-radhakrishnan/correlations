"""Checks of the new forms in Finite_terms_simplified.tex  ->  final_checks.json
 1. closed forms vs the expressions of the file (re-run of ../finite/simp_exact.check, vee/k+ = 1.3 ... 300)
 2. Group II Row I part 2: the collected form (coefficients of lambda, lambda', Lv, log(1-a), 1) vs the pieces 2a-2f
 3. compact (classical) colour factors of the three-rho rows and families vs the colour factors as written
 4. Group I Row IV + Group III Row I, Group II Row I(A) + Group III Row II (classical colour)
 5. sympy: poles, beta_0, running coupling, the finite E-hat terms (exact in a), with the new scheme Lv -> L*
 6. split.json, compact.json, decomp.json are read in"""
import json
import numpy as np
import sympy as sp
import env  # noqa: F401
import fp
import ft
import simp
import simp_exact as SE
from fp import F, G, PI, sq

NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
out = {}

# ---------------------------------------------------------------- 1
print('1. closed forms vs file')
out['closed_forms_vs_file'] = SE.check()

# ---------------------------------------------------------------- 2
print('2. G2-I(2) collected')


def g2_i2_collected(c, kp, V, d=2):
    a, Lv, l1a = SE.av(kp, V)
    br, (A, B, Ap, Bp, zw2, zwp2) = simp.g2_i2_groups(c)
    D = Bp * A - B * Ap
    lam, lamp = np.log((1 - a) * A / B + a), np.log((1 - a) * Ap / Bp + a)
    q = d / 4
    cl = (br['a'] - br['b'] * A / B - br['c'] * B / A + (A * br['d'] - B * br['e']) / (A - B)
          - q * A * B * (Ap - Bp) / ((A - B) * zw2 * zwp2)) / D
    clp = (-br['a'] + br['b'] * Ap / Bp + br['c'] * Bp / Ap + (Ap * br['d'] - Bp * br['e']) / (Bp - Ap)
           + q * Ap * Bp * (A - B) / ((Ap - Bp) * zw2 * zwp2)) / D
    cLv = -br['b'] / (B * Bp) + br['c'] / (A * Ap)
    cl1 = br['b'] / (B * Bp)
    c0 = -q * (1 - a) / (zw2 * zwp2)
    pre = G**4 / (8 * PI**5) / kp
    return pre * (cl * lam + clp * lamp + cLv * Lv + cl1 * l1a + c0) * ft.phase(c, 'w', 'wp') / ft.PI3 * ft.col_G2_I2(c)


r2 = []
for s in (1, 2, 3, 4):
    c = fp.Cfg(s, sym=True, names=NAMES)
    for V in (1.3, 2.5, 10.0, 300.0):
        a_ = g2_i2_collected(c, 1.0, V)
        b_ = sum(SE.g2_i2(c, 1.0, V).values())
        r2.append(float(abs(a_ - b_) / abs(b_)))
out['G2_I2_collected'] = max(r2)
print('   max rel. diff. %.1e' % max(r2))

# ---------------------------------------------------------------- 3, 4
print('3./4. compact colour factors')


def fx(c, A, B, C):
    return np.einsum('abc,a,b,c->', F, A, B, C)


comp = {}
for name in ('G1-II', 'G1-III', 'G1-V(1)', 'K_IV family', 'K_C family', 'G1IV-G3I', 'G2A+G3II'):
    comp[name] = []
for s in (1, 2, 3, 4, 5):
    c = fp.Cfg(s, sym=True, names=NAMES)
    U, r = c.U, c.rho
    CL = ft.CL(c)
    # G1 II: CL^a [U^{ah}(w) f^{hbc}{rho^b(x),rho^c(y)} - f^{abc} U^{bd}(x)U^{ce}(y){rho^d(x),rho^e(y)}]
    a2, b2 = ft.col_G1_II_2a(c)
    comp['G1-II'].append(abs((a2 - b2) - 2 * (fx(c, CL, U['w'] @ r['x'], U['w'] @ r['y']) - fx(c, CL, U['x'] @ r['x'], U['y'] @ r['y'])))
                         / abs(a2 - b2))
    # G1 III: f^{abd} CL^a [(U(z)-U(y)) rho(y)]^b [U(x) rho(x)]^d
    comp['G1-III'].append(abs(ft.col_G1_III(c) - fx(c, CL, (U['z'] - U['y']) @ r['y'], U['x'] @ r['x'])) / abs(ft.col_G1_III(c)))
    # G1 V(1): 2 f^{abc} CL^a [(U(x) rho(y'))^b (U(z) rho(y))^c - (U(y') rho(y'))^b (U(y) rho(y))^c]
    v1 = ft.col_G1_V1(c)
    comp['G1-V(1)'].append(abs(v1 - 2 * (fx(c, CL, U['x'] @ r['yp'], U['z'] @ r['y']) - fx(c, CL, U['yp'] @ r['yp'], U['y'] @ r['y'])))
                           / abs(v1))
    # K_IV family: G1IV - G3I - G3II - G2A/2 = f^{acd} CL^a [(U(w)rho(y))^c (U(w)rho(x))^d - (U(z)rho(y))^c (U(x)rho(x))^d]
    fam = ft.col_G1_IV(c) - ft.col_G3_I(c) - ft.col_G3_II(c) - 0.5 * ft.col_G2_I1(c, 'w', 'wp')
    cmp_ = fx(c, CL, U['w'] @ r['y'], U['w'] @ r['x']) - fx(c, CL, U['z'] @ r['y'], U['x'] @ r['x'])
    comp['K_IV family'].append(abs(fam - cmp_) / abs(cmp_))
    comp['G1IV-G3I'].append(abs(ft.col_G1_IV(c) - ft.col_G3_I(c)) / abs(ft.col_G1_IV(c)))
    comb = np.sum(np.outer(CL, r['y']) * ft.Rone(c, 'w', 'z', 'x'))
    comp['G2A+G3II'].append(abs(-ft.col_G3_II(c) - 0.5 * ft.col_G2_I1(c, 'w', 'wp') - comb) / abs(comb))
    # K_C family: III(1) + 2 V + 2 VI = 2 f^{abc} [(U(y')-U(x'))rho(x')]^a [(U(z)-U(x))rho(x)]^b [(U(w)-U(y))rho(y)]^c
    famC = ft.col_G3_III1(c) + 2 * ft.col_G3_V(c) + 2 * ft.col_G3_VI(c)
    cmpC = 2 * fx(c, (U['yp'] - U['xp']) @ r['xp'], (U['z'] - U['x']) @ r['x'], (U['w'] - U['y']) @ r['y'])
    comp['K_C family'].append(abs(famC - cmpC) / abs(cmpC))
out['compact_colour'] = {k: float(max(v)) for k, v in comp.items()}
for k, v in out['compact_colour'].items():
    print('   %-12s %.1e' % (k, v))

# ---------------------------------------------------------------- 5
print('5. sympy tallies')
a, eps, Nc, al = sp.symbols('a epsilon N_c alpha_s', positive=True)
lX, lXp, lD = sp.symbols('ell_X ell_Xp ell_Delta')
Lv, l1 = -sp.log(a), sp.log(1 - a)
EX, EXp, ED = 2 / eps + lX, 2 / eps + lXp, 2 / eps + lD
b0 = sp.Rational(11, 3) * Nc
rowII = al * Nc / (4 * sp.pi) * ((EX + EXp) * (sp.Rational(11, 3) + 2 * Lv + l1)
                                 + 2 * (sp.Rational(67, 9) - sp.pi**2 / 3 - l1**2 - Lv * l1))
rowV = al * Nc / (2 * sp.pi) * (EX + EXp) * (-2 * Lv - sp.Rational(11, 6))
rowIV = al * Nc / sp.pi * Lv * ED
charge = -al / (2 * sp.pi) * (b0 / 2) * 2 / eps
tot = sp.expand(rowII + rowV + rowIV + charge)
pole = sp.simplify(tot.coeff(eps, -1) / (al / (2 * sp.pi) * (-2) * Nc))
fin = sp.expand(tot - tot.coeff(eps, -1) / eps)
target = (al * Nc / (2 * sp.pi) * Lv * (2 * lD - lX - lXp) + al * Nc / (4 * sp.pi) * l1 * (lX + lXp)
          + al * Nc / (2 * sp.pi) * (sp.Rational(67, 9) - sp.pi**2 / 3 - l1**2 - Lv * l1))
# the large-log part of these rows: coefficient of Lv at a -> 0 (with l1 -> 0)
Lsym = sp.Symbol('L')
largelog = sp.simplify(sp.limit(sp.diff(tot.subs(sp.log(a), -Lsym), Lsym).subs(sp.log(1 - a), 0), a, 0))
out['sympy'] = dict(pole_units_as2pi_m2eps_Nc=str(sp.simplify(sp.expand_log(pole, force=True))),
                    beta0_rows_II_V=str(sp.simplify(sp.expand(rowII + rowV).coeff(eps, -1) - al * Nc / sp.pi * (l1 - 2 * Lv))),
                    finite_minus_expected=str(sp.simplify(fin - target)),
                    coefficient_of_Lv=str(sp.factor(largelog)))
for k, v in out['sympy'].items():
    print('   %s: %s' % (k, v))

# ---------------------------------------------------------------- 6
for f in ('split.json', 'compact.json', 'decomp.json'):
    out[f.split('.')[0]] = json.load(open(f))
json.dump(out, open('final_checks.json', 'w'), indent=1)
