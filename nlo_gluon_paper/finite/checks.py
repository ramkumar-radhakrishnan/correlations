"""All checks quoted in Finite_terms_final.tex -> checks.json.
 1. verify.json (row-by-row check of Finite_terms.tex, run verify.py) and simp_check.json (run simp.py) are read in
 2. log(vee/k+) coefficient of every closed-form row = (region P limit - region T limit)/2pi of the exact integrand
 3. two oscillatory p+ integrals in closed form (mpmath); the C rows at large vee are checked in crows_w.py
 4. classical (symmetrized) three-rho relations: G1-IV + G3-I = 0 (integrand, every p+), G2-I(1)A + G3-II
 5. sympy: the log(vee/k+) coefficients of 2d, 2e, 2f vanish; the pole / beta_0 / log(vee/k+) tallies
"""
import json
import numpy as np
import mpmath as mp
import sympy as sp
import env  # noqa: F401
import fp
import fp1
import fixed23 as fx
import user23 as us
import ft
import simp

NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
out = {}

# ---------------------------------------------------------------- 1
V = json.load(open('verify.json'))
out['verify_spread'] = {k: max(r['spread'] for r in v) for k, v in V.items()}
S = json.load(open('simp_check.json'))
out['simp_rel_at_1e6_1e9'] = {k: [max(e[0] for e in v), max(e[1] for e in v)] for k, v in S.items()}

# ---------------------------------------------------------------- 2
lim = {}
for seed in (1, 2, 3):
    c = fp.Cfg(seed, sym=True, names=NAMES)
    k = 1.0
    ph = ft.phase(c, 'w', 'wp')
    rows = {'G1-IV': (lambda p: fp1.row4(c, k, p), simp.g1_iv(c, k)[0]),
            'G2-I(1)': (lambda p: ph * (fx.g2_1_termA(c, k, p) + fx.g2_1_termB(c, k, p)), simp.g2_i1(c, k)[0]),
            'G2-I(2)': (lambda p: ph * us.g2_2_line62(c, k, p), sum(simp.g2_i2(c, k)[0].values())),
            'G3-I': (lambda p: ph * fx.g3_1_line44(c, k, p), simp.g3_i(c, k)[0]),
            'G3-II': (lambda p: ph * fx.g3_2_line55(c, k, p), simp.g3_ii(c, k)[0])}
    for n, (f, C) in rows.items():
        PT = (1e9 * f(1e9) - 1e-9 * f(1e-9)) / (2 * np.pi)
        lim[n] = max(lim.get(n, 0), float(abs(PT - C) / abs(C)))
out['Lv_coefficient_vs_P_minus_T'] = lim

# ---------------------------------------------------------------- 3
cf = []
for a, k in ((0.7, 1.0), (-1.3, 1.0), (2.5, 0.6)):
    J1 = mp.quadosc(lambda p: mp.exp(-1j * a * p) / (p + k), [0, mp.inf], omega=abs(a))
    J2 = (mp.quad(lambda p: (mp.exp(-1j * a * p) - 1) / p, [0, k])
          + mp.quadosc(lambda p: mp.exp(-1j * a * p) / p, [k, mp.inf], omega=abs(a)))
    cf.append(dict(a=a, k=k, J1=float(abs(J1 - mp.exp(1j * a * k) * mp.e1(1j * a * k))),
                   J2=float(abs(J2 - (-mp.euler - mp.log(1j * a * k))))))
out['closed_forms'] = cf


# ---------------------------------------------------------------- 4
cl = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = fp.Cfg(seed, sym=sym, names=NAMES)
        for p in (0.05, 0.7, 3.0, 40.0):
            g4 = fp1.row4(c, 1.0, p)
            g31 = ft.phase(c, 'w', 'wp') * fx.g3_1_line44(c, 1.0, p)
            cl.append(float(abs(g4 + g31) / abs(g4)))
out['G1IV_plus_G3I_integrand_max_rel'] = max(cl)
comb = []
for seed in (1, 2, 3):
    c = fp.Cfg(seed, sym=True, names=NAMES)
    A, _ = ft.G2_I1(c, 1.0, 37.0, 1e-3)
    K = ft.G3_I_II_kernel(c, 1.0, 37.0, 1e-3)
    pred = ft.phase(c, 'w', 'wp') / ft.PI3 * K * (-np.sum(np.outer(ft.CL(c), c.rho['y']) * ft.Rone(c, 'w', 'z', 'x')))
    comb.append(float(abs(A + ft.G3_II(c, 1.0, 37.0, 1e-3) - pred) / abs(pred)))
out['G2I1A_plus_G3II_combined_max_rel'] = max(comb)
n2 = []
for seed in (1, 2, 3, 4):
    c = fp.Cfg(seed, sym=True, names=NAMES)
    KG3 = ft.G3_I_II_kernel(c, 1.0, 3000.0, 1e-6)
    KIII = sum(ft.G1_III(c, 1.0, 3000.0, 1e-6).values()) / (ft.phase(c, 'w', 'wp') / ft.PI3 * ft.col_G1_III(c))
    n2.append(float(abs(KG3 + KIII) / abs(KG3)))
out['N2_network_G1III_plus_G2I1G3II_rel'] = n2

# ---------------------------------------------------------------- 5
A_, B_, Ap, Bp = sp.symbols('A B Ap Bp', positive=True)
D = Bp * A_ - B_ * Ap
out['Lv_coef_2d'] = str(sp.simplify(-1 / ((A_ - B_) * (Bp - Ap)) + A_ / ((A_ - B_) * D) + Ap / ((Bp - Ap) * D)))
out['Lv_coef_2e'] = str(sp.simplify(1 / ((A_ - B_) * (Bp - Ap)) - B_ / ((A_ - B_) * D) - Bp / ((Bp - Ap) * D)))
out['Lv_coef_2f'] = str(sp.simplify((Bp * B_ - A_ * Ap) / ((A_ - B_) * (Ap - Bp)) - A_ * B_ * (Ap - Bp) / ((A_ - B_) * D)
                                    + Ap * Bp * (A_ - B_) / ((Ap - Bp) * D)))
out['Lv_coef_2b'] = str(sp.simplify(-A_ / (B_ * D) + Ap / (Bp * D) + 1 / (B_ * Bp)))
out['Lv_coef_2c'] = str(sp.simplify((Bp / Ap - B_ / A_) / D - 1 / (A_ * Ap)))
# tallies, relative to LO (alpha_s units), per ordering x 2 (c.c.):  E = 2/eps + ell
g, Nc, Lv, Ehat, l1 = sp.symbols('g N_c L_v Ehat l_1')
alpha = g**2 / (4 * sp.pi)
LOpref = g**2 / (2 * sp.pi**2)
rowII = 2 * (-(g**4 * Nc / (32 * sp.pi**4)) / LOpref) * (-Ehat * (sp.Rational(11, 3) + 2 * Lv + l1) - sp.Rational(67, 9) + sp.pi**2 / 3)
rowV = 2 * ((g**4 * Nc / (16 * sp.pi**4)) / LOpref) * (-2 * Lv - sp.Rational(11, 6)) * Ehat
rowIV_Lv = ((g**2 / (4 * sp.pi**2)) * Nc) * (-Ehat) * (-Lv)          # Group II Row IV, -L_v e^{-ik.(y-y)} W term; E_Delta -> Ehat for the pole tally
beta0 = sp.Rational(11, 3) * Nc
tot = sp.expand(rowII + rowV + rowIV_Lv)
out['tally_beta0_coefficient_of_Ehat'] = str(sp.simplify(tot.coeff(Ehat).subs({Lv: 0, l1: 0}) / alpha))
out['tally_Lv_coefficient_of_Ehat'] = str(sp.simplify(sp.diff(tot, Ehat, Lv) / alpha))
out['tally_constant'] = str(sp.simplify(tot.subs({Ehat: 0}) / (alpha / (2 * sp.pi) * Nc)))
out['rowII_beta0_part'] = str(sp.simplify(sp.expand(rowII).coeff(Ehat).subs({Lv: 0, l1: 0}) / (alpha * beta0 / (4 * sp.pi))))
out['rowV_beta0_part'] = str(sp.simplify(sp.expand(rowV).coeff(Ehat).subs({Lv: 0}) / (alpha * beta0 / (4 * sp.pi))))


def tojson(o):
    if isinstance(o, dict):
        return {k: tojson(v) for k, v in o.items()}
    if isinstance(o, list):
        return [tojson(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


json.dump(tojson(out), open('checks.json', 'w'), indent=1)
for k, v in out.items():
    print(k, ':', v if not isinstance(v, list) or len(str(v)) < 300 else str(v)[:300])
