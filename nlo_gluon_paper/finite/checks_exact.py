"""Checks for Finite_terms_exact.tex (exact in vee, Lambda -> 0)  ->  checks_exact.json
 1. sympy: the delta(1-xi) tally of the poles with log(1-a), a = k+/vee, kept exactly; beta_0; the finite E-hat terms
 2. the four-rho terms (Eq. (4rho) of 4rho.tex): symmetrized sum of the lines with upper limit vee - k+ alone
 3. simp_exact.json (run simp_exact.py) is read in"""
import itertools
import json
import os
import sys
from collections import defaultdict
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'three_rho'))
out = {}

# ---------------------------------------------------------------- 1
a, eps, Nc, al = sp.symbols('a epsilon N_c alpha_s', positive=True)
lX, lXp, lD = sp.symbols('ell_X ell_Xp ell_Delta')
Lv, l1 = -sp.log(a), sp.log(1 - a)
EX, EXp, ED = 2 / eps + lX, 2 / eps + lXp, 2 / eps + lD
b0 = sp.Rational(11, 3) * Nc
# in units of LO (inside the transverse integrals), both orderings
rowII = al * Nc / (4 * sp.pi) * ((EX + EXp) * (sp.Rational(11, 3) + 2 * Lv + l1)
                                 + 2 * (sp.Rational(67, 9) - sp.pi**2 / 3 - l1**2 - Lv * l1))
rowV = al * Nc / (2 * sp.pi) * (EX + EXp) * (-2 * Lv - sp.Rational(11, 6))
rowIV = al * Nc / sp.pi * Lv * ED            # delta(1-xi) part of G2 Row IV after the conversion to the standard plus
charge = -al / (2 * sp.pi) * (b0 / 2) * 2 / eps
tot = sp.expand(rowII + rowV + rowIV + charge)
pole = sp.simplify(tot.coeff(eps, -1) / (al / (2 * sp.pi) * (-2) * Nc))
out['pole_delta_coefficient_units_as2pi_m2eps_Nc'] = str(sp.simplify(sp.expand_log(pole, force=True)))
# Rows II + V: the 11/3 (beta_0) poles cancel; what is left is (alpha_s N_c/pi)(1/eps)(log(1-a) - 2 Lv)
out['beta0_in_rows_II_plus_V'] = str(sp.simplify(sp.expand(rowII + rowV).coeff(eps, -1) - al * Nc / sp.pi * (l1 - 2 * Lv)))
fin = sp.expand(tot - tot.coeff(eps, -1) / eps)
target = (al * Nc / (2 * sp.pi) * Lv * (2 * lD - lX - lXp) + al * Nc / (4 * sp.pi) * l1 * (lX + lXp)
          + al * Nc / (2 * sp.pi) * (sp.Rational(67, 9) - sp.pi**2 / 3 - l1**2 - Lv * l1))
out['finite_minus_expected'] = str(sp.simplify(fin - target))
out['remainder_pole'] = str(sp.simplify(tot.coeff(eps, -1) - (al / (2 * sp.pi)) * (-2) * Nc * sp.Rational(11, 6)))

# ---------------------------------------------------------------- 2
import numerics as nm          # noqa: E402
import fourrho                 # noqa: E402
from four_check import kk4, NC  # noqa: E402


def four_sum(labels, seed):
    ts, f = nm.structure_constants(NC)
    rng = np.random.default_rng(seed)
    ORD4 = ['w', "w'", 'z', 'p1', 'p2', 'p3', 'p4']
    Um = {l: nm.random_adjoint(NC, ts, rng, True) for l in ORD4}
    rv = {l: rng.normal(size=NC * NC - 1) for l in ORD4}
    acc = defaultdict(complex)
    scale = defaultdict(float)
    for c, row, lab, facs, word, kern, cc in fourrho.words():
        if lab.split('.')[0] not in labels:
            continue
        rl = [p for _, p in word]
        for perm in itertools.permutations(['p1', 'p2', 'p3', 'p4']):
            m = dict(zip(rl, perm))
            g = lambda p: m.get(p, p)
            parts, sgn = [], 1
            for (x1, x2), (x3, x4) in kern:
                k, s = kk4(g(x1), g(x2), g(x3), g(x4))
                parts.append(k)
                sgn *= s
            key = tuple(sorted(parts, key=str))
            ff = [('U', u[1], u[2], g(u[3])) for u in facs] + [('rho', i, g(p)) for i, p in word]
            v = c * sgn * nm.eval_colour(ff, f, Um, rv, NC * NC - 1) / 24
            acc[key] += v
            scale[key] += abs(v)
    return float(max(abs(v) for v in acc.values())), float(max(scale.values()))


VEEK = ['L3', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10']
four = {}
for name, labs in (('all lines', ['L%d' % i for i in range(1, 11)]), ('upper limit vee-k+ (L3, L5-L10)', VEEK)):
    four[name] = [four_sum(labs, s) for s in (3, 5)]
    print(name, four[name])
out['four_rho_symmetrized_max_sum_and_scale'] = four

# ---------------------------------------------------------------- 3
out['simp_exact'] = json.load(open(os.path.join(HERE, 'simp_exact.json')))
json.dump(out, open(os.path.join(HERE, 'checks_exact.json'), 'w'), indent=1)
for k in ('pole_delta_coefficient_units_as2pi_m2eps_Nc', 'beta0_in_rows_II_plus_V', 'finite_minus_expected', 'remainder_pole'):
    print(k, ':', out[k])
