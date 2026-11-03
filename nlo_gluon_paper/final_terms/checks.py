r"""Checks for Final_result.tex.   python3 checks.py -> checks.json"""
import difflib
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad

from limblocks import blocks, toks

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# 1. which equations of the new file are identical (token by token) to the verified exact-vee file
N = blocks(os.path.join(HERE, 'input', 'Other_terms.tex'))
A = blocks(os.path.join(HERE, '..', 'neglog', 'corrected', 'Other_terms.tex'))
ident, changed = [], []
for j, (n, b) in enumerate(N):
    cands = [i for i in range(len(A)) if A[i][0] == n]
    r = max([difflib.SequenceMatcher(None, toks(A[i][1]), toks(b), autojunk=False).ratio() for i in cands] or [0])
    (ident if r == 1.0 else changed).append('%d %s (%.3f)' % (j, n, r))
out['identical to verified'] = ident
out['new or changed'] = changed

# 2. Group I Row II part 1: the log(1-Lambda/k+) moved from 2b UV.  2b UV: g^4Nc/(16pi^4) * 2L * [2/eps..] * C(+) ;
#    in the Row II form g^4Nc/(32pi^4) * c * [-2/eps..] * C(-):  c = ?
L, c = sp.symbols('L c')
sol = sp.solve(sp.Eq(sp.Rational(1, 16) * 2 * L, sp.Rational(1, 32) * c), c)
out['coefficient of log(1-Lambda/k+) in Row II part 1 form'] = str(sol[0])
sol = sp.solve(sp.Eq(sp.Rational(1, 16) * sp.Rational(-11, 6), sp.Rational(1, 32) * c), c)
out['beta0 check: -11/6 of 2b UV in Row II form (must be -11/3, cancels +11/3)'] = str(sol[0])

# 3. poles: Row II part 1 + c.c. gives (alpha_s/2pi)(-2/eps)(-Nc log(1-a)) LO; P_gg with the standard plus
#    distribution on [a,1] has the endpoint term 2Nc log(1-a) delta(1-xi); the D_{h/g} counterterm absorbs exactly that.
a = sp.symbols('a', positive=True)
Nc = sp.symbols('N_c', positive=True)
out['pole left over (coefficient of (alpha_s/2pi)(-2/eps) LO)'] = str(sp.simplify(-Nc * sp.log(1 - a)))

# 4. Group III: your plus prescription on [0, vee-k+] + your log(1-a) line = standard plus distribution on [a,1]
kp, vee = 1.0, 37.0
av = kp / vee
F = lambda p: np.cos(0.8 * p / kp) / (1 + 0.3 * p / kp)                      # noqa: E731
G = lambda z: F(kp * (1 - z) / z)                                            # noqa: E731
yours = quad(lambda p: kp / (p * (p + kp)) * (F(p) - F(0.0)), 0, vee - kp, limit=400)[0] + np.log(1 - av) * F(0.0)
std = quad(lambda z: (G(z) - G(1.0)) / (1 - z), av, 1, limit=400)[0] + np.log(1 - av) * G(1.0)
out['G3: yours + log(1-a) line = int_a^1 [1/(1-zeta)]_+ (standard)'] = abs(yours - std)
out['G3: int_0^{vee-k+} dp+/(p++k+) = int_a^1 dzeta/zeta'] = abs(
    quad(lambda p: F(p) / (p + kp), 0, vee - kp, limit=400)[0] - quad(lambda z: G(z) / z, av, 1, limit=400)[0])

# 5. G3 III(2) phase-free line with the upper limit vee - k+
errs = []
for Aq, Bq in ((0.7, 1.9), (2.0, 0.3)):
    num = quad(lambda p: Bq / (Aq * (p + kp) * (p * Aq + kp * Bq)), 0, vee - kp)[0]
    cf = Bq * np.log(Bq / ((1 - av) * Aq + av * Bq)) / (kp * Aq * (Bq - Aq))
    errs.append(abs(num - cf) / abs(cf))
out['G3 III(2) closed form with exact a (rel)'] = max(errs)

# 6. prefactor bookkeeping: N = g^4/((2pi)^3 8pi^5 k+) = alpha_s^2/(4 pi^6 k+)
g, al = sp.symbols('g alpha_s', positive=True)
out['N = alpha_s^2/(4 pi^6 k+)'] = str(sp.simplify((g**4 / ((2 * sp.pi)**3 * 8 * sp.pi**5)).subs(g, sp.sqrt(4 * sp.pi * al)) - al**2 / (4 * sp.pi**6)))
out['NLO/LO prefactor: N k+ / (g^2/((2pi)^3 2pi^2)) = alpha_s/pi^2'] = str(sp.simplify(
    (g**4 / (8 * sp.pi**5)) / (g**2 / (2 * sp.pi**2)) / (al / sp.pi**2)).subs(g, sp.sqrt(4 * sp.pi * al)))

# 7. a -> 0 of every equation: from ../limiting/limcheck.json (difference ~ a)
lc = json.load(open(os.path.join(HERE, '..', 'limiting', 'limcheck.json')))
out['max rel. difference exact-a vs a->0 at a = 1e-2,1e-4,1e-6,1e-8'] = [max(r['reldiff'][i] for v in lc.values() for r in v) for i in range(4)]

json.dump(out, open(os.path.join(HERE, 'checks.json'), 'w'), indent=1, default=str)
for k, v in out.items():
    print(k, ':', v)
