"""The engine's three-rho terms (../bfkl_final/rows.py) written as LaTeX and fed to reorder.py, against reorder.py on
your first form (Group III Row II, Group III Row I, Group I Row IV, Group II Row I): is the difference in the input?"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder as RO
import cmp_first as C
import cmp_reorder as CR
from rows_user import rows
T3, TE = C.T3, C.TE
sys.path.insert(0, os.path.join(HERE, '..', 'bfkl_final'))
import two_total  # noqa  (paths)
import rows as R
J = {lab: t for lab, t, cc in T3.row_terms('J', True)}
U = {r[3]: r[4] for r in rows('BFKL.tex')}
TP = {"x'": "x'", "w'": "w'", "y'": "y'"}
UNIT = r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,e^{-ik(w'-w)}\,"


def kk(a, b, c, d):
    return r"\frac{(%s-%s)\cdot(%s-%s)}{(%s-%s)^{2}\,(%s-%s)^{2}}" % (a, b, c, d, a, b, c, d)


def tex(t):
    k = ''.join(kk(*f[1:]) if f[0] == 'KK' else (_ for _ in ()).throw(ValueError(f)) for f in t['kern'])
    col = ''.join(r"U^{%s%s}(%s)\," % (i, j, p) for i, j, p in t['Us']) + ''.join(r"f^{%s%s%s}\," % f for f in t['fs'])
    col += ''.join(r"\rho^{%s}(%s)\," % (i, p) for i, p in t['rs'])
    return UNIT + k + r"\," + col


def rem_engine(ts, c):
    tot = 0
    for t in ts:
        coef = t['c'] * 3 ** t.get('ncp', 0)
        tot += coef * RO.remainder(tex(t), c)
    return tot


SETS = {'G3 II': ('G3 II', [408], ['G3 II #2'], True), 'G3 I': ('G3 I', [359], ['G3 I #2'], True),
        'G1 IV': ('G1 IV', [188], ['G1 IV #2'], True), 'G2 I first part': ('G2 I(1)', [247, 253], ['G2 I(1) A #2', 'G2 I(1) B #2'], False)}
out = {}
for name, (eng, first, jl, cc) in SETS.items():
    ts = R.THREE[eng]()          # cc already included by build(cc=True)
    res = {}
    for s in (5,):
        c = T3.configs(s)
        e = rem_engine(ts, c) + sum(RO.remainder_cc(J[l], c, cc) for l in jl)
        u = sum(RO.remainder_cc(C.prep(C.plus(U[l])), c, cc) for l in first)
        res[s] = dict(engine_plus_J=str(e), first_form=str(u), rel=abs(e - u) / abs(u))
        print(name, res[s], flush=True)
    out[name] = res
json.dump(out, open(os.path.join(HERE, 'engine_tex.json'), 'w'), indent=1)
