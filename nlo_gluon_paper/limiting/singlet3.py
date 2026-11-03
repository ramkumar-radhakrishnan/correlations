r"""Gaussian projectile: three-rho rows -> two-rho parts (pair rule), and colour-singlet (trace) form of every finite
term of Limiting_cases.tex.  Uses the machinery of ../finite_singlet (project.py: projection to adjoint traces,
dcombine, numerical checks; gsrc.py: general-U colour operators of the rows).

New compared with ../finite_singlet: the y' -> z subtraction of Group III (line 2 of Rows III(1), V, VI and of
Row III part 2) needs the colour operators with U(y') -> U(z): families 'famCz' and 'III2z'.  G1 V 2a + 2b finite
are combined (colour 'V2ab').
Checks: trace form = direct contraction (general SU(3) adjoint U); pair rule = reordering engine; the y' -> z
operators = the original operators with U(y') replaced numerically.
  python3 singlet3.py  -> singlet3.pkl, singlet3.json"""
import copy
import json
import os
import pickle
import sys
from collections import OrderedDict, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FS = os.path.normpath(os.path.join(HERE, '..', 'finite_singlet'))
sys.path.insert(0, FS)
os.chdir(FS)
import env  # noqa: E402,F401
import project as PJ  # noqa: E402
import gsrc  # noqa: E402
import compact as CP  # noqa: E402
import numerics as nm  # noqa: E402
from core import term as cterm, symmetrize  # noqa: E402
os.chdir(HERE)


def zsub(src):
    """the operator with every U(y') -> U(z)"""
    g = copy.deepcopy(src)
    g.subterms = [(s, [('U', u[1], u[2], 'z') if (u[0] == 'U' and u[3] == "y'") else u for u in col], w)
                  for s, col, w in src.subterms]
    return g


SRCS = {}
for key in ('GI.II', 'GI.III', 'GI.IV', 'GI.V1', 'GII.I1', 'GIII.I', 'GIII.II', 'GIII.III', 'GIII.V', 'GIII.VI'):
    SRCS[key] = gsrc.general(key)
    SRCS[key + 'z'] = zsub(SRCS[key])

FAM = OrderedDict([
    ('II', [('GI.II', 1)]),
    ('III', [('GI.III', 1)]),
    ('V1', [('GI.V1', 1)]),
    ('famIV', [('GI.IV', 1), ('GIII.I', -1), ('GIII.II', -1), ('GII.I1', -0.5)]),
    ('famC', [('GIII.III', 1), ('GIII.V', 2), ('GIII.VI', 2)]),
    ('famCz', [('GIII.IIIz', 1), ('GIII.Vz', 2), ('GIII.VIz', 2)]),
])


def groups(name):
    acc = defaultdict(lambda: defaultdict(complex))
    for key, wgt in FAM[name]:
        for pair, d in CP.pair_terms(SRCS[key]).items():
            for k, c in d.items():
                acc[pair][k] += wgt * c
    out = OrderedDict()
    for pair in sorted(acc, key=lambda p: (PJ.PORDER.index(p[0]), PJ.PORDER.index(p[1]))):
        ts = [PJ.T(c, [(u[1], u[2], u[3]) for u in Us], fs, rs) for (fs, Us, rs), c in acc[pair].items() if abs(c) > 1e-12]
        if ts:
            out[pair] = ts
    return out


def engine_check(name, seeds=(1, 2)):
    res = []
    for s in seeds:
        rng = np.random.default_rng(100 + s)
        U = {p: nm.random_adjoint(3, PJ.TS, rng, False) for p in PJ.PORDER}
        r = {p: rng.normal(size=8) for p in PJ.PORDER}
        mine = 0
        for pair, ts in groups(name).items():
            for t in ts:
                fl = [('U', i, j, p) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + [('rho', i, p) for i, p in t['rs']]
                mine += t['c'] * nm.eval_colour(fl, PJ.FF, U, r, 8)
        eng, scale = 0, 0
        for key, wgt in FAM[name]:
            src = SRCS[key]
            for sign, colour, word in src.subterms:
                Us = [(u[1], u[2], u[3]) for u in colour if u[0] == 'U']
                fs = [tuple(u[1:]) for u in colour if u[0] == 'f']
                labs = list(src.rho_labels)
                for w in CP.ordered(src.form, word):
                    for t in symmetrize(cterm(sign, Us, fs, w, [('TAG',) + tuple(labs)], key), (2,)):
                        sig = t['kern'][0][1:]
                        pair = [(labs[i], sig[i]) for i in range(3) if labs[i] != sig[i]][0]
                        keep = min(pair, key=labs.index)
                        mp = lambda p: keep if p in pair else p  # noqa: E731
                        fl = ([('U', i, j, mp(p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']]
                              + [('rho', i, mp(p)) for i, p in t['rs']])
                        v = wgt * t['c'] * 3 ** t.get('ncp', 0) * nm.eval_colour(fl, PJ.FF, U, r, 8)
                        eng += v
                        scale += abs(v)
        res.append(abs(mine - eng) / scale)
    return float(max(res))


def zsub_check(seed=4):
    """classical value of the z-substituted operator = original operator evaluated with U(y') := U(z)"""
    c = gsrc.C(seed, symmetric=False)
    out = {}
    for key in ('GIII.III', 'GIII.V', 'GIII.VI'):
        a = gsrc.classical(SRCS[key + 'z'], c)
        c2 = copy.deepcopy(c)
        c2.U['yp'] = c2.U['z']
        b = gsrc.classical(SRCS[key], c2)
        out[key] = float(abs(a - b) / abs(b))
    return out


# two-rho rows
def subst(terms, old, new):
    return [dict(t, Us=[(i, j, new if p == old else p) for i, j, p in t['Us']]) for t in terms]


TWO = OrderedDict([
    ('LO', PJ.TWO['LO']),
    ('H4', subst(subst(PJ.TWO['LO'], "y'", 'z'), 'y', 'z')),
    ('V2ab', PJ.TWO['V2a'] + PJ.TWO['V2b1']),
    ('V2b2', PJ.TWO['V2b2']),
    ('G2I2', PJ.TWO['G2I2']),
    ('III2', PJ.TWO['III2']),
    ('III2z', subst(PJ.TWO['III2'], "y'", 'z')),
])

if __name__ == '__main__':
    store = OrderedDict(two=OrderedDict(), three=OrderedDict())
    rep = OrderedDict(two=OrderedDict(), three=OrderedDict())
    for name, terms in TWO.items():
        d = PJ.project(terms)
        dc = PJ.dcombine(d)
        chk = PJ.check(terms, d, dc)
        store['two'][name] = dict(raw=d, comb=dc)
        rep['two'][name] = dict(n_raw=len(d), n_comb=len(dc), traces_vs_direct=chk[0], combined_vs_direct=chk[1])
        print(name, rep['two'][name], flush=True)
    for name in FAM:
        g_store, g_rep = OrderedDict(), OrderedDict()
        for pair, terms in groups(name).items():
            d = PJ.project(terms)
            dc = PJ.dcombine(d)
            chk = PJ.check(terms, d, dc)
            g_store[pair] = dict(raw=d, comb=dc)
            g_rep['%s<-%s' % pair] = dict(n_terms=len(terms), n_comb=len(dc), traces_vs_direct=chk[0], combined_vs_direct=chk[1])
        store['three'][name] = g_store
        rep['three'][name] = dict(groups=g_rep, pair_rule_vs_engine=engine_check(name))
        print(name, json.dumps(rep['three'][name]), flush=True)
    rep['zsub_check'] = zsub_check()
    print('zsub', rep['zsub_check'])
    pickle.dump(store, open(os.path.join(HERE, 'singlet3.pkl'), 'wb'))
    json.dump(rep, open(os.path.join(HERE, 'singlet3.json'), 'w'), indent=1)
