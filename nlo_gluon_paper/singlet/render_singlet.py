"""All evolution terms (two rho's, coefficient of N log(vee/Lambda)) in colour-singlet form.

For every source (row), every term is projected (rho^a rho^b -> delta^{ab} <rho^f rho^f>/(N_c^2-1)), written as
adjoint traces (trace.py), the two charge points are named x, y (x attached to w where possible), and equal
(kernel, trace) monomials are combined.  Totals per kernel structure:
  * the sum of all rows, read with your convention U^{bc} = U^{cb} (U^dag -> U, U U -> 1), against JIMWLK x LO
    under the same convention (symbolic);
  * JIMWLK x LO for a general U (from the operator target), and BFKL x LO (from the region-P target).
Output: singlet.pkl (all data), singlet_checks.json."""
import json
import os
import pickle
from collections import OrderedDict, defaultdict
import numpy as np
import sympy as sp
import env  # noqa: F401
from core import kkey, Config
import sets as S2
import target as TG
from parse_user import sec2_terms
import trace as TR
from trace import Nc

HERE = os.path.dirname(os.path.abspath(__file__))
EVO = os.path.join(HERE, '..', 'final_check', 'input', 'Evolution_terms.tex')
ORDN = ['x', 'y', 'w', "w'", 'z']


# ---------------------------------------------------------------- canonical kernel
def vkey(a, b):
    return ((a, b), 1) if ORDN.index(a) < ORDN.index(b) else ((b, a), -1)


def fkey(f):
    a, b, c, d = f[1:]
    v1, s1 = vkey(a, b)
    v2, s2 = vkey(c, d)
    if v1 == v2:
        return ('inv2', v1), s1 * s2
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORDN.index(v[0]), ORDN.index(v[1])))), s1 * s2


def kern_key(kern):
    parts, s = [], 1
    for f in kern:
        assert f[0] == 'KK', f
        k, si = fkey(f)
        parts.append(k)
        s *= si
    return tuple(sorted(parts, key=str)), s


def charges_of(t):
    ch = []
    for _, p in t['rs']:
        if p not in ch:
            ch.append(p)
    return ch


def rename_term(t, mp):
    g = lambda p: mp.get(p, p)
    return dict(t, Us=[(i, j, g(p)) for i, j, p in t['Us']], rs=[(i, g(p)) for i, p in t['rs']],
                kern=[(f[0],) + tuple(g(x) for x in f[1:]) for f in t['kern']])


LOKEY = ('KK', ('x', 'w'), ('y', "w'"))


def monomial(t, sym=False):
    """-> list of (kernel key, trace key, coefficient) (one entry; empty if zero)"""
    ch = charges_of(t)
    assert len(ch) == 2, t
    best = None
    for mp in ({ch[0]: 'x', ch[1]: 'y'}, {ch[0]: 'y', ch[1]: 'x'}):
        tt = rename_term(t, mp)
        kk, ks = kern_key(tt['kern'])
        c, ws = TR.singlet(tt)
        if sym:
            c, ws = symconv(c, ws)
        cand = (0 if LOKEY in kk else 1, str(kk), str(ws), kk, ws, c * ks)
        if best is None or cand[:3] < best[:3]:
            best = cand
    _, _, _, kk, ws, c = best
    return kk, ws, sp.expand(c)


def symconv(c, ws):
    """your convention U^{bc} = U^{cb}: U^dag = U and U(p) U(p) = 1"""
    if c == 0:
        return c, ws
    ws2 = [[('U', l[1]) if l[0] == 'Ud' else l for l in w] for w in ws]
    out = []
    for w in ws2:
        changed = True
        while changed and w:
            changed = False
            n = len(w)
            for k in range(n):
                a, b = w[k], w[(k + 1) % n]
                if n >= 2 and a[0] == 'U' and b[0] == 'U' and a[1] == b[1]:
                    w = [] if n == 2 else [l for q, l in enumerate(w) if q not in {k, (k + 1) % n}]
                    changed = True
                    break
        out.append(w)
    c2, ws3 = TR.simplify(c, out)
    if ws3 is None:
        return sp.Integer(0), ()
    return TR.canon(c2, ws3)


def collect(terms, sym=False):
    acc = defaultdict(lambda: defaultdict(lambda: sp.Integer(0)))
    for t in terms:
        kk, ws, c = monomial(t, sym)
        if c == 0:
            continue
        acc[kk][ws] += c
    out = OrderedDict()
    for kk in sorted(acc, key=kernel_order):
        d = {ws: sp.factor(sp.expand(c)) for ws, c in acc[kk].items() if sp.expand(c) != 0}
        if d:
            out[kk] = OrderedDict(sorted(d.items(), key=lambda it: trace_order(it[0])))
    return out


def kernel_order(kk):
    lo = [f for f in kk if f[0] == 'KK' and 'z' not in (f[1][1], f[2][1])]
    return (0 if LOKEY in kk else 1, str(lo), str(kk))


def trace_order(ws):
    return (len(ws), sum(len(w) for w in ws), str(ws))


# ---------------------------------------------------------------- sources
def sources():
    A = S2.A_sets()
    S = OrderedDict()
    S[('Genuine', 'Group I, Row II')] = A['g:I.II']
    S[('Genuine', 'Group I, Row V')] = A['g:I.V']
    S[('Genuine', 'Group II, Row I')] = A['g:II.I']
    S[('Genuine', 'Group II, Row IV')] = A['g:II.IV']
    S[('Genuine', 'Group III, Row III')] = A['g:III.III']
    sec2 = sec2_terms(EVO)
    names = OrderedDict([
        ("From $4\\, \\rho's$||", ('From four $\\rho$\'s', 'your Sec.~2')),
        ('Group I (with c.c.)|Row II|', ('From three $\\rho$\'s', 'Group I, Row II')),
        ('Group I (with c.c.)|Row III|N1', ('From three $\\rho$\'s', 'Group I, Row III, $N_1$')),
        ('Group I (with c.c.)|Row III|N2', ('From three $\\rho$\'s', 'Group I, Row III, $N_2$')),
        ('Group I (with c.c.)|Row IV|N1', ('From three $\\rho$\'s', 'Group I, Row IV, $N_1$')),
        ('Group I (with c.c.)|Row IV|N2', ('From three $\\rho$\'s', 'Group I, Row IV, $N_2$')),
        ('Group I (with c.c.)|Row V|', ('From three $\\rho$\'s', 'Group I, Row V')),
        ('Group II|Row I|', ('From three $\\rho$\'s', 'Group II, Row I')),
        ('Group III\\, (with c.c.)|Row I|N1', ('From three $\\rho$\'s', 'Group III, Row I, $N_1$')),
        ('Group III\\, (with c.c.)|Row I|N2', ('From three $\\rho$\'s', 'Group III, Row I, $N_2$')),
        ('Group III\\, (with c.c.)|Row II|N1', ('From three $\\rho$\'s', 'Group III, Row II, $N_1$')),
        ('Group III\\, (with c.c.)|Row II|N2', ('From three $\\rho$\'s', 'Group III, Row II, $N_2$')),
        ('Group III\\, (with c.c.)|Row III|', ('From three $\\rho$\'s', 'Group III, Row III')),
        ('Group III\\, (with c.c.)|Row V|', ('From three $\\rho$\'s', 'Group III, Row V')),
        ('Group III\\, (with c.c.)|Row VI|', ('From three $\\rho$\'s', 'Group III, Row VI')),
    ])
    by = defaultdict(list)
    for lab, t in sec2:
        by[lab].append(t)
    assert set(by) == set(names), set(by) ^ set(names)
    S[('From four $\\rho$\'s', '\\texttt{4rho.tex} lists (missing in \\texttt{Evolution\\_terms.tex}, M1)')] = A['4rho']
    for lab, key in names.items():
        S[key] = by[lab]
    return S


# ---------------------------------------------------------------- comparison helpers
def diff(a, b):
    """a - b for two collected dicts -> dict of nonzero"""
    out = {}
    for kk in set(a) | set(b):
        for ws in set(a.get(kk, {})) | set(b.get(kk, {})):
            d = sp.expand(a.get(kk, {}).get(ws, 0) - b.get(kk, {}).get(ws, 0))
            if d != 0:
                out[(kk, ws)] = d
    return out


def numeric_value(coll, cfg):
    """sum over monomials at the points of cfg, symmetrized over (x, y) -> (p1, p2), (p2, p1); general U"""
    return 0.5 * (_numeric(coll, cfg, {'x': 'p1', 'y': 'p2'}) + _numeric(coll, cfg, {'x': 'p2', 'y': 'p1'}))


def _numeric(coll, cfg, mp):
    P = {k: v for k, v in cfg.P.items()}
    tot = 0.0
    for kk, d in coll.items():
        kv = 1.0
        for f in kk:
            if f[0] == 'inv2':
                (a, b), = f[1:]
                v = P[mp.get(a, a)] - P[mp.get(b, b)]
                kv *= 1.0 / (v @ v)
            else:
                (a, b), (c, d_) = f[1:]
                A = P[mp.get(a, a)] - P[mp.get(b, b)]
                B = P[mp.get(c, c)] - P[mp.get(d_, d_)]
                kv *= (A @ B) / ((A @ A) * (B @ B))
        for ws, c in d.items():
            ws2 = [[(l[0], mp.get(l[1], l[1])) if l[0] != 'T' else l for l in w] for w in ws]
            tot += kv * TR.eval_words(c, ws2, cfg.U, cfg.f)
    return tot


if __name__ == '__main__':
    import refP
    S = sources()
    rows = OrderedDict((k, collect(ts)) for k, ts in S.items())
    allterms = [t for ts in S.values() for t in ts]
    total_written = collect(allterms)
    total_sym = collect(allterms, sym=True)
    J = TG.jimwlk()
    J_gen = collect(J)
    J_sym = collect(J, sym=True)
    B = refP.bfkl()
    B_gen = collect(B)
    checks = {}
    dsym = diff(total_sym, J_sym)
    checks['rows_minus_JIMWLK_symbolic_with_UT=U'] = len(dsym)
    dgen = diff(total_written, J_gen)
    checks['rows_as_written_minus_JIMWLK_general_U_monomials'] = len(dgen)
    # numerical checks of the collected forms against the projected targets (general U)
    num = []
    from appF import singlet_value, appF, lobfkl
    for seed in (1, 2, 3):
        cfg = Config(seed, False)
        jn = numeric_value(J_gen, cfg)
        jt = 0.5 * sum(sum(appF(cfg, a, b)) for a, b in (('p2', 'p1'), ('p1', 'p2')))
        bn = numeric_value(B_gen, cfg)
        bt = -3 * 0.5 * sum(lobfkl(cfg, a, b) for a, b in (('p2', 'p1'), ('p1', 'p2')))
        num.append(dict(seed=seed, J_collected=complex(jn).real, appF_corrected=complex(jt).real,
                        B_collected=complex(bn).real, appE_times_minusNc=complex(bt).real))
    checks['numeric'] = num
    sizes = {' / '.join(k): [len(S[k]), sum(len(d) for d in rows[k].values())] for k in S}
    checks['terms_and_monomials_per_source'] = sizes
    checks['JIMWLK_general_monomials'] = sum(len(d) for d in J_gen.values())
    checks['BFKL_general_monomials'] = sum(len(d) for d in B_gen.values())
    checks['kernel_structures'] = dict(rows=len(total_written), JIMWLK=len(J_gen), BFKL=len(B_gen))
    pickle.dump(dict(rows=rows, total_written=total_written, total_sym=total_sym, J_gen=J_gen, J_sym=J_sym,
                     B_gen=B_gen, dgen=dgen), open(os.path.join(HERE, 'singlet.pkl'), 'wb'))
    json.dump(checks, open(os.path.join(HERE, 'singlet_checks.json'), 'w'), indent=1, default=str)
    for k, v in checks.items():
        print(k, ':', v)
