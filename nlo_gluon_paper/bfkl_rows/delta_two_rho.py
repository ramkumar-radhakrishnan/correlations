r"""What changes in the BFKL list with the exact split (no large-p+ step), in the reordering engine of ../two_rho/core.py
(transcriptions of ../bfkl_final/rows.py):

  Delta = + the log(vee/k+) parts of line 2 of Group III Rows III(1), III(2), V, VI and Group II Row IV (split exactly)
          - (B) of Group III Rows III(1), III(2), V, VI                                     (large-p+ step, withdrawn)
          - H4 of Group II Row IV                                                            (large-p+ step, withdrawn)

1. Every engine form equals the LaTeX term it stands for, pointwise (classical charges, k = 0, log(vee/k+) = 1);
   the pole term of Group II Row IV is compared coefficient by coefficient in 1/eps, 1.
2. Two-rho terms of Delta (genuine two-rho terms and the commutator terms of the three-rho ones, all c.c.'s), key by key
   (kernel monomial at random points; z-free colour integrated over z in dim. reg.).
The earlier list (with (B), H4) had two-rho sum BFKL x LO; the exact list has BFKL x LO + two-rho(Delta).
Output: delta_two_rho.json"""
import importlib.util
import json
import os
import sys
from collections import defaultdict
import numpy as np
import rows as R
import check_rows as CR
import texeval as TE

HERE = R.HERE
BF = os.path.normpath(os.path.join(HERE, '..', 'bfkl_final'))
for sub in ('two_rho', 'three_rho', 'groups23', 'finite', 'finite_v2', 'bfkl_check'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.append(p)
from core import symmetrize, accumulate_N, accumulate_dr, config_pair  # noqa: E402

_spec = importlib.util.spec_from_file_location('bf_rows', os.path.join(BF, 'rows.py'))
BR = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BR)
NC, I, XP, WP, YP = BR.NC, BR.I, BR.XP, BR.WP, BR.YP

UB, UO, EX = CR.UB, CR.UO, CR.EX
NAME = {"x'": 'xp', "w'": 'wp', "y'": 'yp'}
UNIT = TE.fp.G**4 / (16 * np.pi**5) / (2 * np.pi)**3          # N without log(vee/k+), k+ = 1


# ---------------------------------------------------------------- engine forms
def neg(ts):
    return [dict(t, c=-t['c']) for t in ts]


def g2_iv_line2(rename=True):
    """+log(vee/k+) x line 2 of Group II Row IV = minus the 'eps' term of ../bfkl_final/rows.py"""
    return BR.build(2 * NC, [(1, [BR.KK('x', 'y', XP, YP), ('ELL', 'y', YP)])], BR._lo4('y', YP), 'G2 IV line 2',
                    mp={'y': 'w', YP: WP} if rename else None)


def h4(rename=True):
    return BR.build(2 * NC, [(1, [BR.KK('w', 'z', WP, 'z'), BR.KK('x', 'z', XP, 'z')])], BR._lo4('z', 'z'), 'G2 IV H4')


# the line-2 terms are minus the Group III rows of ../bfkl_final/rows.py (same kernels and colour; checked in validate())
LINE2 = {n: (lambda f=f: neg(f())) for n, f in (('G3 III(1)', BR.g3_iii1), ('G3 III(2)', BR.g3_iii2),
                                                 ('G3 V', BR.g3_v), ('G3 VI', BR.g3_vi))}
BTERM = {'G3 III(1)': BR.g3_iii1_P, 'G3 III(2)': lambda: BR.g3_iii2_P(-1), 'G3 V': BR.g3_v_P, 'G3 VI': BR.g3_vi_P}


def plain(fn):
    """the engine row without renaming and without its c.c. (the term as written)"""
    old = BR.build

    def b(coef, kern, colour, tag, mp=None, cc=False):
        return old(coef, kern, colour, tag, None, False)
    BR.build = b
    try:
        return fn()
    finally:
        BR.build = old


# ---------------------------------------------------------------- 1. pointwise validation against the LaTeX
def engine_value(terms, c, eps):
    tot = 0.0
    for t in terms:
        k = 1.0
        for f in t['kern']:
            if f[0] == 'KK':
                a, b, cc_, d = (c.pt[NAME.get(p, p)] for p in f[1:])
                A, B = a - b, cc_ - d
                k *= (A @ B) / ((A @ A) * (B @ B))
            elif f[0] == 'ELL':
                r = c.pt[NAME.get(f[1], f[1])] - c.pt[NAME.get(f[2], f[2])]
                k *= -np.pi * (2 / eps + 2 * np.euler_gamma + np.log((r @ r) * 1.7 / 4))
            else:
                raise ValueError(f)
        fl = [('U', i, j, NAME.get(p, p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + \
             [('r', i, NAME.get(p, p)) for i, p in t['rs']]
        tot += t['c'] * k * TE.Col([(1.0, fl)]).value(c) * NC ** t.get('ncp', 0)
    return tot * UNIT


def text_value(tex, c, eps):
    return TE.term_value(tex, TE.Env(c, KP=1.0, V=np.e, L=1e-9, EPS=eps, MU2=1.7))


def full(B, t):
    return B.prefix + ' (' + t + ') ' + B.suffix


def hidden_B_text():
    """the (B) terms as they were written (../bfkl_final/static.tex, 'The new terms')"""
    import re
    S = open(os.path.join(BF, 'static.tex')).read()
    i = S.index(r'\paragraph{The new terms}')
    blk = S[S.index(r'\begin{align}', i) + len(r'\begin{align}'):S.index(r'\end{align}', i)]
    parts = re.split(r'\\text\{Row ([IV()12]+):\}\\quad&', blk)[1:]
    out = {}
    for k in range(0, len(parts), 2):
        t = parts[k + 1].replace(r'\Lv', r'\log{\left(\frac{\vee}{k^{+}}\right)}').replace(r'\nc', '')
        # (B) has y' -> z in the colour and K(y'-z) -> K(w'-z): the remaining y' of the text is w'
        out['G3 ' + parts[k]] = re.sub(r'[,.]\s*(\\\\)?\s*$', '', t.strip())
    return out


def validate():
    HB = hidden_B_text()
    pairs = [('line 2 ' + n, plain(f), full(UO[n], EX[n][1]['B'])) for n, f in LINE2.items()]
    pairs += [('(B) ' + n, plain(f), HB[n].replace("y'", "w'")) for n, f in BTERM.items()]
    pairs += [('H4 (BFKL.tex G2 IV line 2)', plain(h4), full(UB['G2 IV'], UB['G2 IV'].terms[1]))]
    out = {}
    for lab, eng, tex in pairs:
        rel = []
        for seed in (1, 2, 3):
            c = TE.cfg(seed)
            c.kv = np.zeros(2)
            a, b = engine_value(eng, c, 0.37), text_value(tex, c, 0.37)
            rel.append(float(abs(a - b) / max(abs(a), abs(b))))
        out[lab] = max(rel)
    # Group II Row IV line 2: pole and constant (d_perp = 2 - eps in the text)
    eng = plain(lambda: g2_iv_line2(False))
    tex = full(UO['G2 IV'], EX['G2 IV'][1]['B'])
    t0 = full(UB['G2 IV'], UB['G2 IV'].terms[0])
    es = (0.02, 0.05, 0.1)
    M = np.array([[1 / e, 1.0, e] for e in es])
    rel, rel0 = [], []
    for seed in (1, 2, 3):
        c = TE.cfg(seed)
        c.kv = np.zeros(2)
        a = np.linalg.solve(M, np.array([engine_value(eng, c, e) for e in es]))
        b = np.linalg.solve(M, np.array([text_value(tex, c, e) for e in es]))
        z = np.linalg.solve(M, np.array([text_value(t0, c, e) for e in es]))
        rel.append(float(max(abs(a[0] - b[0]) / abs(b[0]), abs(a[1] - b[1]) / abs(b[1]))))
        rel0.append(float(max(abs(b[0] + z[0]) / abs(b[0]), abs(b[1] + z[1]) / abs(b[1]))))
    out['line 2 G2 IV (1/eps and constant)'] = max(rel)
    out['line 2 G2 IV + your BFKL.tex G2 IV line 1 (1/eps and constant)'] = max(rel0)
    return out


# ---------------------------------------------------------------- 2. two-rho terms, key by key
def has_ell(t):
    return any(k[0] == 'ELL' for k in t['kern'])


def keyed(groups, seed):
    c1, c2 = config_pair(seed, True)
    acc = defaultdict(lambda: defaultdict(complex))
    for lab, ts in groups.items():
        accumulate_N([t for t in ts if not has_ell(t)], c1, c2, acc, label=lab)
        accumulate_dr([t for t in ts if has_ell(t)], c1, acc, label=lab)
    return acc


def sym2(ts):
    return [s for t in ts for s in symmetrize(t, keep=(2,))]


def delta_groups():
    G = {}
    for n, f in LINE2.items():
        G['+line 2 ' + n] = sym2(f())
    G['+line 2 G2 IV'] = sym2(g2_iv_line2())
    for n, f in BTERM.items():
        G['-(B) ' + n] = neg(sym2(f()))
    G['-H4'] = neg(sym2(h4()))
    return G


def nonzero(acc, tol=1e-9):
    bad = []
    for key, d in acc.items():
        if key[0] == 'ell0':
            continue
        tot = sum(d.values())
        sc = sum(abs(v) for v in d.values())
        if abs(tot) > tol * max(1.0, sc):
            bad.append((key, tot, {k: complex(v) for k, v in d.items() if abs(v) > 1e-12}))
    return bad


if __name__ == '__main__':
    res = dict(validation=validate())
    for k, v in res['validation'].items():
        print('%-62s %.1e' % (k, v))
    G = delta_groups()
    parts = {'Delta': G,
             'line 2 terms only': {k: v for k, v in G.items() if k.startswith('+')},
             '(B) and H4 only': {k: v for k, v in G.items() if k.startswith('-')}}
    for name, groups in parts.items():
        r = []
        for seed in (1, 2, 3):
            acc = keyed(groups, seed)
            bad = nonzero(acc)
            r.append(dict(n_keys=sum(1 for k in acc if k[0] != 'ell0'), n_nonzero=len(bad),
                          nonzero=[dict(key=str(k), total=[float(np.real(t)), float(np.imag(t))],
                                        sources={l: [float(np.real(v)), float(np.imag(v))] for l, v in d.items()})
                                   for k, t, d in sorted(bad, key=lambda b: str(b[0]))]))
        res[name] = r
        print('%-20s structures %s, non-zero %s' % (name, [x['n_keys'] for x in r], [x['n_nonzero'] for x in r]))
    for e in res['Delta'][0]['nonzero']:
        print('   %s  %s' % (e['key'][:110], np.round(e['total'], 4)))
    json.dump(res, open(os.path.join(HERE, 'delta_two_rho.json'), 'w'), indent=1)
