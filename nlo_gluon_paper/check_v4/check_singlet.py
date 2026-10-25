"""The colour-singlet forms of Secs. 3, 4 of the new Evolution_terms.tex.
 (a) against JIMWLK x LO and BFKL x LO in singlet form (Evolution_singlet.pdf Eqs. (3.1), (3.2), general U);
 (b) against the singlet projection of your own collected forms of Secs. 3, 4 (rho^a rho^b -> delta^{ab});
 (c) the same for U = U^T.
x -> p1, y (or x') -> p2, symmetrized p1 <-> p2.  Output: check_singlet.json"""
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SING = os.path.normpath(os.path.join(HERE, '..', 'singlet'))
sys.path.insert(0, SING)
import env  # noqa: E402,F401  (singlet/env.py)
sys.path.append(HERE)
import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
import trace as TR  # noqa: E402  (singlet/trace.py)
from core import Config, kernel_value, placements, sub_kern  # noqa: E402
from render_singlet import numeric_value  # noqa: E402
from parse_v4 import equations, parse_traces  # noqa: E402

NC = 3
D = pickle.load(open(os.path.join(SING, 'singlet.pkl'), 'rb'))


def user_value(mons, cfg, mp):
    tot = 0.0
    P = cfg.P
    for c, ncp, kern, traces in mons:
        g = lambda p: mp.get(p, p)
        kv = kernel_value(sub_kern(kern, g), P)
        # product over traces of the sum over alternatives
        alts = [[]]
        for tr in traces:
            alts = [a + [(cc, w)] for a in alts for cc, w in tr]
        for a in alts:
            co = c * NC ** ncp
            words = []
            for cc, w in a:
                co *= cc
                words.append([(l[0], g(l[1])) if l[0] != 'T' else l for l in w])
            tot += kv * TR.eval_words(sp.Integer(1) * co, words, cfg.U, cfg.f).real
    return tot


def user_sym(mons, cfg):
    return 0.5 * (user_value(mons, cfg, {'x': 'p1', 'y': 'p2', "x'": 'p2'}) +
                  user_value(mons, cfg, {'x': 'p2', 'y': 'p1', "x'": 'p1'}))


def projected(terms, cfg):
    """singlet projection of two-rho terms (rho rho -> delta), symmetrized over the charge placements"""
    tot = 0.0
    for t in terms:
        c, ws = TR.singlet(t)
        if c == 0:
            continue
        maps, wgt = placements(t)
        kk = t['kern']
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            kv = kernel_value(sub_kern(kk, g), cfg.P)
            words = [[(l[0], g(l[1])) if l[0] != 'T' else l for l in w] for w in ws]
            tot += wgt * kv * TR.eval_words(c, words, cfg.U, cfg.f).real
    return tot


if __name__ == '__main__':
    E = equations()
    Js = parse_traces([e for e in E if e['kind'] == 'singlet' and e['sec'].startswith('JIMWLK')][0]['text'])
    Bs = parse_traces([e for e in E if e['kind'] == 'singlet' and e['sec'].startswith('BFKL')][0]['text'])
    Jc = [e for e in E if e['sec'].startswith('JIMWLK') and e['kind'] == 'terms'][2]['terms']
    Bc = [e for e in E if e['sec'].startswith('BFKL') and e['kind'] == 'terms'][1]['terms']
    out = {}
    rows = []
    for sym in (False, True):
        for seed in (1, 2, 3):
            cfg = Config(seed, sym)
            r = dict(sym=sym, seed=seed,
                     J_user_singlet=user_sym(Js, cfg), J_note=numeric_value(D['J_gen'], cfg), J_proj_user_sec3=projected(Jc, cfg),
                     B_user_singlet=user_sym(Bs, cfg), B_note=numeric_value(D['B_gen'], cfg), B_proj_user_sec4=projected(Bc, cfg))
            rows.append({k: (float(np.real(v)) if not isinstance(v, bool) else v) for k, v in r.items()})
            print(rows[-1])
    out['values'] = rows
    rel = lambda a, b: abs(a - b) / max(abs(a), abs(b))
    for sym in (False, True):
        R = [r for r in rows if r['sym'] == sym]
        tag = 'U=U^T' if sym else 'general U'
        out[tag] = dict(J_singlet_vs_note=max(rel(r['J_user_singlet'], r['J_note']) for r in R),
                        J_singlet_vs_proj_sec3=max(rel(r['J_user_singlet'], r['J_proj_user_sec3']) for r in R),
                        B_singlet_vs_note=max(rel(r['B_user_singlet'], r['B_note']) for r in R),
                        B_singlet_vs_proj_sec4=max(rel(r['B_user_singlet'], r['B_proj_user_sec4']) for r in R))
        print(tag, out[tag])
    out['n_monomials'] = dict(J=len(Js), B=len(Bs))
    json.dump(out, open(os.path.join(HERE, 'check_singlet.json'), 'w'), indent=1)
