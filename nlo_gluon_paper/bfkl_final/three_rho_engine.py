"""Symmetrized three-rho terms of bfkl.tex (all rows, all c.c.'s) with the reordering engine.
Pointwise test: the sum over all rows, with the three charge positions placed on p1, p2, p3 in all 6 ways (weight 1/6;
the charges are integration variables), vanishes at random points?  As written (K_f in G1 III) and with fix F3."""
import itertools
import json
import numpy as np
import env  # noqa: F401
from core import symmetrize, sub_kern, EXT
import numerics as nm
import rows as R

TS, F = nm.structure_constants(3)
PTS = ['p1', 'p2', 'p3', 'w', "w'", 'z']


class Cfg3:
    def __init__(self, seed, sym=True):
        rng = np.random.default_rng(seed)
        self.P = {l: rng.normal(size=2) for l in PTS}
        self.U = {l: nm.random_adjoint(3, TS, rng, sym) for l in PTS}
        self.r = {l: rng.normal(size=8) for l in PTS}


def kval(kern, P):
    v = 1.0
    for f in kern:
        a, b, c, d = (P[p] for p in f[1:])
        A, B = a - b, c - d
        v *= (A @ B) / ((A @ A) * (B @ B))
    return v


def charges(t):
    ch = []
    for _, p in t['rs']:
        if p not in ch:
            ch.append(p)
    return ch


def value3(terms, cfg):
    tot, scale = 0.0, 0.0
    for t in terms:
        ch = charges(t)
        assert len(ch) == 3, (t['tag'], ch)
        for perm in itertools.permutations(['p1', 'p2', 'p3']):
            mp = dict(zip(ch, perm))
            g = lambda p, mp=mp: mp.get(p, p)
            fl = [('U', i, j, g(p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + \
                 [('rho', i, g(p)) for i, p in t['rs']]
            v = t['c'] * kval(sub_kern(t['kern'], g), cfg.P) * nm.eval_colour(fl, F, cfg.U, cfg.r, 8) / 6
            tot += v
            scale += abs(v)
    return tot, scale


def sym3(F3=True):
    out = {}
    for name, f in R.THREE.items():
        ts = f(F3=F3) if name == 'G1 III' else f()
        out[name] = [s for t in ts for s in symmetrize(t, keep=(3,))]
    return out


if __name__ == '__main__':
    res = {}
    for label, F3 in (('as written', False), ('with F3', True)):
        S = sym3(F3)
        worst, per = [], {}
        for seed in (1, 2, 3, 4, 5):
            cfg = Cfg3(seed)
            tot, sc = 0.0, 0.0
            for name, ts in S.items():
                v, s = value3(ts, cfg)
                tot += v
                sc += s
                per.setdefault(name, []).append(abs(v))
            worst.append(abs(tot) / sc)
        res[label] = dict(rel_residual=worst, n_terms=sum(len(v) for v in S.values()))
        print('%-11s  |sum| / sum|terms| = %s' % (label, ', '.join('%.1e' % x for x in worst)))
    json.dump(res, open('three_rho_engine.json', 'w'), indent=1)
