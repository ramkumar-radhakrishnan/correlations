"""Pointwise fingerprints for symmetrized (classical) three-rho terms.
Charges are placed on p1, p2, p3 in all 3! ways (weight 1/6); kernel keys are canonical monomials."""
import env  # noqa
import itertools
import numpy as np
from collections import defaultdict
import numerics as nm
from core import NC, charges, facs, sub_kern

ORD3 = ['p1', 'p2', 'p3', 'p4', 'w', "w'", 'z']


def vkey(a, b):
    return ((a, b), 1) if ORD3.index(a) < ORD3.index(b) else ((b, a), -1)


def kkey(a, b, c, d):
    v1, s1 = vkey(a, b)
    v2, s2 = vkey(c, d)
    if v1 == v2:
        return ('inv2', v1), s1 * s2
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORD3.index(v[0]), ORD3.index(v[1])))), s1 * s2


class Config3:
    def __init__(self, seed, symmetric=True):
        rng = np.random.default_rng(seed)
        self.ts, self.f = nm.structure_constants(NC)
        self.U = {l: nm.random_adjoint(NC, self.ts, rng, symmetric) for l in ORD3}
        self.r = {l: rng.normal(size=NC * NC - 1) for l in ORD3}
        self.ncol = NC * NC - 1

    def colour(self, fl):
        return nm.eval_colour(fl, self.f, self.U, self.r, self.ncol)


def placements3(t):
    ch = charges(t)
    n = len(ch)
    slots = ['p1', 'p2', 'p3', 'p4'][:n]
    maps = [dict(zip(ch, perm)) for perm in itertools.permutations(slots)]
    return maps, 1.0 / len(maps)


def accumulate3(terms, cfg, acc=None, label=None):
    if acc is None:
        acc = defaultdict(lambda: defaultdict(complex))
    for t in terms:
        lab = label or t['tag']
        maps, wgt = placements3(t)
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            kern = sub_kern(t['kern'], g)
            col = cfg.colour(facs(t, mp)) * t['c'] * wgt * NC ** t.get('ncp', 0)
            parts, sign = [], 1
            for f in kern:
                if f[0] == 'PHI':
                    k, s = ('PHI', f[1], f[3]), 1
                else:
                    k, s = kkey(*f[1:])
                parts.append(k)
                sign *= s
            acc[('pt',) + tuple(sorted(parts, key=str))][lab] += sign * col
    return acc


def report3(acc, tol=1e-9, top=5, show=True):
    bad = []
    for k, d in acc.items():
        tot = sum(d.values())
        sc = max([abs(v) for v in d.values()] + [1e-300])
        if abs(tot) > tol * max(1.0, sc):
            bad.append((abs(tot), sc, k, d))
    bad.sort(key=lambda x: -x[0])
    if show:
        print('keys: %d  nonzero residual: %d' % (len(acc), len(bad)))
        for a, sc, k, d in bad[:top]:
            print('  %9.4f %9.4f %s' % (a, sc, k))
            print('            ', {kk: round(v.real, 4) + 1j * round(v.imag, 4) for kk, v in d.items() if abs(v) > 1e-12})
    return bad
