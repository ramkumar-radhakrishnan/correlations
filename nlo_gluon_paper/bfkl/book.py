"""Exact bookkeeping by kernel structure and colour network, for two- or three-rho (classical) terms
with arbitrary products of KK factors.

Charge names: x, x', y.  For every term all injective namings of its charges are tried; the naming
with the best score (charge tied to w called x, charge tied to w' called x') and the smallest key is
used.  If several namings tie, the term is averaged over them (exact, weight 1/m), so that the sums
per network are independent of the choice."""
import itertools
import numpy as np
import sympy as sp
from collections import OrderedDict, defaultdict
import env  # noqa
import numerics as nm
from engine import simplify, canonical_rename, Nc
from core import NC

LABS = ['x', "x'", 'y', 'w', "w'", 'z']
ORDK = LABS


def vkey(a, b):
    return ((a, b), 1) if ORDK.index(a) < ORDK.index(b) else ((b, a), -1)


def kkey(f):
    if f[0] == 'PHI':
        return ('PHI', f[1], f[3]), 1
    v1, s1 = vkey(*f[1:3])
    v2, s2 = vkey(*f[3:5])
    if v1 == v2:
        return ('inv2', v1), s1 * s2
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORDK.index(v[0]), ORDK.index(v[1])))), s1 * s2


def kern_key(kern):
    parts, sign = [], 1
    for f in kern:
        k, s = kkey(f)
        parts.append(k)
        sign *= s
    return tuple(sorted(parts, key=str)), sign


def score(key):
    s = 0
    for k in key:
        pairs = []
        if k[0] == 'KK':
            pairs = [k[1], k[2]]
        elif k[0] == 'inv2':
            pairs = [k[1]]
        for a, b in pairs:
            if (a, b) == ('x', 'w'):
                s += 4
            if (a, b) == ("x'", "w'"):
                s += 4
            if 'y' in (a, b):
                s -= 1
            if (a, b) == ('x', "x'"):
                s += 1
    return s


def charges(t):
    ch = []
    for _, p in t['rs']:
        if p not in ch:
            ch.append(p)
    return ch


def namings(t):
    ch = charges(t)
    best = None
    out = []
    for labs in itertools.permutations(['x', "x'", 'y'], len(ch)):
        mp = dict(zip(ch, labs))
        g = lambda p, mp=mp: mp.get(p, p)
        kern = [(f[0],) + tuple(g(x) for x in f[1:]) for f in t['kern']]
        key, sign = kern_key(kern)
        rank = (-score(key), str(key))
        if best is None or rank < best:
            best, out = rank, [(mp, key, sign)]
        elif rank == best:
            out.append((mp, key, sign))
    return out


def tex_kkey(key):
    out = []
    for k in key:
        if k[0] == 'inv2':
            out.append(r'\frac{1}{(%s-%s)^2}' % k[1])
        elif k[0] == 'PHI':
            out.append(r'\Phi_{%s}(%s)' % (k[2], k[1]))
        else:
            (a, b), (c, d) = k[1], k[2]
            out.append(r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d))
    return r'\,'.join(out)


def to_exact(c):
    from fractions import Fraction
    c = complex(c)
    re_ = Fraction(c.real).limit_denominator(100000)
    im_ = Fraction(c.imag).limit_denominator(100000)
    assert abs(float(re_) - c.real) < 1e-9 and abs(float(im_) - c.imag) < 1e-9, c
    return sp.Rational(re_.numerator, re_.denominator) + sp.I * sp.Rational(im_.numerator, im_.denominator)


class Book:
    def __init__(self, seeds=(5, 9)):
        self.cfgs = []
        ts, f = nm.structure_constants(NC)
        self.f = f
        for s in seeds:
            rng = np.random.default_rng(s)
            U = {l: nm.random_adjoint(NC, ts, rng, True) for l in LABS}
            r = {l: rng.normal(size=NC * NC - 1) for l in LABS}
            self.cfgs.append((U, r))
        self.data = OrderedDict()
        self.cache = {}

    def fp(self, fl):
        return [nm.eval_colour(fl, self.f, U, r, NC * NC - 1) for U, r in self.cfgs]

    def add(self, t, group):
        nm_ = namings(t)
        w = sp.Rational(1, len(nm_))
        for mp, key, sign in nm_:
            g = lambda p, mp=mp: mp.get(p, p)
            fl = [('U', i, j, g(p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + [('rho', i, g(p)) for i, p in t['rs']]
            ck = repr(fl)
            if ck in self.cache:
                c, fl2 = self.cache[ck]
            else:
                c, fl2, rules = simplify(sp.Integer(1), fl)
                if c != 0 and fl2:
                    fl2 = canonical_rename(fl2)
                self.cache[ck] = (c, fl2)
            if c == 0 or not fl2:
                continue
            coef = to_exact(t['c'] * sign) * c * w * Nc ** t.get('ncp', 0)
            v = self.fp(fl2)
            if abs(v[0]) < 1e-12 and abs(v[1]) < 1e-12:
                continue
            nets = self.data.setdefault(key, [])
            for net in nets:
                j = 0 if abs(net['fp'][0]) > 1e-12 else 1
                r = v[j] / net['fp'][j]
                if abs(abs(r) - 1) < 1e-8:
                    s = int(round(r.real))
                    if all(abs(a - s * b) < 1e-8 * (1 + abs(a)) for a, b in zip(v, net['fp'])):
                        net['groups'][group] += coef * s
                        return_flag = True
                        break
            else:
                nets.append(dict(col=fl2, fp=v, groups=defaultdict(lambda: sp.Integer(0), {group: coef})))

    def finalize(self):
        for key, nets in self.data.items():
            for net in nets:
                net['groups'] = {g: sp.nsimplify(sp.expand(v)) for g, v in net['groups'].items()
                                 if sp.nsimplify(sp.expand(v)) != 0}
        return self

    def dump(self):
        return {k: [dict(col=n['col'], groups=dict(n['groups'])) for n in v] for k, v in self.data.items()}
