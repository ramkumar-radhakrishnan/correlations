"""Exact bookkeeping of the two-rho terms by kernel structure and colour network.

Every two-rho term is renamed by roles:  x = charge tied to w, x' = charge tied to w' (both in the
z-free KK factor); if one charge is tied to both, it is called x and the other y.  The kernel key is
then canonical, the colour is simplified (engine.simplify, U^{bc} = U^{cb}), and colour structures are
grouped into networks (equal up to a sign, numerically checked with two random configurations)."""
import numpy as np
import sympy as sp
from collections import OrderedDict, defaultdict
import numerics as nm
from engine import simplify, canonical_rename, tex_colour, Nc
from core import facs, NC

LABS = ['x', "x'", 'y', 'w', "w'", 'z']


def roles(t):
    ch = []
    for _, p in t['rs']:
        if p not in ch:
            ch.append(p)
    lo = [f for f in t['kern'] if f[0] == 'KK' and 'z' not in f[1:]]
    assert len(lo) == 1, (t['tag'], t['kern'])
    a, b, c, d = lo[0][1:]
    tied = {}
    for p, q in ((a, b), (b, a), (c, d), (d, c)):
        if q in ('w', "w'") and p in ch:
            tied.setdefault(q, p)
    cw, cwp = tied['w'], tied["w'"]
    if cw != cwp:
        mp = {cw: 'x', cwp: "x'"}
    else:
        other = [c_ for c_ in ch if c_ != cw]
        mp = {cw: 'x'}
        if other:
            mp[other[0]] = 'y'
    return mp


ORDK = ['x', "x'", 'y', 'w', "w'", 'z']


def vkey(a, b):
    return ((a, b), 1) if ORDK.index(a) < ORDK.index(b) else ((b, a), -1)


def kkey(f):
    if f[0] != 'KK':
        raise ValueError(f)
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


def tex_kkey(key):
    out = []
    for k in key:
        if k[0] == 'inv2':
            out.append(r'\frac{1}{(%s-%s)^2}' % k[1])
        else:
            (a, b), (c, d) = k[1], k[2]
            out.append(r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d))
    return r'\,'.join(out)


def to_exact(c):
    from fractions import Fraction
    c = complex(c)
    re_ = Fraction(c.real).limit_denominator(10000)
    im_ = Fraction(c.imag).limit_denominator(10000)
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
        self.data = OrderedDict()      # kernel key -> list of networks

    def fp(self, fl):
        return [nm.eval_colour(fl, self.f, U, r, NC * NC - 1) for U, r in self.cfgs]

    def add(self, t, group):
        mp = roles(t)
        g = lambda p: mp.get(p, p)
        kern = [(f[0],) + tuple(g(x) for x in f[1:]) for f in t['kern']]
        key, sign = kern_key(kern)
        fl = [('U', i, j, g(p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + [('rho', i, g(p)) for i, p in t['rs']]
        c, fl2, rules = simplify(sp.Integer(1), fl)
        if c == 0 or not fl2:
            return
        fl2 = canonical_rename(fl2)
        coef = to_exact(t['c'] * sign) * c * Nc ** t.get('ncp', 0)
        v = self.fp(fl2)
        if abs(v[0]) < 1e-12:
            return
        nets = self.data.setdefault(key, [])
        for net in nets:
            r = v[0] / net['fp'][0]
            if abs(abs(r) - 1) < 1e-8:
                s = int(round(r.real))
                assert all(abs(a - s * b) < 1e-8 * (1 + abs(a)) for a, b in zip(v, net['fp']))
                net['groups'][group] += coef * s
                if len(fl2) < len(net['col']):
                    net['col'], net['fp'], flip = fl2, v, s
                    for gk in net['groups']:
                        net['groups'][gk] *= flip
                return
        nets.append(dict(col=fl2, fp=v, groups=defaultdict(lambda: sp.Integer(0), {group: coef})))

    def finalize(self):
        for key, nets in self.data.items():
            for net in nets:
                net['groups'] = {g: sp.nsimplify(sp.expand(v)) for g, v in net['groups'].items() if sp.nsimplify(sp.expand(v)) != 0}
        return self
