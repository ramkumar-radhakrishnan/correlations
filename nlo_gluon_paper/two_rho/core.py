"""Core of the two-rho comparison.

A term is a dict
    c     : complex coefficient, in units of N = (1/(2pi)^3) g^4/(16 pi^5) (1/k+) log(vee/Lambda)
    Us    : list of (i, j, pos)      U^{ij}(pos)
    fs    : list of (a, b, c)        f^{abc}
    rs    : list of (i, pos)         rho^i(pos) (ordered product; after symmetrization: classical)
    kern  : list of factors
              ('KK', a, b, c, d)     KK(a-b).KK(c-d),  KK(X) = X/X^2
              ('PHI', x, z, w)       -[(x-w)^2+(x-z)^2]/[(x-w)^2-(x-z)^2]^2   (edge kernel of Row III)
              ('ELL', a, b)          int_z K(a,b,z) in dimensional regularization (already integrated)
              ('BG1', x, w)          -2/eps - 2 gamma - log((x-w)^2 mu^2/4)  (Group I Row II, already integrated)
    tag   : provenance label
External positions: 'w' (measured gluon, amplitude), "w'" (conjugate), 'z' (soft gluon).
Every other position is a charge position (integrated)."""
import itertools
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
THREE = os.path.join(HERE, '..', 'three_rho')
for p in (os.path.join(HERE, 'reference'), os.path.join(THREE, 'reference'), THREE):
    p = os.path.normpath(p)
    if p not in sys.path:
        sys.path.insert(0, p)
from weyl import T as WT, weyl as wweyl  # noqa: E402
import numerics as nm  # noqa: E402

EXT = ('w', "w'", 'z')
NC = 3


def term(c, Us, fs, rs, kern, tag='', ncp=0):
    """ncp: explicit power of Nc multiplying the coefficient (kept symbolic for the tables)"""
    return dict(c=complex(c), Us=list(Us), fs=list(fs), rs=list(rs), kern=list(kern), tag=tag, ncp=ncp)


def chase(mp, p):
    n = 0
    while p in mp and mp[p] != p and n < 20:
        p = mp[p]
        n += 1
    return p


def sub_kern(kern, g):
    out = []
    for f in kern:
        out.append((f[0],) + tuple(g(x) for x in f[1:]))
    return out


def symmetrize(t, keep=(2,)):
    """Weyl decomposition with [rho^a(p), rho^b(q)] = i f^{abc} rho^c(p) delta(p-q);
    returns the terms with len(rs) in keep (rho's are then classical)."""
    out = []
    for w in wweyl(WT(t['c'], t['Us'], t['fs'], t['rs']), 2):
        if len(w.rs) not in keep:
            continue
        g = lambda p, mp=w.mp: chase(mp, p)
        out.append(term(w.c, [(i, j, g(p)) for i, j, p in w.Us], w.fs, [(i, g(p)) for i, p in w.rs],
                        sub_kern(t['kern'], g), t['tag'], t.get('ncp', 0)))
    return out


def rename(t, mp):
    g = lambda p: mp.get(p, p)
    return term(t['c'], [(i, j, g(p)) for i, j, p in t['Us']], t['fs'], [(i, g(p)) for i, p in t['rs']],
                sub_kern(t['kern'], g), t['tag'], t.get('ncp', 0))


def conj_swap(t, P1="w'", P0='w'):
    """complex conjugate of an ordered term followed by P1 <-> P0 (restores the phase)."""
    s = term(np.conj(t['c']), t['Us'], t['fs'], t['rs'][::-1], t['kern'], t['tag'] + '*', t.get('ncp', 0))
    return rename(s, {P1: P0, P0: P1})


def facs(t, posmap=None):
    g = (lambda p: posmap.get(p, p)) if posmap else (lambda p: p)
    out = [('U', i, j, g(p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']]
    out += [('rho', i, g(p)) for i, p in t['rs']]
    return out


def charges(t):
    ch = []
    for _, p in t['rs']:
        if p not in ch:
            ch.append(p)
    for i, j, p in t['Us']:
        if p not in EXT and p not in ch:
            raise ValueError('Wilson line at a non-charge, non-external point: %s in %s' % (p, t['tag']))
    for f in t['kern']:
        for x in f[1:]:
            if x not in EXT and x not in ch:
                raise ValueError('kernel point %s is neither external nor a charge (%s)' % (x, t['tag']))
    return ch


# ---------------------------------------------------------------- canonical kernel keys
ORD = ['p1', 'p2', 'w', "w'", 'z']


def vkey(a, b):
    return ((a, b), 1) if ORD.index(a) < ORD.index(b) else ((b, a), -1)


def kkey(a, b, c, d):
    v1, s1 = vkey(a, b)
    v2, s2 = vkey(c, d)
    if v1 == v2:
        return ('inv2', v1), s1 * s2
    return ('KK',) + tuple(sorted([v1, v2], key=lambda v: (ORD.index(v[0]), ORD.index(v[1])))), s1 * s2


def factor_key(f):
    if f[0] == 'KK':
        return kkey(*f[1:])
    if f[0] == 'PHI':
        return ('PHI', f[1], f[3]), 1
    if f[0] in ('ELL', 'BG1'):
        a, b = sorted(f[1:], key=ORD.index)
        return (f[0], a, b), 1
    raise ValueError(f)


def soft_pair(kern):
    """if z appears in exactly one KK factor of the form KK(a-z).KK(b-z) (any orientation),
    return (index, a, b, sign) with factor = sign * K(a,b,z); else None."""
    hits = [k for k, f in enumerate(kern) if 'z' in f[1:]]
    if len(hits) != 1 or kern[hits[0]][0] != 'KK':
        return None
    k = hits[0]
    a, b, c, d = kern[k][1:]
    sign = 1
    if a == 'z':
        a, b, sign = b, a, -sign
    if c == 'z':
        c, d, sign = d, c, -sign
    if b != 'z' or d != 'z' or 'z' in (a, c):
        return None
    return k, a, c, sign


class Config:
    """random points, symmetric (or general) adjoint Wilson lines, random classical charges"""
    def __init__(self, seed, symmetric=True):
        rng = np.random.default_rng(seed)
        self.ts, self.f = nm.structure_constants(NC)
        self.P = {l: rng.normal(size=2) for l in ORD}
        self.U = {l: nm.random_adjoint(NC, self.ts, rng, symmetric) for l in ORD}
        self.r = {l: rng.normal(size=NC * NC - 1) for l in ORD}
        self.ncol = NC * NC - 1

    def colour(self, fl):
        return nm.eval_colour(fl, self.f, self.U, self.r, self.ncol)


def placements(t):
    ch = charges(t)
    if len(ch) == 1:
        return [{ch[0]: 'p1'}], 1.0
    if len(ch) == 2:
        return [{ch[0]: 'p1', ch[1]: 'p2'}, {ch[0]: 'p2', ch[1]: 'p1'}], 0.5
    raise ValueError('expected 1 or 2 charges, got %s (%s)' % (ch, t['tag']))


def accumulate(terms, cfg, acc=None, integrate=False, label=None):
    """acc[key][label] += value.  key = canonical kernel monomial (+ sector info).
    integrate=True: terms whose colour does not depend on z and whose kernel contains z only in
    one factor K(a,b,z) are integrated over z (dim. reg.): int_z K(a,b,z) = pi[P - log (a-b)^2], 0 if a=b."""
    from collections import defaultdict
    if acc is None:
        acc = defaultdict(lambda: defaultdict(complex))
    for t in terms:
        lab = label or t['tag']
        maps, wgt = placements(t)
        zcol = any(p == 'z' for _, _, p in t['Us'])
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            kern = sub_kern(t['kern'], g)
            col = cfg.colour(facs(t, mp)) * t['c'] * wgt * NC ** t.get('ncp', 0)
            sp_ = soft_pair(kern) if (integrate and not zcol) else None
            if sp_ is None:
                parts, sign = [], 1
                for f in kern:
                    k, s = factor_key(f)
                    parts.append(k)
                    sign *= s
                key = ('pt',) + tuple(sorted(parts, key=str))
                acc[key][lab] += sign * col
            else:
                k0, a, b, s0 = sp_
                rest = [f for k, f in enumerate(kern) if k != k0]
                parts, sign = [], s0
                for f in rest:
                    k, s = factor_key(f)
                    parts.append(k)
                    sign *= s
                rk = tuple(sorted(parts, key=str))
                if a == b:
                    acc[('ell0',) + rk + (a,)][lab] += sign * col
                else:
                    aa, bb = sorted([a, b], key=ORD.index)
                    acc[('P',) + rk][lab] += np.pi * sign * col
                    acc[('LOG', aa, bb) + rk][lab] += -np.pi * sign * col
    return acc


def expand_dimreg(terms):
    """('ELL',a,b) -> pi P - pi log; ('BG1',x,w) -> BG1 pole const - log(x-w)^2, recorded as keys by accumulate_dr"""
    return terms


def accumulate_dr(terms, cfg, acc, label=None):
    """already-integrated genuine terms: kernel = [LO-like KK factor] x ELL or BG1"""
    for t in terms:
        lab = label or t['tag']
        maps, wgt = placements(t)
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            kern = sub_kern(t['kern'], g)
            col = cfg.colour(facs(t, mp)) * t['c'] * wgt * NC ** t.get('ncp', 0)
            special = [f for f in kern if f[0] in ('ELL', 'BG1')]
            rest = [f for f in kern if f[0] not in ('ELL', 'BG1')]
            assert len(special) == 1
            parts, sign = [], 1
            for f in rest:
                k, s = factor_key(f)
                parts.append(k)
                sign *= s
            rk = tuple(sorted(parts, key=str))
            sp_ = special[0]
            aa, bb = sorted(sp_[1:], key=ORD.index)
            if sp_[0] == 'ELL':
                acc[('P',) + rk][lab] += np.pi * sign * col
                acc[('LOG', aa, bb) + rk][lab] += -np.pi * sign * col
            else:   # BG1 = -2/eps - 2 gamma - log(mu^2/4) - log (x-w)^2
                acc[('PG1',) + rk][lab] += sign * col
                acc[('LOG', aa, bb) + rk][lab] += -sign * col
    return acc


def kernel_value(kern, P):
    v = 1.0
    for f in kern:
        if f[0] == 'KK':
            a, b, c, d = f[1:]
            A, B = P[a] - P[b], P[c] - P[d]
            v *= (A @ B) / ((A @ A) * (B @ B))
        elif f[0] == 'PHI':
            x, z, w = (P[l] for l in f[1:])
            a, b = (x - w) @ (x - w), (x - z) @ (x - z)
            v *= -(a + b) / (a - b) ** 2
        else:
            raise ValueError(f)
    return v


def pointwise_total(terms, cfg):
    """numerical value of the sum at the configuration's points (no keys, no integration)."""
    tot, scale = 0.0, 0.0
    for t in terms:
        maps, wgt = placements(t)
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            P = {l: cfg.P[l] for l in ORD}
            v = kernel_value(sub_kern(t['kern'], g), P) * cfg.colour(facs(t, mp)) * t['c'] * wgt * NC ** t.get('ncp', 0)
            tot += v
            scale += abs(v)
    return tot, scale


def zfree(t, cfg, cfg2):
    """numerically: does the colour factor of t depend on U(z)?  cfg2 differs from cfg only in U(z)."""
    maps, _ = placements(t)
    mp = maps[0]
    a = cfg.colour(facs(t, mp))
    b = cfg2.colour(facs(t, mp))
    return abs(a - b) < 1e-10 * (1 + abs(a))


def accumulate_N(terms, cfg, cfg2, acc, label=None):
    """like accumulate(integrate=True) but the z-dependence of the colour is detected numerically."""
    from collections import defaultdict
    for t in terms:
        lab = label or t['tag']
        integ = zfree(t, cfg, cfg2)
        maps, wgt = placements(t)
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            kern = sub_kern(t['kern'], g)
            col = cfg.colour(facs(t, mp)) * t['c'] * wgt * NC ** t.get('ncp', 0)
            sp_ = soft_pair(kern) if integ else None
            if sp_ is None:
                parts, sign = [], 1
                for f in kern:
                    k, s = factor_key(f)
                    parts.append(k)
                    sign *= s
                acc[('pt',) + tuple(sorted(parts, key=str))][lab] += sign * col
            else:
                k0, a, b, s0 = sp_
                rest = [f for k, f in enumerate(kern) if k != k0]
                parts, sign = [], s0
                for f in rest:
                    k, s = factor_key(f)
                    parts.append(k)
                    sign *= s
                rk = tuple(sorted(parts, key=str))
                if a == b:
                    acc[('ell0', a) + rk][lab] += sign * col
                else:
                    aa, bb = sorted([a, b], key=ORD.index)
                    acc[('P',) + rk][lab] += np.pi * sign * col
                    acc[('LOG', aa, bb) + rk][lab] += -np.pi * sign * col
    return acc


def config_pair(seed, symmetric=True):
    c1 = Config(seed, symmetric)
    c2 = Config(seed, symmetric)
    rng = np.random.default_rng(seed + 1000)
    c2.U['z'] = nm.random_adjoint(NC, c2.ts, rng, symmetric)
    return c1, c2
