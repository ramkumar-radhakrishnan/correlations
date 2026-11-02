r"""Transverse UV scan of every term of your Other_terms.tex.

Wilson lines and charges are smooth functions of the transverse position:
  U(p) = adjoint of V(p) = W(p) diag(1,-1,-1) W(p)^dagger,  W(p) = exp(i theta_a(p) t^a)  (SU(3), U = U^T),
  rho^a(p), theta_a(p): sums of sines with O(1) wave numbers.
For every pair of points (P, Q) of a term, P = Q + r n(phi): the integrand is averaged over 16 equally spaced angles
(exact for angular harmonics < 16) and the average A(r) is followed for r = 1e-1 ... 1e-5.  In two dimensions
int d^2r A(r) converges at r -> 0 iff r^2 A(r) -> 0.  We report the exponent s in A ~ r^s (s <= -2: divergent).
Longitudinal integrals are kept as parameters (p+ = 0.37 k+, 0.83 k+ and, for Group III, 2.3 k+); a plus prescription
[f]_+ is evaluated with its subtraction at p+ = 0.
  python3 uvscan.py            -> uvscan.json"""
import itertools
import json
import os
import re
import sys

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'three_rho'))
import terms as TM  # noqa: E402
import texeval as TE  # noqa: E402
import texfix as TF  # noqa: E402
import numerics as nm  # noqa: E402

TS, F = nm.structure_constants(3)
TE.F = F
KP, VEE = 1.0, 3.0
NAMES = ['x', 'y', 'z', 'w', 'xp', 'wp', 'yp']


class Field:
    def __init__(self, seed):
        rng = np.random.default_rng(seed)
        self.q = rng.normal(size=(8, 3, 2)) * 0.9
        self.ph = rng.uniform(0, 2 * np.pi, size=(8, 3))
        self.c = rng.normal(size=(8, 3)) * 0.8
        self.qr = rng.normal(size=(8, 3, 2)) * 0.9
        self.phr = rng.uniform(0, 2 * np.pi, size=(8, 3))
        self.cr = rng.normal(size=(8, 3))

    def U(self, p):
        key = (float(p[0]), float(p[1]))
        cache = self.__dict__.setdefault('_uc', {})
        if key in cache:
            return cache[key]
        cache[key] = self._U(p)
        return cache[key]

    def _U(self, p):
        th = (self.c * np.sin(self.q @ p + self.ph)).sum(1)
        W = expm(1j * sum(th[a] * TS[a] for a in range(8)))
        V = W @ np.diag([1, -1, -1]) @ W.conj().T
        return nm.adjoint(V, TS)

    def rho(self, p):
        return (self.cr * np.sin(self.qr @ p + self.phr)).sum(1) + 0.3


class Cfg:
    def __init__(self, field, pts, kv):
        self.pt = {k: np.asarray(v, float) for k, v in pts.items()}
        self.U = {k: field.U(v) for k, v in self.pt.items()}
        self.rho = {k: field.rho(v) for k, v in self.pt.items()}
        self.kv = kv


class Term:
    def __init__(self, tex):
        self.tex = tex
        self.code = TM.code_of(tex)
        self.cc = compile(self.code, '<term>', 'eval')
        self.pts = sorted(set(re.findall(r'"(xp|yp|wp|x|y|z|w)"', self.code)) & set(NAMES), key=NAMES.index)
        self.plus = [p for p in ('PLUSKP', 'PLUSKQ', 'PLUS1Q', 'PLUSK') if p in self.code]
        self.long = 'PP' in self.code or 'QQ' in self.code or bool(self.plus)

    def raw(self, cfg, PP=0.0, QQ=0.0, plus=0.0):
        env = TE.Env(cfg, KP=KP, V=VEE, L=0.0, EPS=0.0)

        def names_of(ix):
            d = env.names(ix)
            d.update(PP=PP, QQ=QQ, PLUSK=plus, PLUSKQ=plus, PLUS1Q=plus, PLUSKP=plus)
            return d
        return TF.eval_code(self.code, names_of, post=lambda v: TE.scalar_or_col(v, cfg))

    def value(self, cfg, pl=0.37):
        """integrand at the longitudinal parameter pl (in units of k+), plus prescriptions resolved"""
        if not self.long:
            return self.raw(cfg)
        p = pl * KP
        if not self.plus:
            return self.raw(cfg, PP=p, QQ=p)
        kind = self.plus[0]
        fac = {'PLUSK': KP / p, 'PLUSKQ': KP / p, 'PLUS1Q': 1.0 / p, 'PLUSKP': KP / (p * (p + KP))}[kind]
        A = self.raw(cfg, PP=p, QQ=p, plus=0.0)
        B = self.raw(cfg, PP=p, QQ=p, plus=1.0) - A
        A0 = self.raw(cfg, PP=0.0, QQ=0.0, plus=0.0)
        B0 = self.raw(cfg, PP=0.0, QQ=0.0, plus=1.0) - A0
        return A + fac * (B - B0)


def base_points(seed):
    rng = np.random.default_rng(1000 + seed)
    return {n: rng.normal(size=2) * 0.9 for n in NAMES}, rng.normal(size=2)


RS = [10.0**(-k) for k in (1.5, 2.25, 3, 3.75, 4.5)]
NANG = 16


def pair_scan(evalf, field, pts, kv, P, Q):
    """A(r) for P -> Q (angle average)"""
    out = []
    for r in RS:
        acc = 0.0
        for k in range(NANG):
            ph = 2 * np.pi * (k + 0.5) / NANG
            q = dict(pts)
            q[P] = pts[Q] + r * np.array([np.cos(ph), np.sin(ph)])
            acc += evalf(Cfg(field, q, kv))
        out.append(acc / NANG)
    return np.array(out)


def exponent(A):
    """slope of log|A| vs log r over the three smallest r (A ~ r^s)"""
    a = np.abs(np.real(A)) + np.abs(np.imag(A))
    if np.all(a[-3:] < 1e-300):
        return 99.0
    lr = np.log10(RS[-3:])
    la = np.log10(np.maximum(a[-3:], 1e-300))
    return float(np.polyfit(lr, la, 1)[0])


def scan_eq(terms, field, pts, kv, pl):
    pts_all = sorted(set().union(*[set(t.pts) for t in terms]), key=NAMES.index)
    res = {}
    for P, Q in itertools.combinations(pts_all, 2):
        per = []
        for t in terms:
            if P in t.pts and Q in t.pts:
                per.append(pair_scan(lambda c, t=t: t.value(c, pl), field, pts, kv, P, Q))
            else:
                per.append(None)
        tot = sum(a for a in per if a is not None) if any(a is not None for a in per) else None
        res['%s-%s' % (P, Q)] = {'terms': [None if a is None else exponent(a) for a in per],
                                 'sum': None if tot is None else exponent(tot),
                                 'sum_r2A': None if tot is None else [float(abs(v) * r * r) for v, r in zip(tot, RS)]}
    return res


def main(only=None):
    eqs = TM.all_equations()
    out = {}
    for name, n, texs in eqs:
        key = '%s #%d' % (name, n)
        if only and only not in key:
            continue
        try:
            ts = [Term(t) for t in texs]
        except Exception as ex:
            out[key] = {'skipped': str(ex)[:200]}
            print(key, 'skipped', str(ex)[:80], flush=True)
            continue
        pls = [0.37, 0.83] + ([2.3] if 'Group III' in name else []) if any(t.long for t in ts) else [None]
        out[key] = {}
        for seed in (0, 1):
            field = Field(seed)
            pts, kv = base_points(seed)
            for pl in pls:
                r = scan_eq(ts, field, pts, kv, pl)
                out[key]['seed%d pl=%s' % (seed, pl)] = r
                bad = {k: v for k, v in r.items() if v['sum'] is not None and v['sum'] < -1.7}
                badt = {k: [i for i, e in enumerate(v['terms']) if e is not None and e < -1.7] for k, v in r.items()}
                badt = {k: v for k, v in badt.items() if v}
                print(key, 'seed', seed, 'pl', pl, '| SUM divergent:', {k: round(v['sum'], 2) for k, v in bad.items()},
                      '| terms divergent:', badt, flush=True)
        json.dump(out, open(os.path.join(HERE, 'uvscan%s.json' % ('' if not only else '_' + re.sub(r'\W+', '_', only))), 'w'), indent=1)
    return out



# ---- speed: memoize the colour contractions per configuration
_orig_value = TE.Col.value


def _cached_value(self, cfg):
    cache = cfg.__dict__.setdefault('_cc', {})
    tot = 0.0
    for c, facs in self.m:
        if c == 0:
            continue
        key = tuple(facs)
        if key not in cache:
            cache[key] = _contract(facs, cfg)
        tot += c * cache[key]
    return tot


def _contract(facs, cfg):
    labels, ops, subs = {}, [], []
    for fa in facs:
        if fa[0] == 'U':
            ops.append(cfg.U[fa[3]]); idx = fa[1:3]
        elif fa[0] == 'r':
            ops.append(cfg.rho[fa[2]]); idx = fa[1:2]
        else:
            ops.append(F); idx = fa[1:4]
        sub = ''
        for i in idx:
            if i not in labels:
                labels[i] = 'abcdefghijklmnopqrstuvwxyz'[len(labels)]
            sub += labels[i]
        subs.append(sub)
    return np.einsum(','.join(subs) + '->', *ops, optimize='greedy')


TE.Col.value = _cached_value


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
