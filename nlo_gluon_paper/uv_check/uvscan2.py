r"""Transverse UV scan, with principal-value angular averages and whole-row sums.

For every pair (P, Q) of points of an equation (and of a whole row), P = Q + r n(phi), and the integrand is averaged
over phi at fixed r.  If one of the singular curves of surfscan.py (D1, D2, D3, the A_0 segment) crosses the circle,
the phi integral is done as a principal value: Gauss-Legendre between the crossings, and symmetric pairs
f(phi0+t) + f(phi0-t) around each crossing.  A(r) ~ r^s: s <= -2 means int d^2r diverges at P -> Q.
  python3 uvscan2.py [filter]   -> uvscan2[_filter].json"""
import itertools
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uvscan as UV  # noqa: E402
import terms as TM  # noqa: E402
import surfscan as SS  # noqa: E402

RS = [10.0**(-k) for k in (1.5, 2.25, 3, 3.75, 4.5)]
GL_X, GL_W = np.polynomial.legendre.leggauss(8)
SURF_F = [v[0] for v in SS.SURF.values()]
ACTIVE = []


def crossings(pts, P, Q, r):
    """angles where a singular curve crosses the circle |P - Q| = r"""
    phis = np.linspace(0, 2 * np.pi, 721)
    out = []
    for F in ACTIVE:
        def f(ph):
            q = dict(pts)
            q[P] = pts[Q] + r * np.array([np.cos(ph), np.sin(ph)])
            try:
                return F(q)
            except KeyError:
                return 1.0
        vals = [f(ph) for ph in phis]
        for a, b, fa, fb in zip(phis, phis[1:], vals, vals[1:]):
            if fa * fb < 0:
                lo, hi = a, b
                for _ in range(60):
                    m = 0.5 * (lo + hi)
                    if f(lo) * f(m) <= 0:
                        hi = m
                    else:
                        lo = m
                out.append(0.5 * (lo + hi))
    return sorted(out)


def ang_avg(evalf, field, pts, kv, P, Q, r):
    def g(ph):
        q = dict(pts)
        q[P] = pts[Q] + r * np.array([np.cos(ph), np.sin(ph)])
        return evalf(UV.Cfg(field, q, kv))
    cr = crossings(pts, P, Q, r)
    if not cr:
        n = UV.NANG
        return sum(g(2 * np.pi * (k + 0.5) / n) for k in range(n)) / n, 0
    # principal value: half-width h around each crossing, Gauss-Legendre elsewhere
    cr = np.array(cr)
    gaps = np.diff(np.concatenate([cr, [cr[0] + 2 * np.pi]]))
    h = 0.25 * min(gaps.min(), 0.4)
    tot = 0.0
    for c in cr:     # symmetric pairs
        for x, w in zip(GL_X, GL_W):
            t = 0.5 * h * (x + 1)
            tot += 0.5 * h * w * (g(c + t) + g(c - t))
    for c, gap in zip(cr, gaps):   # regular segments
        a, b = c + h, c + gap - h
        for x, w in zip(GL_X, GL_W):
            tot += 0.5 * (b - a) * w * g(0.5 * (a + b) + 0.5 * (b - a) * x)
    return tot / (2 * np.pi), len(cr)


def scan_group(terms, field, pts, kv, pl, pairs=None):
    pts_all = sorted(set().union(*[set(t.pts) for t in terms]), key=UV.NAMES.index)
    res = {}
    for P, Q in itertools.combinations(pts_all, 2):
        if pairs and (P, Q) not in pairs:
            continue
        rel = [t for t in terms if P in t.pts and Q in t.pts]
        if not rel:
            continue
        A, ncr = [], 0
        for r in RS:
            v, n = ang_avg(lambda c: sum(t.value(c, pl) for t in rel), field, pts, kv, P, Q, r)
            A.append(v)
            ncr = max(ncr, n)
        A = np.array(A)
        res['%s-%s' % (P, Q)] = {'s': round(UV.exponent(A), 2), 'r2A': [float(abs(a) * r * r) for a, r in zip(A, RS)],
                                 'crossings': ncr}
    return res


def groups():
    """equations, and rows with more than one equation"""
    eqs = TM.all_equations()
    out = []
    rows = {}
    for name, n, texs in eqs:
        out.append(('%s #%d' % (name, n), texs))
        row = ' / '.join(name.split(' / ')[:2])
        rows.setdefault(row, []).append((n, name, texs))
    for row, lst in rows.items():
        if len(lst) > 1:
            out.append(('%s ROW TOTAL' % row, sum((t for _, _, t in lst), [])))
    return out


def main(only=None):
    out = {}
    fn = os.path.join(HERE, 'uvscan2%s.json' % ('' if not only else '_' + re.sub(r'\W+', '_', only)))
    for key, texs in groups():
        if only and only not in key:
            continue
        try:
            ts = [UV.Term(t) for t in texs]
        except Exception as ex:
            out[key] = {'skipped': str(ex)[:200]}
            print(key, 'skipped', str(ex)[:80], flush=True)
            continue
        longi = any(t.long for t in ts)
        ACTIVE[:] = [f for f, _, pats in SS.SURF.values() if any(SS.has(t, pats) for t in texs)]
        pls = ([0.37, 0.83] + ([2.3] if 'Group III' in key else [])) if longi else [None]
        out[key] = {}
        for seed in (0, 1):
            field = UV.Field(seed)
            pts, kv = UV.base_points(seed)
            for pl in pls:
                r = scan_group(ts, field, pts, kv, pl)
                out[key]['seed%d pl=%s' % (seed, pl)] = r
                bad = {k: v['s'] for k, v in r.items() if v['s'] < -1.7}
                print(key, '| seed', seed, 'pl', pl, '| divergent pairs (s):', bad, flush=True)
        json.dump(out, open(fn, 'w'), indent=1)
    return out


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)
