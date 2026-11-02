r"""Singular surfaces (not coincident points) of your terms.

Some denominators vanish on whole curves of the transverse plane:
  D1 = (x-z)^2 - (x-w)^2                       (Group I Row III; Group II Row I 2d, 2e)
  D2 = (x'-w')^2 (x-z)^2 - (x-w)^2 (x'-z)^2     (Group II Row I second part)
  D3 = (x'-w')^2 - (x'-z)^2                    (Group II Row I 2d, 2e, 2f)
  S  = w on the segment [x, z]                  (A_0 of Group I Row III, ~ pi/|(x-w)x(x-z)|)
For each equation that contains the structure, a point is moved along a line through the surface (root t0), and the
integrand g is sampled at t0 +- delta.  E = (g+ + g-)/2, O = (g+ - g-)/2, E ~ delta^eE, O ~ delta^eO.
Across the surface:  eE > -1 and eO > -1: integrable;  eO = -1 (eE > -1): integrable only as a principal value;
eE <= -1 or eO <= -2: not integrable (divergent).
  python3 surfscan.py  -> surfscan.json"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uvscan as UV  # noqa: E402
import terms as TM  # noqa: E402

DEL = [10.0**(-k) for k in (2, 2.5, 3, 3.5, 4, 4.5, 5)]


def sq(v):
    return float(v @ v)


SURF = {
    'D1': (lambda P: sq(P['x'] - P['z']) - sq(P['x'] - P['w']), 'z', ['(x-z)^{2}-(x-w)^{2}', '(x-w)^{2}-(x-z)^{2}']),
    'D2': (lambda P: sq(P['xp'] - P['wp']) * sq(P['x'] - P['z']) - sq(P['x'] - P['w']) * sq(P['xp'] - P['z']), 'x',
           ["(x'-w')^{2}\\, (x-z)^{2}-(x-w)^{2}\\, (x'-z)^{2}", "(x'-w')^{2}\\,(x-z)^{2} - (x-w)^{2}\\,(x'-z)^{2}",
            "(x'-w')^{2}\\, (x-z)^{2} - (x-w)^{2}\\, (x'-z)^{2}", "(x'-w')^{2}\\,(x-z)^{2}\\,(x-w)^{2} - (x-w)^{4}",
            "(x'-w')^{4}\\,(x-z)^{2}", "(x'-w')^{2}\\, (x - z)^{2} - (x-w)^{2}\\,(x'-z)^{2}"]),
    'D3': (lambda P: sq(P['xp'] - P['wp']) - sq(P['xp'] - P['z']), 'z', ["(x'-w')^{2}-(x'-z)^{2}", "(x'-z)^{2} - (x'-w')^{2}", "(x'-z)^{2}-(x'-w')^{2}"]),
    'S': (lambda P: float((P['x'] - P['z'])[0] * (P['w'] - P['z'])[1] - (P['x'] - P['z'])[1] * (P['w'] - P['z'])[0]), 'w', [r'\mathcal{A}_{0}']),
}


def has(tex, pats):
    t = tex.replace(' ', '')
    return any(p.replace(' ', '') in t for p in pats)


def base_for(name, seed):
    pts, kv = UV.base_points(seed)
    if name == 'S':           # put w near the segment [x, z]
        pts['w'] = 0.6 * pts['x'] + 0.4 * pts['z'] + 0.05 * np.array([1.0, 0.3])
    return pts, kv


def root(fun, mover, pts, u):
    def f(t):
        q = dict(pts)
        q[mover] = pts[mover] + t * u
        return fun(q)
    ts = np.linspace(-1.5, 1.5, 601)
    vals = [f(t) for t in ts]
    best = None
    for a, b, fa, fb in zip(ts, ts[1:], vals, vals[1:]):
        if fa == 0 or fa * fb < 0:
            lo, hi = a, b
            for _ in range(80):
                m = 0.5 * (lo + hi)
                if f(lo) * f(m) <= 0:
                    hi = m
                else:
                    lo = m
            t0 = 0.5 * (lo + hi)
            if best is None or abs(t0) < abs(best):
                best = t0
    return best


def slope(vals):
    a = np.abs(np.asarray(vals))
    lr = np.log10(DEL[-4:])
    la = np.log10(np.maximum(a[-4:], 1e-300))
    return float(np.polyfit(lr, la, 1)[0])


def probe(evalf, field, pts, kv, mover, u, t0):
    E, O = [], []
    for d in DEL:
        vals = []
        for s in (+1, -1):
            q = dict(pts)
            q[mover] = pts[mover] + (t0 + s * d) * u
            vals.append(evalf(UV.Cfg(field, q, kv)))
        E.append(0.5 * (vals[0] + vals[1]))
        O.append(0.5 * (vals[0] - vals[1]))
    return slope(E), slope(O), [abs(complex(e)) for e in E], [abs(complex(o)) for o in O]


def classify(eE, eO):
    if eE <= -0.85 or eO <= -1.85:
        return 'DIVERGENT'
    if eO <= -0.85:
        return 'principal value'
    return 'integrable'


def main():
    out = {}
    import uvscan2 as U2
    for key, texs in U2.groups():
        try:
            ts = [UV.Term(t) for t in texs]
        except Exception as ex:
            print(key, 'skipped', str(ex)[:60])
            continue
        for sname, (fun, mover, pats) in SURF.items():
            idx = [i for i, t in enumerate(texs) if has(t, pats)]
            if not idx:
                continue
            res = {}
            for seed in (0, 1):
                field = UV.Field(seed)
                pts, kv = base_for(sname, seed)
                rng = np.random.default_rng(seed + 7)
                u = rng.normal(size=2)
                u /= np.linalg.norm(u)
                t0 = root(fun, mover, pts, u)
                if t0 is None:
                    continue
                pl = 0.37 if any(t.long for t in ts) else None
                per = {}
                for i in idx:
                    eE, eO, _, _ = probe(lambda c, t=ts[i]: t.value(c, pl or 0.37), field, pts, kv, mover, u, t0)
                    per[i] = [round(eE, 2), round(eO, 2), classify(eE, eO)]
                eE, eO, Ev, Ov = probe(lambda c: sum(t.value(c, pl or 0.37) for t in ts), field, pts, kv, mover, u, t0)
                res['seed%d' % seed] = {'terms': per, 'sum': [round(eE, 2), round(eO, 2), classify(eE, eO)]}
            out['%s | %s' % (key, sname)] = res
            print(key, sname, res, flush=True)
    json.dump(out, open(os.path.join(HERE, 'surfscan.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
