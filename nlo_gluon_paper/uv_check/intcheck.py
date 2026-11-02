r"""Group I Row V d^3N_2b finite at x -> z with the p+ / q+ integrals done (Gauss-Legendre on [0, k+], Lambda -> 0),
keeping the orientation of each integral and its 1/(2 pi).  The lines have different longitudinal variables and the
11/6 line has none, so they can only cancel after the integrals.   python3 intcheck.py -> intcheck.json"""
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uvscan as UV  # noqa: E402
import uvscan2 as U2  # noqa: E402
import terms as TM  # noqa: E402

XN, XW = np.polynomial.legendre.leggauss(16)
LIM = re.compile(r'\\int([_^])\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}([_^])\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}')


def orientation(tex):
    m = LIM.search(tex)
    if not m or 'frac{d' not in tex[m.end():m.end() + 30]:
        return None
    lo, hi = (m.group(2), m.group(4)) if m.group(1) == '_' else (m.group(4), m.group(2))
    low_is_zero = lo.strip() in ('0', r'\Lambda')
    return +1 if low_is_zero else -1


class IntTerm:
    def __init__(self, tex):
        self.t = UV.Term(tex)
        self.sign = orientation(tex)
        self.pts = self.t.pts

    def value(self, cfg, pl=None):
        if self.sign is None:
            return self.t.raw(cfg)
        tot = 0.0
        for x, w in zip(XN, XW):
            p = 0.5 * (x + 1)
            tot += 0.5 * w * self.t.value(cfg, p)
        return self.sign * tot / (2 * np.pi)


def main():
    for name, n, texs in TM.all_equations():
        if 'Row V / Second part' in name and n == 3:
            break
    ts = [IntTerm(t) for t in texs]
    print('orientations', [t.sign for t in ts], flush=True)
    U2.ACTIVE[:] = []
    out = {}
    for seed in (0, 1):
        field = UV.Field(seed)
        pts, kv = UV.base_points(seed)
        res = {'sum': U2.scan_group(ts, field, pts, kv, None, pairs=[('x', 'z')])['x-z']}
        res['lines'] = [U2.scan_group([t], field, pts, kv, None, pairs=[('x', 'z')])['x-z'] for t in ts]
        out['seed%d' % seed] = res
        print(seed, 'sum s=', res['sum']['s'], 'r^2A=', ['%.2e' % v for v in res['sum']['r2A']],
              '| lines s=', [l['s'] for l in res['lines']], flush=True)
        json.dump(out, open(os.path.join(HERE, 'intcheck.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
