"""Fully pointwise test (no z integration at all): the two integrated genuine terms are replaced by
z-integrands:  Row II:  bracket -> (1/pi) [K(x,w,z) - K(w,w,z)];  Group II Row IV: ELL(w,w') -> K(w,w',z) (+ a K(w,w,z), K(w',w',z))."""
import numpy as np
from collections import defaultdict
from core import Config, accumulate
import target as TG
from compare import report
import rowcmp


def g1_pointwise(ts, alpha=1.0):
    out = []
    for t in ts:
        f = [f for f in t['kern'] if f[0] == 'BG1'][0]
        rest = [g for g in t['kern'] if g[0] != 'BG1']
        x, w = f[1], f[2]
        c = t['c'] / np.pi
        out.append(dict(t, c=c, kern=rest + [('KK', x, 'z', w, 'z')]))
        if alpha:
            out.append(dict(t, c=-alpha * c, kern=rest + [('KK', w, 'z', w, 'z')]))
    return out


def g2iv_pointwise(ts, beta=0.0):
    out = []
    for t in ts:
        f = [f for f in t['kern'] if f[0] == 'ELL'][0]
        rest = [g for g in t['kern'] if g[0] != 'ELL']
        a, b = f[1], f[2]
        out.append(dict(t, kern=rest + [('KK', a, 'z', b, 'z')]))
        if beta:
            out.append(dict(t, c=-beta * t['c'], kern=rest + [('KK', a, 'z', a, 'z')]))
            out.append(dict(t, c=-beta * t['c'], kern=rest + [('KK', b, 'z', b, 'z')]))
    return out


if __name__ == '__main__':
    S = rowcmp.user_sets()
    base = [k for k in S if not k.startswith(('A1', 'A2', 'B1', 'B2', 'D', '-', 'G:G1_II', 'G:G2_IV', 'G:G3_III'))]
    sets = base + ['D', '-RowIII2', '-RowIIaU']
    for beta in (0.0, 0.5):
        cfg = Config(11)
        acc = defaultdict(lambda: defaultdict(complex))
        for k in sets:
            accumulate(S[k], cfg, acc, False, label=k)
        accumulate([dict(t, c=-t['c']) for t in S['G:G3_III']], cfg, acc, False, label='G3_III (sign fixed)')
        accumulate(g1_pointwise(S['G:G1_II']), cfg, acc, False, label='G1_II pointwise')
        accumulate(g2iv_pointwise(S['G:G2_IV'], beta), cfg, acc, False, label='G2_IV pointwise')
        accumulate([dict(t, c=-t['c']) for t in TG.jimwlk()], cfg, acc, False, label='-JIMWLK')
        bad = report(acc, show=False)
        print('beta', beta, ': keys', len(acc), 'nonzero', len(bad))
        for key, tot, sc, d in bad:
            print('   %8.4f %s' % (abs(tot), key), {k: round(v.real, 3) for k, v in d.items() if abs(v) > 1e-9})
