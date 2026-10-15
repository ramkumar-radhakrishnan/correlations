"""Compare the new list (v3) with the corrected list of ../three_rho_v2 (which cancels), group by group,
on the symmetrised three-rho part (pointwise, random U symmetric)."""
import sys, os
from collections import defaultdict
import numpy as np
import env  # noqa
from core import symmetrize
import acc3 as A3
import parse3 as P3
V2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'three_rho_v2')
sys.path.insert(0, V2)
import parse as P2
import fixes as FX

GRP = lambda b: {'4rho.a': '4rho', '4rho.b': '4rho', 'I.II.1': 'I.II', 'I.II.2': 'I.II', 'D': 'I.III', 'I.V.1': 'I.V', 'I.V.2': 'I.V'}.get(b, b)


def acc_of(entries, as_term, seed, sym=True):
    cfg = A3.Config3(seed, sym)
    out = defaultdict(lambda: defaultdict(complex))
    for e in entries:
        acc = defaultdict(lambda: defaultdict(complex))
        A3.accumulate3(symmetrize(as_term(e), (3,)), cfg, acc, label='x')
        for k, d in acc.items():
            out[GRP(e['block'])][(k, e['cls'])] += d['x']
    return out


if __name__ == '__main__':
    a3 = acc_of(P3.parse(), P3.as_term, 5)
    a2 = acc_of(FX.corrected(), P2.as_term, 5)
    for g in ['4rho', 'I.II', 'I.III', 'I.IV', 'I.V', 'II.I', 'III.I', 'III.II', 'III.III', 'III.V', 'III.VI']:
        for cls in ('alpha', "alpha'", 'beta', "beta'"):
            keys = set(k for k in a3[g] if k[1] == cls) | set(k for k in a2[g] if k[1] == cls)
            d = max([abs(a3[g][k] - a2[g][k]) for k in keys] + [0])
            s = max([abs(a2[g][k]) for k in keys] + [abs(a3[g][k]) for k in keys] + [0])
            if d > 1e-9:
                # try ratio
                rs = [a3[g][k] / a2[g][k] for k in keys if abs(a2[g][k]) > 1e-9]
                print('%-8s %-7s differs: max diff %.3f (scale %.3f) ratios %s' % (g, cls, d, s, np.round(rs, 3)))
            else:
                print('%-8s %-7s same' % (g, cls))
