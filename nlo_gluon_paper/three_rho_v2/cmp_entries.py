"""Term-by-term comparison of the list as sent with the earlier entries (ThreeRho_cancellation, as written),
class by class, in common variables (symmetrized charges)."""
import sys
import numpy as np
import env
from core import symmetrize
import numerics as nm
import parse as P
import cancel_data as CD
from sources import SOURCES

NC = 3
SRC = {s.key: s for s in SOURCES}
LABS = ['x', "x'", 'y', 'w', "w'", 'z']


class Cfg:
    def __init__(self, seed):
        rng = np.random.default_rng(seed)
        self.ts, self.f = nm.structure_constants(NC)
        self.U = {l: nm.random_adjoint(NC, self.ts, rng, True) for l in LABS}
        self.r = {l: rng.normal(size=NC * NC - 1) for l in LABS}

    def val(self, fl):
        return nm.eval_colour(fl, self.f, self.U, self.r, NC * NC - 1)


CFG = [Cfg(3), Cfg(8)]


def fp_user(e):
    fl = [('U', i, j, p) for i, j, p in e['Us']] + [('f',) + tuple(f) for f in e['fs']] + [('rho', i, p) for i, p in e['rs']]
    c = e['pref'] * e['sign'] * e['mult']
    return np.array([c * cf.val(fl) for cf in CFG])


def fp_entry(en):
    c = complex(en['coeff'])
    return np.array([c * cf.val(en['col']) for cf in CFG])


def compare(blocks, keys, show=True):
    U = [e for e in P.parse() if e['block'] in blocks]
    O = []
    for k in keys:
        O += CD.source_entries(SRC[k] if k in SRC else CD.D_SOURCE)[0]
    out = []
    for cls in ['alpha', "alpha'", 'beta', "beta'", 'gamma', "gamma'"]:
        u = [(e, fp_user(e)) for e in U if e['cls'] == cls]
        o = [(en, fp_entry(en)) for en in O if en['cls'] == cls]
        tu = sum((f for _, f in u), np.zeros(2))
        to = sum((f for _, f in o), np.zeros(2))
        if not u and not o:
            continue
        used = set()
        unmatched_u = []
        for e, f in u:
            hit = None
            for j, (en, g) in enumerate(o):
                if j not in used and np.allclose(f, g, atol=1e-10):
                    hit = j
                    break
            if hit is None:
                unmatched_u.append((e, f))
            else:
                used.add(hit)
        unmatched_o = [(en, g) for j, (en, g) in enumerate(o) if j not in used]
        out.append((cls, tu, to, unmatched_u, unmatched_o))
        if show:
            print('  class %-7s user %2d terms, earlier %2d entries; class totals differ by %.2e' % (cls, len(u), len(o), np.abs(tu - to).max()))
            for e, f in unmatched_u:
                print('     user term %s.%d not matched: %s' % (e['block'], e['n'], e['text'][:150].replace('\n', ' ')))
                for en, g in unmatched_o:
                    r = f / g
                    if np.allclose(r, r[0], rtol=1e-8):
                        print('         = %s x earlier entry %s' % (np.round(r[0], 6), en['label']))
            for en, g in unmatched_o:
                print('     earlier entry %s not matched: coeff %s colour %s' % (en['label'], en['coeff'], en['col']))
    return out


if __name__ == '__main__':
    for blocks, keys in ((['III.I'], ['GIII.I']), (['III.II'], ['GIII.II'])):
        print(blocks, 'vs', keys)
        compare(blocks, keys)
