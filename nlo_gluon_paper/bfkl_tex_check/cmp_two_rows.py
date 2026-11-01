r"""Two-rho rows of BFKL.tex: the 'written as' form against your first form (texeval, pointwise, U = U^T, the c.c. added
for Group III; epsilon = 0.37).  Group II Row I second part, Group II Row IV, Group III Row III second part."""
import json, os, re, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmp_first as C
T3, TE = C.T3, C.TE
sys.path.insert(0, os.path.join(HERE, '..', 'evol_check'))
from cmp_three import kdot
from rows_user import rows

PAIRS = {'G2 I second part': ([295, 300, 304, 309, 313], [318], False),
         'G2 IV': ([345], [351], False),
         'G3 III second part': ([489], [498], True)}


def val(t, c, cc):
    v = TE.term_value(t, T3.env(c))
    if cc:
        v += TE.term_value(t, T3.env(T3.swapped_w(c))).conjugate()
    return v


def sym2(t, c, cc):
    """sum over the two placements of the charges x, x' (the measure is symmetric)"""
    out = 0
    for p in (('x', 'y', 'xp'), ('xp', 'y', 'x')):
        out += val(t, T3.permuted(c, p), cc)
    return out


if __name__ == '__main__':
    U = {r[3]: r[4] for r in rows('BFKL.tex')}
    out = {}
    for name, (first, written, cc) in PAIRS.items():
        f = [C.prep(C.plus(U[l])) for l in first]
        w = [kdot(C.prep(U[l])) for l in written]
        r = []
        for s in (1, 2, 3):
            c = T3.configs(s)
            a = sum(sym2(t, c, cc) for t in f)
            b = sum(sym2(t, c, False) for t in w)
            r.append((abs(a - b) / abs(a), complex(b / a)))
        out[name] = dict(rel=max(x[0] for x in r), ratio=[str(x[1]) for x in r])
        print(name, out[name])
    json.dump(out, open(os.path.join(HERE, 'cmp_two_rows.json'), 'w'), indent=1)
