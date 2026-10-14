"""Two-rho terms: literal charge order of the list (corrected) vs operator order (original expressions, c.c. reversed).
Pointwise comparison row by row, and the total with the genuine and four-rho two-rho terms vs JIMWLK x LO."""
import json
import sys
from collections import OrderedDict
import env
from core import symmetrize, Config, pointwise_total
import parse as P
import fixes as FX
import target as TG
sys.path.insert(0, env.ALL2)
import sets as S2   # ../two_rho_all/sets.py: operator-ordered two-rho terms, one label per source

ROW = {'4rho.a': 'A,B', '4rho.b': 'A,B', 'I.II.1': '3:I.II', 'I.II.2': '3:I.II', 'I.III': '3:I.III', 'D': 'D',
       'I.IV': '3:I.IV', 'I.V.1': '3:I.V', 'I.V.2': '3:I.V', 'II.I': '3:II.I', 'III.I': '3:III.I', 'III.II': '3:III.II',
       'III.III': '3:III.III', 'III.V': '3:III.V', 'III.VI': '3:III.VI'}


def literal_sets():
    out = OrderedDict()
    for e in FX.corrected():
        out.setdefault(ROW[e['block']], []).extend(symmetrize(P.as_term(e), (2,)))
    return out


if __name__ == '__main__':
    L = literal_sets()
    A = S2.A_sets()
    res = dict(rows={}, totals=[])
    cfg = Config(3, True)
    for k in L:
        a = pointwise_total(L[k], cfg)[0].real
        b = pointwise_total(A[k], cfg)[0].real
        res['rows'][k] = dict(literal=a, operator=b)
        print('%-10s literal %10.5f   operator order %10.5f   %s' % (k, a, b, 'same' if abs(a - b) < 1e-9 * (1 + abs(b)) else 'DIFFERENT'))
    other = [t for k in ('g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III', '4rho') for t in A[k]]
    J = TG.jimwlk()
    for seed in (3, 7, 11):
        cfg = Config(seed, True)
        lit = sum(pointwise_total(ts, cfg)[0] for ts in L.values()).real + pointwise_total(other, cfg)[0].real
        op = sum(pointwise_total(A[k], cfg)[0] for k in L).real + pointwise_total(other, cfg)[0].real
        j = pointwise_total(J, cfg)[0].real
        res['totals'].append(dict(seed=seed, literal=lit, operator=op, J=j))
        print('seed %d: total two-rho (literal) %.6f   (operator order) %.6f   JIMWLK x LO %.6f' % (seed, lit, op, j))
    json.dump(res, open('two_cmp.json', 'w'), indent=1)
