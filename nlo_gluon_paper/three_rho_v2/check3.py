"""Symmetrized three-rho parts of the list as sent: do they cancel?  Row-by-row comparison with the earlier sources."""
import sys
import numpy as np
from collections import defaultdict, OrderedDict
import env
from core import symmetrize
import acc3 as A3
import parse as P

sys.path.insert(0, env.COMB)
import allsets as AS  # earlier rows (combined note)

E = P.parse()
ROWMAP = {'4rho.a': 'A,B', '4rho.b': 'A,B', 'I.II.1': 'I.II', 'I.II.2': 'I.II', 'I.III': 'I.III', 'I.IV': 'I.IV',
          'I.V.1': 'I.V', 'I.V.2': 'I.V', 'II.I': 'II.I', 'III.I': 'III.I', 'III.II': 'III.II', 'III.III': 'III.III',
          'III.V': 'III.V', 'III.VI': 'III.VI'}


def user_sets(entries=E, kern='class'):
    out = OrderedDict()
    for e in entries:
        out.setdefault(ROWMAP[e['block']], []).extend(symmetrize(P.as_term(e, kern), (3,)))
    return out


def totals(sets, cfg):
    acc = defaultdict(lambda: defaultdict(complex))
    for k, ts in sets.items():
        A3.accumulate3(ts, cfg, acc, label=k)
    return acc


if __name__ == '__main__':
    U = user_sets()
    old_w = AS.A_sets((3,), False)
    old_c = AS.A_sets((3,), True)
    for seed in (5, 17):
        cfg = A3.Config3(seed, True)
        acc = totals(U, cfg)
        bad = A3.report3(acc, show=False)
        print('seed', seed, ': kernel keys', len(acc), ' non-cancelling keys', len(bad))
        for a, sc, k, d in bad:
            print('   %10.5f  %s' % (a, k))
    cfg = A3.Config3(5, True)
    print('\nrow by row (sum over kernel keys of |user - earlier|):')
    for row in ['A,B', 'I.II', 'I.III', 'I.IV', 'I.V', 'II.I', 'III.I', 'III.II', 'III.III', 'III.V', 'III.VI', 'D']:
        a = totals({'u': U.get(row, [])}, cfg)
        for lab, old in (('as written', old_w), ('corrected', old_c)):
            b = totals({'o': old.get(row, [])}, cfg)
            keys = set(a) | set(b)
            diff = sum(abs(sum(a[k].values()) - sum(b[k].values())) for k in keys)
            size = sum(abs(sum(b[k].values())) for k in keys)
            print('  %-8s vs earlier %-10s  diff %.3e   (size %.3e)' % (row, lab, diff, size))
