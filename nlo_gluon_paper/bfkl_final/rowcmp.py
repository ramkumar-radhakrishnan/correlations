"""Row by row: the two-rho terms of bfkl.tex (with F3) against the earlier verified log(vee/k+) set
(../bfkl/lsets.py, corrected: region-P limit minus the row's log(vee/Lambda) terms), key by key."""
import json
import os
import sys
from collections import defaultdict
import numpy as np
import env  # noqa: F401
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('bfkl', os.path.join('bfkl', 'reference')):
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', p)))
import lsets  # noqa: E402
import two_rho_engine as TE  # noqa: E402

MAP = {'G1 II 2a': 'I.II', 'G1 II(1)': 'I.II', 'G1 III': 'I.III', 'G1 IV': 'I.IV', 'G1 V(1)': 'I.V', 'G1 V(2a)': 'I.V',
       'G1 V(2b UV)': 'I.V', 'G1 V(2b fin)': 'I.V', 'G2 I(1)': 'II.I', 'G2 I(2)': 'II.I', 'G2 IV': 'II.IV',
       'G3 I': 'III.I', 'G3 II': 'III.II', 'G3 III(1)': 'III.III', 'G3 III(2)': 'III.III', 'G3 V': 'III.V', 'G3 VI': 'III.VI'}


def old_rows():
    S = lsets.sets(keep=(2,), corrected=True)
    out = defaultdict(list)
    for lab, ts in S.items():
        out[lab.rsplit(' ', 1)[0]] += ts
    return out


def new_rows(F3=True):
    FT = TE.file_terms(F3)
    out = defaultdict(list)
    for k, ts in FT.items():
        out[MAP[k]] += ts
    return out


if __name__ == '__main__':
    O, N = old_rows(), new_rows()
    J, B4, Bref = TE.targets()
    tot = {'old': [t for ts in O.values() for t in ts], 'B-J': Bref + TE.neg(J[-1]['terms'])}
    res = {}
    for seed in (1, 2):
        acc = TE.keyed({'old': tot['old'], '-(B-J)': TE.neg(tot['B-J'])}, seed)
        res.setdefault('old total vs B-J (bad keys)', []).append(len(TE.residual(acc)))
    print(res)
    for row in sorted(set(O) | set(N)):
        bad_s = []
        for seed in (1, 2):
            acc = TE.keyed({'new': N.get(row, []), 'old': TE.neg(O.get(row, []))}, seed)
            bad = TE.residual(acc)
            bad_s.append(len(bad))
            if seed == 1 and bad:
                for k, a, s, d in sorted(bad, key=lambda b: -b[1])[:6]:
                    print('      %-8s %.3e %s  %s' % (row, a, k, {l: np.round(v, 3) for l, v in d.items()}))
        res[row] = bad_s
        print('%-8s bad keys new vs old: %s' % (row, bad_s))
    json.dump(res, open('rowcmp.json', 'w'), indent=1)
