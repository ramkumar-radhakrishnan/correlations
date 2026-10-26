"""Independent cross-check of the extraction.  For every row with no hidden logarithm:
   (coefficient of log(vee/k+)) - (coefficient of log(k+/Lambda))  [new grouping, from the J and B parts generated here]
 = coefficient of log(vee/k+) in the old grouping (log(vee/Lambda) + log(vee/k+)), computed earlier from the rows
   (../finite_v2/split.py, C_row; checked there against the numerical vee dependence of every row).
Output: cross_check.json"""
import json
import os
import sys
import numpy as np
import blocks as BL
import gen as GN
import texeval as TE
sys.path.insert(0, os.path.normpath(os.path.join(TE.HERE, '..', 'finite_v2')))
sys.path.insert(0, os.path.normpath(os.path.join(TE.HERE, '..', 'finite')))
import env  # noqa
import split as SP

BLK, SPL = GN.all_splits()
E = np.e


def coef(names, part, c, kp):
    """sum over the blocks of the J (or B) parts, with the large log set to 1"""
    env = TE.Env(c, KP=kp, V=E * kp, L=kp / E)
    tot = 0.0
    for name in names:
        B = BLK[name]
        for n, s in enumerate(SPL[name]):
            t = s.get(part)
            if t is None:
                continue
            tot += TE.term_value(B.prefix + ' (' + t + ') ' + B.suffix, env)
    return tot


ROWS = {
    'G1 III (as written)': (['G1 III N1', 'G1 III N2', 'G1 III N3a', 'G1 III N3b'], SP.C_g1_iii),
    'G1 IV': (['G1 IV'], SP.C_g1_iv),
    'G2 I(1)': (['G2 I(1) A', 'G2 I(1) B'], lambda c, kp: sum(np.atleast_1d(SP.C_g2_i1(c, kp)))),
    'G2 I(2)': (['G2 I(2a)', 'G2 I(2b)', 'G2 I(2c)', 'G2 I(2d)', 'G2 I(2e)', 'G2 I(2f)'], SP.C_g2_i2),
    'G3 I': (['G3 I'], SP.C_g3_i),
    'G3 II': (['G3 II'], SP.C_g3_ii),
}

if __name__ == '__main__':
    out = {}
    for row, (names, C) in ROWS.items():
        rel = []
        for seed in (1, 2, 3):
            c = SP.cfg(seed)
            kp = 1.0
            J, B = coef(names, 'J', c, kp), coef(names, 'B', c, kp)
            ref = C(c, kp)
            if isinstance(ref, dict):
                ref = sum(ref.values())
            ref = complex(np.sum(ref))
            rel.append(abs((B - J) - ref) / max(abs(ref), 1e-30))
        out[row] = max(rel)
        print('%-22s  (B - J) vs C_row: %.1e' % (row, max(rel)))
    json.dump(out, open('cross_check.json', 'w'), indent=1)
