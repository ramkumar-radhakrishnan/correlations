r"""Pointwise check: every equation of Limiting_cases.tex = the limit k+/vee -> 0 (Lambda -> 0) of the checked
exact-vee equation (../neglog/corrected/Other_terms.tex).  Each equation is evaluated as a sum of its lines (uv_check
machinery: smooth SU(3) Wilson lines and charges, longitudinal variable fixed at p+ = 0.37 k+, plus prescriptions
resolved), the old one at vee = k+/a for a = 1e-2 ... 1e-8, the new one as written.
  python3 limcheck.py  -> limcheck.json"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'uv_check'))
import uvscan as UV  # noqa: E402
import terms as TM  # noqa: E402
from limblocks import blocks  # noqa: E402

OLD = blocks(os.path.join(HERE, '..', 'neglog', 'corrected', 'Other_terms.tex'))
NEW = blocks(os.path.join(HERE, 'input', 'Limiting_cases.tex'))
# (old index, new index, label)
PAIRS = [(2, 3, 'G1 III part 1'), (3, 4, 'G1 III part 2'), (4, 5, 'G1 III part 3a'), (5, 6, 'G1 III part 3b'),
         (7, 8, 'G1 IV'), (8, 9, 'G1 V first part'), (9, 10, 'G1 V 2a'), (11, 11, 'G1 V 2b finite'),
         (12, 12, 'G2 I 1a'), (13, 13, 'G2 I 1b'), (14, 14, 'G2 I 2a'), (15, 15, 'G2 I 2b'), (16, 16, 'G2 I 2c'),
         (17, 17, 'G2 I 2d'), (18, 18, 'G2 I 2e'), (19, 19, 'G2 I 2f'),
         (20, 21, 'G3 I'), (21, 22, 'G3 II'), (22, 23, 'G3 III part 1'), (23, 24, 'G3 III part 2'),
         (24, 25, 'G3 V'), (25, 26, 'G3 VI')]
AS = [1e-2, 1e-4, 1e-6, 1e-8]


def eq_terms(block):
    # texeval's LOG is log|.|, so the modulus bars only need to become brackets
    block = block.replace(r'\left|', r'\left(').replace(r'\right|', r'\right)')
    return [UV.Term(t) for t in TM.equation_terms(block)]


def total(ts, cfg, pl):
    return sum(np.asarray(t.value(cfg, pl), dtype=complex) for t in ts)


def main():
    out = {}
    for io, jn, lab in PAIRS:
        to, tn = eq_terms(OLD[io][1]), eq_terms(NEW[jn][1])
        rows = []
        for seed in (0, 1, 2):
            field = UV.Field(seed)
            pts, kv = UV.base_points(seed)
            cfg = UV.Cfg(field, pts, kv)
            UV.VEE = 1.0e30
            vn = total(tn, cfg, 0.37)
            d = []
            for a in AS:
                UV.VEE = UV.KP / a
                vo = total(to, cfg, 0.37)
                d.append(float(np.abs(vo - vn).max() / max(np.abs(vn).max(), 1e-300)))
            rows.append({'seed': seed, 'new': float(np.abs(vn).max()), 'reldiff': d})
        out[lab] = rows
        print('%-16s' % lab, ' | '.join(' '.join('%.1e' % x for x in r['reldiff']) for r in rows), flush=True)
    json.dump(out, open(os.path.join(HERE, 'limcheck.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
