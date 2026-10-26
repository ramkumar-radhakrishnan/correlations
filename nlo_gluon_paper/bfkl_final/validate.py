"""Check of the transcription in rows.py: every row (no renaming, no c.c., classical charges) against the independent
numerical transcription of the same row (../finite_v2/split.py, ../bfkl_check/corrected_lvk.py, user_lvk.py)."""
import json
import numpy as np
import env  # noqa: F401
import split as SP
import corrected_lvk as CO
import user_lvk as UL
import three_rho_test as TT
from fp import G, PI
from ft import PI3
import rows as R

NAME = {"x'": 'xp', "w'": 'wp', "y'": 'yp'}


def kval(f, P):
    a, b, c, d = (P[NAME.get(p, p)] for p in f[1:])
    A, B = a - b, c - d
    return (A @ B) / ((A @ A) * (B @ B))


def value(terms, cfg):
    tot = 0.0
    for t in terms:
        k = 1.0
        for f in t['kern']:
            k *= kval(f, cfg.pt)
        fl = [('U', i, j, NAME.get(p, p)) for i, j, p in t['Us']] + [('f',) + tuple(f) for f in t['fs']] + \
             [('rho', i, NAME.get(p, p)) for i, p in t['rs']]
        tot += t['c'] * k * nm.eval_colour(fl, F, cfg.U, cfg.rho, 8)
    return tot


import numerics as nm  # noqa: E402
TS, F = nm.structure_constants(3)
UNIT = G**4 / (16 * PI**5) / PI3          # N without log, k+ = 1


def no_rename(fn):
    old = R.build

    def b(coef, kern, colour, tag, mp=None, cc=False):
        return old(coef, kern, colour, tag, None, False)
    R.build = b
    try:
        return fn()
    finally:
        R.build = old


REF = {
    'G1 II 2a': (R.g1_ii_2a, lambda c: TT.g1_ii_2a(c)),
    'G1 III (F3)': (R.g1_iii, lambda c: TT.g1_iii_F3(c)),
    'G1 III (as written)': (lambda: R.g1_iii(F3=False), lambda c: SP.C_g1_iii(c, 1.0)),
    'G1 IV': (R.g1_iv, lambda c: SP.C_g1_iv(c, 1.0)),
    'G1 V(1)': (R.g1_v1, lambda c: SP.C_g1_v1(c, 1.0)),
    'G2 I(1)': (R.g2_i1, lambda c: SP.C_g2_i1(c, 1.0)),
    'G3 I': (R.g3_i, lambda c: SP.C_g3_i(c, 1.0)),
    'G3 II': (R.g3_ii, lambda c: SP.C_g3_ii(c, 1.0)),
    'G3 III(1)': (R.g3_iii1, lambda c: SP.ROWS['G3-III(1)'][1](c, 1.0)),
    'G3 V': (R.g3_v, lambda c: SP.ROWS['G3-V'][1](c, 1.0)),
    'G3 VI': (R.g3_vi, lambda c: SP.ROWS['G3-VI'][1](c, 1.0)),
    'G1 V(2a)': (R.g1_v2a, lambda c: UL.g1_v2a(c, 1.0)),
    'G1 V(2b fin)': (R.g1_v2b_fin, lambda c: UL.g1_v2b_fin(c, 1.0)),
    'G2 I(2)': (R.g2_i2, lambda c: CO.g2_i2(c, 1.0)),
    'G3 III(2)': (R.g3_iii2, lambda c: CO.g3_iii2(c, 1.0)),
}

if __name__ == '__main__':
    out = {}
    for name, (eng, ref) in REF.items():
        ts = no_rename(eng)
        rel = []
        for s in (1, 2, 3):
            c = SP.cfg(s)
            c.kv = np.zeros(2)
            a = value(ts, c) * UNIT
            b = complex(ref(c))
            rel.append(abs(a - b) / abs(b))
        out[name] = max(rel)
        print('%-22s engine vs numerical row: %.1e' % (name, max(rel)))
    json.dump(out, open('validate.json', 'w'), indent=1)
