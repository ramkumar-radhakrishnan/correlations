"""Symmetrized three-rho terms of bfkl.tex (= the corrected log(vee/k+) list), pointwise in z.

Direct set: every three-rho term whose colour is CL^a(x',w') x (two charges on the amplitude side): Group I rows, Group III
rows, Group II Row I part 1 block A.  (The c.c.'s and block B are the mirror images and cancel in the same way.)
All rows are mapped to the common variables x', w' (conjugate LO), w (measured gluon), x, y (charges), z (soft gluon);
the phase is then e^{-ik(w'-w)} for all of them (kv = 0 below).  Classical charges give the symmetric part S(O).
Group I Row II 2a is written as a z integral:  X'.F(w;x,y) = (1/pi) int_z F_z."""
import copy
import json
import numpy as np
import env  # noqa: F401
from fp import G, PI, sq, F as FF
import ft
import split as SP
from ft import Kv, PI3

CAN = ('x', 'y', 'xp', 'wp', 'z', 'w')


def mapped(c, m):
    """configuration for a row: row variable r takes the point/U/rho of canonical variable m[r]"""
    d = copy.deepcopy(c)
    for r, cn in m.items():
        d.pt[r], d.U[r], d.rho[r] = c.pt[cn].copy(), c.U[cn].copy(), c.rho[cn].copy()
    return d


def swap_xy(c):
    return mapped(c, {'x': 'y', 'y': 'x'})


ID = {n: n for n in CAN}
kp = 1.0
p8 = 1j * G**4 / (8 * PI**5) / kp
p16 = p8 / 2
p32 = p8 / 4


def Fz(c):
    x, y, z, w = (c.pt[n] for n in ('x', 'y', 'z', 'w'))
    K = lambda v: v / sq(v)
    Xp = K(c.pt['xp'] - c.pt['wp'])
    return ((Xp @ K(x - w)) * (K(y - z) @ (K(x - z) - K(w - z)))
            - (Xp @ K(y - w)) * (K(x - z) @ (K(y - z) - K(w - z))))


def g1_ii_2a(c):
    ca, cb = ft.col_G1_II_2a(c)
    return 1j * G**4 / (32 * PI**5) / kp * Fz(c) * (-ca + cb) / PI3


def g1_v1_t1(c):
    # term 1, phase e^{-ik(w'-z)}: row z = measured gluon, row x = soft gluon, row y' = charge
    d = mapped(c, {'z': 'w', 'x': 'z', 'yp': 'x', 'y': 'y', 'xp': 'xp', 'wp': 'wp'})
    kk = ft.G1_V1_kernels(d)
    return p32 / PI3 * kk['b'] * ft.col_G1_V1(d)


def g1_v1_t2(c):
    # term 2, phase e^{-ik(w'-x)}: row x = measured gluon, row z = soft gluon, row y' = charge
    d = mapped(c, {'x': 'w', 'z': 'z', 'yp': 'x', 'y': 'y', 'xp': 'xp', 'wp': 'wp'})
    kk = ft.G1_V1_kernels(d)
    return p32 / PI3 * kk['q'] * ft.col_G1_V1(d)


def g2_i1_A(c):
    return -p16 * SP.Ck_iv(c, 'w', 'x', 'wp') * ft.col_G2_I1(c, 'w', 'wp') / PI3


def crow(name):
    """the c.c. of Group III Rows III(1), V, VI: row y' -> w (amplitude gluon), row w -> w', row y <-> x'"""
    f = SP.ROWS[name][1]
    def g(c):
        d = mapped(c, {'yp': 'w', 'w': 'wp', 'y': 'xp', 'xp': 'y', 'x': 'x', 'z': 'z'})
        return np.conj(f(d, kp))
    return g


def col_III_w(c):
    """O_III with U^{bc}(z) -> U^{bc}(w)"""
    U, r = c.U, c.rho
    return np.einsum('abd,a,b,d->', FF, ft.CL(c), U['w'] @ r['y'] - U['y'] @ r['y'], U['x'] @ r['x'])


def kerns(c):
    K = lambda v: v / sq(v)
    x, y, z, w = (c.pt[n] for n in ('x', 'y', 'z', 'w'))
    Xp, Y, Zw, Xw, Xz, Yw = K(c.pt['xp'] - c.pt['wp']), K(y - z), K(z - w), K(x - w), K(x - z), K(y - w)
    return dict(Kh=(Xp @ Zw) * (Y @ Xz), Kb=(Xp @ Xw) * (Y @ Xz), Kf=(Xp @ Y) / sq(z - w), D=(Xp @ Yw) * (Xz @ K(w - z)))


def g1_iii_simpl(c):
    k = kerns(c)
    return p8 * (k['Kh'] - k['Kb'] + k['Kf']) * ft.col_G1_III(c) / PI3


def g1_iii_F3(c):
    """F3: the K_f O_III term replaced by -D:  p8[(K_h - K_b) O_III + X'.K(y-w) K(x,w,z) O_III(z->w)]"""
    k = kerns(c)
    return p8 * ((k['Kh'] - k['Kb']) * ft.col_G1_III(c) + k['D'] * col_III_w(c)) / PI3


def row_fn(name):
    return lambda c: SP.ROWS[name][1](c, kp)


DIRECT = {
    'G1 II 2a': g1_ii_2a,
    'G1 III': row_fn('G1-III'),
    'G1 IV': row_fn('G1-IV'),
    'G1 V(1) t1': g1_v1_t1,
    'G1 V(1) t2': g1_v1_t2,
    'G2 I(1) A': g2_i1_A,
    'G3 I': row_fn('G3-I'),
    'G3 II': row_fn('G3-II'),
    'G3 III(1) c.c.': crow('G3-III(1)'),
    'G3 V c.c.': crow('G3-V'),
    'G3 VI c.c.': crow('G3-VI'),
}


def cfg(seed):
    c = SP.cfg(seed)
    c.kv = np.zeros(2)
    return c


def total(c, rows):
    vals = {k: complex(f(c)) + complex(f(swap_xy(c))) for k, f in rows.items()}
    return sum(vals.values()), sum(abs(v) for v in vals.values()), vals


if __name__ == '__main__':
    out = {}
    # 0. the G1 III kernel of the file simplifies to p8 [K_h - K_b + K_f]
    rel = max(abs(g1_iii_simpl(cfg(s)) - SP.C_g1_iii(cfg(s), kp)) / abs(SP.C_g1_iii(cfg(s), kp)) for s in range(1, 6))
    out['G1 III = p8(Kh - Kb + Kf)'] = rel
    print('G1 III kernel = p8 (K_h - K_b + K_f): rel. diff %.1e' % rel)
    for label, rows in (('as in bfkl.tex', DIRECT), ('G1 III with F3 (K_f O_III -> -D)', dict(DIRECT, **{'G1 III': g1_iii_F3}))):
        res = []
        for s in range(1, 9):
            t, scale, _ = total(cfg(s), rows)
            res.append(abs(t) / scale)
        out[label] = res
        print('%-36s  |sum|/sum|terms| over 8 configurations: max %.1e  min %.1e' % (label, max(res), min(res)))
    json.dump(out, open('three_rho_test.json', 'w'), indent=1)
