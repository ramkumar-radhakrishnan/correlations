"""Simplified (Lambda -> 0, vee -> infinity) forms of the closed-form finite rows, with the large logarithm
L_v = log(vee/k+) separated:   finite_row(vee, Lambda) = C_row * L_v + R_row + O(k+/vee, Lambda/k+).
Each function returns (C, R) at the transverse points of c (same normalization, colour and phase as ft.py).
check() compares with ft.py at vee/k+ = 1e6..1e9, Lambda/k+ = 1e-7..1e-10."""
import numpy as np
import env  # noqa: F401
from fp import G, PI, I2, sq
import ft
from ft import vec, Kv, ij, br4, brlog, AB, PI3


def g1_iv(c, kp):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    A, B = AB(c)
    p8 = 1j * G**4 / (8 * PI**5) / kp
    C = p8 * (ij(Xp, Yz, br4(c)) - 0.5 * (Xp @ Yz) / sq(vec(c, 'z', 'w')))
    R = p8 * ij(Xp, Yz, br4(c)) * np.log(A / B)
    pc = ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_IV(c)
    return C * pc, R * pc


def g3_kernel(c, kp):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    A, B = AB(c)
    p8 = 1j * G**4 / (8 * PI**5) / kp
    C = -p8 * (ij(Xp, Yz, br4(c)) - 0.5 * (Xp @ Yz) / sq(vec(c, 'z', 'w')))
    R = -p8 * ij(Xp, Yz, br4(c)) * np.log(A / B)
    return C, R


def g3_i(c, kp):
    C, R = g3_kernel(c, kp)
    pc = ft.phase(c, 'w', 'wp') / PI3 * ft.col_G3_I(c)
    return C * pc, R * pc


def g3_ii(c, kp):
    C, R = g3_kernel(c, kp)
    pc = ft.phase(c, 'w', 'wp') / PI3 * ft.col_G3_II(c)
    return C * pc, R * pc


def g2_i1(c, kp):
    p16 = 1j * G**4 / (16 * PI**5) / kp
    Xp, Yz, Xpw = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'xp', 'w')
    A, B = AB(c, 'x', 'w')
    A2, B2 = AB(c, 'x', 'wp')
    CA = -p16 * (ij(Xp, Yz, br4(c, 'w')) - 0.5 * (Xp @ Yz) / sq(vec(c, 'z', 'w')))
    RA = -p16 * ij(Xp, Yz, br4(c, 'w')) * np.log(A / B)
    CB = p16 * (ij(Xpw, Yz, br4(c, 'wp')) - 0.5 * (Xpw @ Yz) / sq(vec(c, 'z', 'wp')))
    RB = p16 * ij(Xpw, Yz, br4(c, 'wp')) * np.log(A2 / B2)
    ph = ft.phase(c, 'w', 'wp') / PI3
    ca, cb = ft.col_G2_I1(c, 'w', 'wp'), ft.col_G2_I1(c, 'wp', 'w')
    return ph * (CA * ca + CB * cb), ph * (RA * ca + RB * cb)


def g1_iii(c, kp):
    A, B = AB(c)
    xz, xw, zw = vec(c, 'x', 'z'), vec(c, 'x', 'w'), vec(c, 'z', 'w')
    Xp, Yz, Zw, Xw, Xz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'z', 'w'), Kv(c, 'x', 'w'), Kv(c, 'x', 'z')
    p16, p8 = 1j * G**4 / (16 * PI**5) / kp, 1j * G**4 / (8 * PI**5) / kp
    XY = Xp @ Yz
    lAB, lA, lB = np.log(A / B), np.log(np.abs(A / (A - B))), np.log(np.abs((B - A) / B))
    # N1
    C1 = -p8 * (Xp @ Xw) * (Xz @ Yz)
    R1 = p16 * XY * (A + B) / (A - B)**2 * lAB - p8 * (Xp @ Xw) * (Xz @ Yz) * lAB
    # N2 = -n2
    S6 = ((Xp @ Zw) * (Yz @ xw) / (B - A) - (Xp @ Xw) * (Yz @ Zw) * A / (B - A) + XY * (zw @ xw) * A / (sq(zw) * (B - A)**2)
          - (Xp @ Zw) * (Yz @ Xz) * B / (B - A) + (Xp @ xz) * (Yz @ Zw) / (B - A) - XY * (zw @ xz) * B / (sq(zw) * (B - A)**2))
    C2 = -(-p8 * (Xp @ Xw) * (Yz @ Zw) + p8 * S6)
    R2 = -(p8 * S6 * lA)
    # N3a
    T3 = (XY / sq(zw) + XY * (zw @ xw) * A / (sq(zw) * (B - A)**2) - XY * (zw @ xz) * B / (sq(zw) * (B - A)**2))
    C3 = p8 * T3
    root = np.sqrt(A * B - (xw @ xz)**2)
    A0 = (np.arctan((A - xw @ xz) / root) - np.arctan(-(xw @ xz) / root)) / root
    R3 = (p8 * lB * (-XY * (zw @ xw) * A / (sq(zw) * (B - A)**2) + XY * (zw @ xz) * B / (sq(zw) * (B - A)**2)
                     - (Xp @ Zw) * (Yz @ xw) / (B - A) + (Xp @ Zw) * (Yz @ Xz) * B / (B - A)
                     + (Yz @ Zw) * (Xp @ Xw) * A / (B - A) - (Yz @ Zw) * (Xp @ xz) / (B - A))
          + p8 * A0 * ((Xp @ Zw) * (Yz @ xw) + (Yz @ Zw) * (Xp @ xz) - (xw @ xz) * (Xp @ Zw) * (Yz @ xz) / A))
    # N3b
    L4 = np.log(sq(zw) / B)
    R4 = (-p8 * A0 * (xw @ xz) * (Yz @ Zw) * (Xp @ xw) / B
          - p16 * L4 * (Xp @ Zw) * (Yz @ xz) / A + p16 * L4 * (Yz @ Zw) * (Xp @ xw) / B - p16 * L4 * XY / sq(zw))
    pc = ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_III(c)
    return (C1 + C2 + C3) * pc, (R1 + R2 + R3 + R4) * pc


def g2_i2_groups(c):
    """br2a .. br2e as in ft.G2_I2"""
    x, xp, z, w, wp = (c.pt[n] for n in ('x', 'xp', 'z', 'w', 'wp'))
    xz, xw, zw, xpz, xpwp, zwp = x - z, x - w, z - w, xp - z, xp - wp, z - wp
    A, B, Ap, Bp, zw2, zwp2 = sq(xz), sq(xw), sq(xpz), sq(xpwp), sq(zw), sq(zwp)
    br = {}
    br['a'] = ((xw @ zwp) * (xpz @ zw) / (zw2 * zwp2) - (xw @ xpwp) * (zw @ xpz) / (2 * zw2 * Bp)
               + (xz @ xpz) * (xw @ zwp) / (2 * A * zwp2) - (xpz @ xz) * (xpwp @ xw) / (4 * A * Bp)
               + (xpwp @ zw) * (xz @ zwp) / (zw2 * zwp2) + (xpwp @ zw) * (xz @ xpz) / (2 * zw2 * Ap)
               - (xz @ zwp) * (xpwp @ xw) / (2 * B * zwp2) - (xpwp @ xw) * (xpz @ xz) / (4 * B * Ap))
    br['b'] = ((xpwp @ xw) * (zw @ zwp) / (zw2 * zwp2) + (xpwp @ xw) * (xpz @ zw) / (2 * zw2 * Ap)
               + (xw @ xpwp) * (xz @ zwp) / (2 * A * zwp2) + (xz @ xpz) * (xpwp @ xw) / (4 * A * Ap))
    br['c'] = ((xpz @ xz) * (zw @ zwp) / (zw2 * zwp2) - (xpz @ xz) * (xpwp @ zw) / (2 * zw2 * Bp)
               - (xw @ zwp) * (xz @ xpz) / (2 * B * zwp2) + (xz @ xpz) * (xpwp @ xw) / (4 * B * Bp))
    br['d'] = ((xw @ zw) * (Ap - Bp) / (2 * zw2 * zwp2) + (xz @ xw) * (Ap - Bp) / (4 * A * zwp2)
               + (xpwp @ zwp) * (A - B) / (2 * zw2 * zwp2) + (xpz @ xpwp) * (A - B) / (4 * zw2 * Ap))
    br['e'] = ((xz @ zw) * (Ap - Bp) / (2 * zw2 * zwp2) - (xz @ xw) * (Ap - Bp) / (4 * B * zwp2)
               + (xpz @ zwp) * (A - B) / (2 * zw2 * zwp2) - (xpz @ xpwp) * (A - B) / (4 * zw2 * Bp))
    return br, (A, B, Ap, Bp, zw2, zwp2)


def g2_i2(c, kp, d=2):
    br, (A, B, Ap, Bp, zw2, zwp2) = g2_i2_groups(c)
    D = Bp * A - B * Ap
    lX, lXp = np.log(A / B), np.log(Ap / Bp)
    pre = G**4 / (8 * PI**5) / kp
    pf = G**4 * d / (32 * PI**5)
    C = {'2a': 0, '2b': -pre * br['b'] / (B * Bp), '2c': pre * br['c'] / (A * Ap), '2d': 0, '2e': 0, '2f': 0}
    R = {'2a': pre * br['a'] / D * (np.log(A / Ap) - np.log(B / Bp)),
         '2b': pre * br['b'] * (-A * lX / (B * D) + Ap * lXp / (Bp * D)),
         '2c': pre * br['c'] * (Bp / Ap * lXp - B / A * lX) / D,
         '2d': pre * br['d'] * (A * lX / ((A - B) * D) + Ap * lXp / ((Bp - Ap) * D)),
         '2e': pre * br['e'] * (-B * lX / ((A - B) * D) - Bp * lXp / ((Bp - Ap) * D)),
         '2f': pf / kp / (zw2 * zwp2) * (-A * B * (Ap - Bp) / ((A - B) * D) * lX + Bp * Ap * (A - B) / ((Ap - Bp) * D) * lXp - 1)}
    pc = ft.phase(c, 'w', 'wp') / PI3 * ft.col_G2_I2(c)
    return {k: C[k] * pc for k in C}, {k: R[k] * pc for k in R}


ROWS = {
    'G1-III': (g1_iii, lambda c, k, V, L: sum(ft.G1_III(c, k, V, L).values())),
    'G1-IV': (g1_iv, ft.G1_IV),
    'G2-I(1)': (g2_i1, lambda c, k, V, L: sum(ft.G2_I1(c, k, V, L))),
    'G2-I(2)': (lambda c, k: tuple(sum(d.values()) for d in g2_i2(c, k)), lambda c, k, V, L: sum(ft.G2_I2(c, k, V, L).values())),
    'G3-I': (g3_i, ft.G3_I),
    'G3-II': (g3_ii, ft.G3_II),
}


def check(seeds=(1, 2, 3, 4)):
    import fp
    out = {}
    for name, (s, f) in ROWS.items():
        worst = []
        for seed in seeds:
            c = fp.Cfg(seed, sym=True, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp'))
            kp = 1.0
            C, R = s(c, kp)
            errs = []
            for V, L in ((1e6, 1e-7), (1e9, 1e-10)):
                fin = f(c, kp, V, L)
                errs.append(abs(fin - (C * np.log(V / kp) + R)) / abs(fin))
            worst.append(errs)
        out[name] = worst
        print('%-8s rel. difference finite - (C L_v + R) at vee/k = 1e6, 1e9:' % name,
              ['%.1e/%.1e' % tuple(e) for e in worst])
    return out


if __name__ == '__main__':
    import json
    json.dump(check(), open('simp_check.json', 'w'), indent=1)
