"""The corrected log(vee/k+) list (Finite_terms_1_1_corrected.tex), transcribed term by term as written there,
and checked row by row against the exact coefficient of log(vee/k+) of the same row of Finite_terms_1.tex
(C_row of ../finite_v2/split.py, and d/dlog(vee) of the file itself at vee/k+ = 1e7..1e9).

Same normalization as user_lvk.py: coefficient of L_v = log(vee/k+), with 1/(2 pi)^3, the prefactor, the phase and
the colour, without c.c.  Every kernel is written out explicitly here (no helper brackets), so that this file checks
the LaTeX, not the closed forms of ../finite_v2.  Output: corrected_check.json"""
import json
import numpy as np
import env  # noqa: F401
from fp import G, PI, sq
import ft
import split as SP
import user_lvk as UL
import rowcheck as RC

P8 = UL.P8
P16 = UL.P16
P32 = UL.P32


def d(a, b):
    return a @ b


def pts(c):
    return {k: c.pt[k] for k in c.pt}


def T_contract(c, Xp, Y, w='w'):
    """X'^i Y^j [delta^{ij}/(2(z-w)^2) + (x-z)^j (z-w)^i/((z-w)^2 (x-z)^2) - (x-z)^j (x-w)^i/((x-w)^2 (x-z)^2)
                 - (x-w)^i (z-w)^j/((z-w)^2 (x-w)^2)]"""
    p = pts(c)
    x, z, W = p['x'], p['z'], p[w]
    xz, xw, zw = x - z, x - W, z - W
    return (d(Xp, Y) / (2 * sq(zw)) + d(Xp, zw) * d(Y, xz) / (sq(zw) * sq(xz))
            - d(Xp, xw) * d(Y, xz) / (sq(xw) * sq(xz)) - d(Xp, xw) * d(Y, zw) / (sq(zw) * sq(xw)))


# ---------------------------------------------------------------- Group I
def g1_iii(c, kp):
    p = pts(c)
    x, y, z, w, xp, wp = p['x'], p['y'], p['z'], p['w'], p['xp'], p['wp']
    xpwp, yz, xw, xz, zw = xp - wp, y - z, x - w, x - z, z - w
    A, B = sq(xz), sq(xw)
    # d^3N_1: blue  -(i g^4 f/8pi^5)(1/k+) L_v (x'-w').(x-w) (x-z).(y-z) / ((x'-w')^2 (x-w)^2 (x-z)^2 (y-z)^2)
    n1 = -P8(kp) * d(xpwp, xw) * d(xz, yz) / (sq(xpwp) * B * A * sq(yz))
    # d^3N_2 = -(1/(2pi)^3) int { ... }: inside the braces
    br2 = (-P8(kp) * d(xpwp, xw) * d(yz, zw) / (sq(xpwp) * sq(yz) * sq(zw) * B)            # yours
           + P8(kp) * (d(xpwp, zw) * d(yz, xw) / (sq(xpwp) * sq(yz) * sq(zw) * (B - A))  # blue, from log[(..k+ - .. vee)/(..)]
                       - d(xpwp, xw) * d(yz, zw) * A / (sq(xpwp) * sq(yz) * sq(zw) * (B - A) * B)
                       + d(xpwp, yz) * d(zw, xw) * A / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2)
                       - d(xpwp, zw) * d(yz, xz) * B / (sq(xpwp) * sq(yz) * sq(zw) * (B - A) * A)
                       + d(xpwp, xz) * d(yz, zw) / (sq(xpwp) * sq(yz) * sq(zw) * (B - A))
                       - d(xpwp, yz) * d(zw, xz) * B / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2)))
    n2 = -br2
    # d^3N_3a: yours, log(k+/vee) = -L_v
    n3 = -(-P8(kp) * d(xpwp, yz) / (sq(xpwp) * sq(yz) * sq(zw))
           - P8(kp) * d(xpwp, yz) * d(zw, xw) * A / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2)
           + P8(kp) * d(xpwp, yz) * d(zw, xz) * B / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2))
    return UL.ph(c) * (n1 + n2 + n3) * ft.col_G1_III(c)


def g1_iv(c, kp):
    p = pts(c)
    Xp, Y = p['xp'] - p['wp'], p['y'] - p['z']
    Xp, Y = Xp / sq(Xp), Y / sq(Y)
    yours = -P16(kp) * d(Xp, Y) / sq(p['z'] - p['w'])
    blue = P8(kp) * T_contract(c, Xp, Y, 'w')
    return UL.ph(c) * (yours + blue) * ft.col_G1_IV(c)


# ---------------------------------------------------------------- Group II
def g2_i1(c, kp):
    p = pts(c)
    Y = (p['y'] - p['z']) / sq(p['y'] - p['z'])
    XpA = (p['xp'] - p['wp']) / sq(p['xp'] - p['wp'])
    XpB = (p['xp'] - p['w']) / sq(p['xp'] - p['w'])
    A_ = (P32(kp) * d(XpA, Y) / sq(p['z'] - p['w']) - P16(kp) * T_contract(c, XpA, Y, 'w')) * ft.col_G2_I1(c, 'w', 'wp')
    B_ = (-P32(kp) * d(XpB, Y) / sq(p['z'] - p['wp']) + P16(kp) * T_contract(c, XpB, Y, 'wp')) * ft.col_G2_I1(c, 'wp', 'w')
    return UL.ph(c) * (A_ + B_)


def g2_i2(c, kp):
    p = pts(c)
    x, xp, z, w, wp = p['x'], p['xp'], p['z'], p['w'], p['wp']
    A, B, Ap, Bp = sq(x - z), sq(x - w), sq(xp - z), sq(xp - wp)
    zw2, zwp2 = sq(z - w), sq(z - wp)
    b2b = (d(xp - wp, x - w) * d(z - w, z - wp) / (zw2 * zwp2) + d(xp - wp, x - w) * d(xp - z, z - w) / (2 * zw2 * Ap)
           + d(x - w, xp - wp) * d(x - z, z - wp) / (2 * A * zwp2) + d(x - z, xp - z) * d(xp - wp, x - w) / (4 * A * Ap))
    b2c = (d(xp - z, x - z) * d(z - w, z - wp) / (zw2 * zwp2) - d(xp - z, x - z) * d(xp - wp, z - w) / (2 * zw2 * Bp)
           - d(x - w, z - wp) * d(x - z, xp - z) / (2 * B * zwp2) + d(x - z, xp - z) * d(xp - wp, x - w) / (4 * B * Bp))
    pre = G**4 / (8 * PI**5) / kp
    return UL.ph(c) * pre * (-b2b / (B * Bp) + b2c / (A * Ap)) * ft.col_G2_I2(c)


# ---------------------------------------------------------------- Group III
def g3_i(c, kp):
    p = pts(c)
    Xp, Y = (p['xp'] - p['wp']) / sq(p['xp'] - p['wp']), (p['y'] - p['z']) / sq(p['y'] - p['z'])
    return UL.ph(c) * (P16(kp) * d(Xp, Y) / sq(p['z'] - p['w']) - P8(kp) * T_contract(c, Xp, Y, 'w')) * ft.col_G3_I(c)


def g3_ii(c, kp):
    p = pts(c)
    Xp, Y = (p['xp'] - p['wp']) / sq(p['xp'] - p['wp']), (p['y'] - p['z']) / sq(p['y'] - p['z'])
    return UL.ph(c) * (P16(kp) * d(Xp, Y) / sq(p['z'] - p['w']) - P8(kp) * T_contract(c, Xp, Y, 'w')) * ft.col_G3_II(c)


def _t2(c):
    p = pts(c)
    ypz, xz, xpyp, yw = p['yp'] - p['z'], p['x'] - p['z'], p['xp'] - p['yp'], p['y'] - p['w']
    return d(ypz, xz) * d(xpyp, yw) / (sq(ypz) * sq(xpyp) * sq(yw) * sq(xz))


def _ph3(c):
    return np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) / ft.PI3


def crow(coef, colf):
    def f(c, kp):
        return _ph3(c) * coef / kp * _t2(c) * colf(c)
    return f


def g3_iii2(c, kp):
    p = pts(c)
    x, z, w, xp, yp = p['x'], p['z'], p['w'], p['xp'], p['yp']
    br = (d(x - w, xp - yp) * d(z - w, yp - z) / (sq(z - w) * sq(x - w) * sq(yp - z) * sq(xp - yp))
          + d(x - z, yp - z) * d(x - w, xp - yp) / (2 * sq(x - z) * sq(x - w) * sq(yp - z) * sq(xp - yp)))
    return _ph3(c) * G**4 / (8 * PI**5) / kp * br * ft.col_G3_III2(c)


ROWS = {
    'G1-II(2a)': UL.g1_ii_2a,
    'G1-III': g1_iii,
    'G1-IV': g1_iv,
    'G1-V(1)': UL.g1_v1,
    'G1-V(2a)': UL.g1_v2a,
    'G1-V(2b fin)': UL.g1_v2b_fin,
    'G2-I(1)': g2_i1,
    'G2-I(2)': g2_i2,
    'G3-I': g3_i,
    'G3-II': g3_ii,
    'G3-III(1)': crow(1j * G**4 / (16 * PI**5), ft.col_G3_III1),
    'G3-III(2)': g3_iii2,
    'G3-V': crow(1j * G**4 / (8 * PI**5), ft.col_G3_V),
    'G3-VI': crow(1j * G**4 / (8 * PI**5), ft.col_G3_VI),
}


def run(seeds=(1, 2, 3, 4, 5)):
    out = {}
    for name, f in ROWS.items():
        rel, rel_d = [], []
        for s in seeds:
            c = SP.cfg(s)
            kp = 1.0
            v = complex(f(c, kp))
            C = complex(RC.CROW[name](c, kp))
            rel.append(abs(v - C) / abs(C))
            if name in RC.FILE and s <= 3:
                osc = SP.ROWS[name][2] if name in SP.ROWS else None
                if osc is None:
                    fl = RC.FILE[name]
                else:   # Group III Rows III, V, VI: remove the bounded e^{-i kappa (1-a)/a} terms (zero after int_y')
                    fl = (lambda F, O: (lambda c_, k_, V_: F(c_, k_, V_) - O(c_, k_, V_)))(RC.FILE[name], osc)
                dd = complex(RC.dlog(fl, c, kp))
                rel_d.append(abs(v - dd) / abs(dd))
        out[name] = dict(vs_exact=max(rel), vs_dlog_file=max(rel_d) if rel_d else None)
        print('%-13s corrected vs exact C_row: %.1e   vs d/dlog(vee) of the file: %s'
              % (name, max(rel), '%.1e' % max(rel_d) if rel_d else '-'))
    return out


if __name__ == '__main__':
    json.dump(run(), open('corrected_check.json', 'w'), indent=1)
