"""The terms of input/Finite_terms_1_1.tex (your log(vee/k+) list), transcribed term by term.

Every function returns the coefficient of L_v = log(vee/k+) (log(k+/vee) = -L_v) of one row, at the transverse
points of a configuration c (classical charges, U = U^T as in your files), with 1/(2 pi)^3, the prefactor, the
phase and the colour, without c.c. -- the same normalization as ../finite/ft.py and ../finite_v2/split.py.
Rows marked NA in the file return 0."""
import numpy as np
import env  # noqa: F401
from fp import G, PI, sq, NC
import ft
from ft import vec, Kv, ij, PI3, AB
import simp_exact as SE

P8 = lambda kp: 1j * G**4 / (8 * PI**5) / kp
P16 = lambda kp: 1j * G**4 / (16 * PI**5) / kp
P32 = lambda kp: 1j * G**4 / (32 * PI**5) / kp


def ph(c, a='w', b='wp'):
    return ft.phase(c, a, b) / PI3


# ---------------------------------------------------------------- Group I
def g1_ii_2a(c, kp):
    """-(1/(2pi)^3) i g^4 f^{hbc}/(32 pi^4) (1/k+) L_v [K(x'-w')^i F^i] CL^a U^{ah}(w){rho^b(x),rho^c(y)}
       +(1/(2pi)^3) i g^4 f^{abc}/(32 pi^4) (1/k+) L_v [K(x'-w')^i F^i] CL^a U^{bd}(x)U^{ce}(y){rho^d(x),rho^e(y)}
    (the int_z of the file is spurious: nothing depends on z)"""
    x, y, w = c.pt['x'], c.pt['y'], c.pt['w']
    F = Kv(c, 'x', 'w') * np.log(sq(w - y) / sq(x - y)) - Kv(c, 'y', 'w') * np.log(sq(w - x) / sq(x - y))
    k = Kv(c, 'xp', 'wp') @ F
    ca, cb = ft.col_G1_II_2a(c)
    pre = ph(c) * 1j * G**4 / (32 * PI**4) / kp * k
    return -pre * ca + pre * cb


def g1_iii(c, kp):
    A, B = AB(c)
    xz, xw, zw = vec(c, 'x', 'z'), vec(c, 'x', 'w'), vec(c, 'z', 'w')
    Xp, Yz, Zw, Xw = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'z', 'w'), Kv(c, 'x', 'w')
    XY = Xp @ Yz
    n1 = -P16(kp) * XY * (A + B) / (A - B)**2                       # d^3N_1 (coefficient of L_v)
    n2 = +P8(kp) * (Xp @ Xw) * (Yz @ Zw)                           # d^3N_2 = -(...)(-i g^4/8pi^5)
    # d^3N_3a: three terms with log(k+/vee) = -L_v
    n3 = (P8(kp) * XY / sq(zw)
          + P8(kp) * XY * (zw @ xw) * A / (sq(zw) * (B - A)**2)
          - P8(kp) * XY * (zw @ xz) * B / (sq(zw) * (B - A)**2))
    return ph(c) * (n1 + n2 + n3) * ft.col_G1_III(c)


def g1_iv(c, kp):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    return ph(c) * (-P16(kp)) * (Xp @ Yz) / sq(vec(c, 'z', 'w')) * ft.col_G1_IV(c)


def g1_v1(c, kp):
    kk = ft.G1_V1_kernels(c)
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    # -(1/(2pi)^3) i g^4/(32 pi^5 k+) log(k+/vee) [...]  =  +(1/(2pi)^3) i g^4/(32 pi^5 k+) L_v [...]
    return P32(kp) / PI3 * (ez * kk['b'] + ex * kk['q']) * ft.col_G1_V1(c)


def _v2_user(c):
    """the two brackets of 2a as written (contracted)"""
    x, y, z = c.pt['x'], c.pt['y'], c.pt['z']
    Xp, Zx = Kv(c, 'xp', 'wp'), Kv(c, 'z', 'x')
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    bz = ((Xp @ (y - z)) * (Zx @ (x - z)) / sq(x - z) + (Zx @ (y - x)) * (Xp @ (y - z)) / (2 * sq(y - x))) / sq(y - z)
    bx = ((Xp @ (y - x)) * (Zx @ (x - z)) / sq(x - z) - (Xp @ (y - x)) * (Zx @ (y - z)) / (2 * sq(y - z))) / sq(y - x)
    return ez * bz + ex * bx


def g1_v2a(c, kp):
    # -(1/(2pi)^3) g^4/(16 pi^5 k+) log(k+/vee) [...] = +(...) L_v
    return G**4 / (16 * PI**5) / kp / PI3 * _v2_user(c) * SE.col_V2a(c)


def g1_v2b_fin(c, kp):
    x, y, z = c.pt['x'], c.pt['y'], c.pt['z']
    Xp, Zx = Kv(c, 'xp', 'wp'), Kv(c, 'z', 'x')
    ez = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['z']))
    ex = np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['x']))
    yx, yz, xz = y - x, y - z, x - z
    xpwp = c.pt['xp'] - c.pt['wp']
    zx = z - x
    t_x = ((yx @ xpwp) * (xz @ zx) / (sq(xz) * sq(xpwp) * sq(zx) * sq(yx))
           - (yx @ xpwp) * (yz @ zx) / (2 * sq(yz) * sq(xpwp) * sq(zx) * sq(yx)))
    t_z = ((yz @ xpwp) * (xz @ zx) / (sq(xz) * sq(xpwp) * sq(zx) * sq(yz))
           + (yx @ zx) * (yz @ xpwp) / (2 * sq(yx) * sq(xpwp) * sq(zx) * sq(yz)))
    t3 = (xpwp @ yz) / (sq(xpwp) * sq(yz)) * (yx @ yz) / (sq(yx) * sq(xz))
    pre = G**4 / (16 * PI**5) / kp / PI3
    return pre * (ex * t_x + ez * t_z) * SE.col_V2b(c, 'x') + 2 * pre * ez * t3 * SE.col_V2b(c, 'z')


# ---------------------------------------------------------------- Group II
def g2_i1(c, kp):
    Xp, Yz, Xpw = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z'), Kv(c, 'xp', 'w')
    A_ = P32(kp) * (Xp @ Yz) / sq(vec(c, 'z', 'w')) * ft.col_G2_I1(c, 'w', 'wp')
    B_ = -P32(kp) * (Xpw @ Yz) / sq(vec(c, 'z', 'wp')) * ft.col_G2_I1(c, 'wp', 'w')
    return ph(c) * (A_ + B_)


def g2_i2(c, kp, d=2):
    """2d + 2e + 2f as written"""
    x, xp, z, w, wp = (c.pt[n] for n in ('x', 'xp', 'z', 'w', 'wp'))
    xz, xw, zw, xpz, xpwp, zwp = x - z, x - w, z - w, xp - z, xp - wp, z - wp
    A, B, Ap, Bp, zw2, zwp2 = sq(xz), sq(xw), sq(xpz), sq(xpwp), sq(zw), sq(zwp)
    zwp_ = z - wp
    br_d = ((xw @ zw) * (Ap - Bp) / (2 * zw2 * zwp2) + (xz @ xw) * (Ap - Bp) / (4 * A * zwp2)
            + (xpwp @ zwp_) * (A - B) / (2 * zw2 * zwp2) + (xpz @ xpwp) * (A - B) / (4 * zw2 * Ap))
    br_e = ((xz @ zw) * (Ap - Bp) / (2 * zw2 * zwp2) - (xz @ xw) * (Ap - Bp) / (4 * B * zwp2)
            + (xpz @ zwp_) * (A - B) / (2 * zw2 * zwp2) - (xpz @ xpwp) * (A - B) / (4 * zw2 * Bp))
    pre = G**4 / (8 * PI**5) / kp
    t2d = -pre / ((A - B) * (Bp - Ap)) * br_d
    t2e = +pre / ((A - B) * (Bp - Ap)) * br_e
    t2f = G**4 * d / (32 * PI**5) / kp * (Bp * B - A * Ap) / (zw2 * zwp2 * (A - B) * (Ap - Bp))
    return ph(c) * (t2d + t2e + t2f) * ft.col_G2_I2(c)


# ---------------------------------------------------------------- Group III
def g3_i(c, kp):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    return ph(c) * P16(kp) * (Xp @ Yz) / sq(vec(c, 'z', 'w')) * ft.col_G3_I(c)


def g3_ii(c, kp):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    return ph(c) * P16(kp) * (Xp @ Yz) / sq(vec(c, 'z', 'w')) * ft.col_G3_II(c)


ZERO = lambda c, kp: 0.0
ROWS = {
    'G1-II(2a)': g1_ii_2a,
    'G1-III': g1_iii,
    'G1-IV': g1_iv,
    'G1-V(1)': g1_v1,
    'G1-V(2a)': g1_v2a,
    'G1-V(2b fin)': g1_v2b_fin,
    'G2-I(1)': g2_i1,
    'G2-I(2)': g2_i2,
    'G3-I': g3_i,
    'G3-II': g3_ii,
    'G3-III(1)': ZERO,
    'G3-III(2)': ZERO,
    'G3-V': ZERO,
    'G3-VI': ZERO,
}

# coefficients of the two dim-reg rows (analytic; units of LO, with c.c.):
#   G1 Row II part 1: (alpha_s N_c/4pi) c_II L_v LO[E_X + E_X'],  c_II = 1 (your list), 2 (the file: log(vee/k+) + log(vee/k+ - 1))
#   G1 Row V 2b UV:   -(alpha_s N_c/pi) L_v LO[E_X + E_X']            (your list = the file)
#   G2 Row IV:        NA in your list;  the file: -(alpha_s/2pi) 2 N_c L_v LO[W], W = -E_Delta, + the small-xi end
