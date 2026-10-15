"""Transcriptions of the integrated lines (after p+ / w' integration) of Group_II.tex and Group_III.tex."""
import numpy as np
from fp import F, NC, G, PI, I2, K, sq
import user23 as us

pi3 = us.pi3


def logs(c, k, Lam, vee, X, W):
    A, B = sq(c.pt[X] - c.pt['z']), sq(c.pt[X] - c.pt[W])
    return np.log(((vee - k) * A + k * B) / (Lam * A + k * B))


def br4(c, W, X='x'):
    """[d/(2(z-w)^2) + (x-z)^j(z-w)^i/((z-w)^2(x-z)^2) - (x-z)^j(x-w)^i/((x-w)^2(x-z)^2) - (x-w)^i(z-w)^j/((z-w)^2(x-w)^2)] [i,j]"""
    x, z, w = c.pt[X], c.pt['z'], c.pt[W]
    xz, xw, zw = x - z, x - w, z - w
    return (I2 / (2 * sq(zw)) + np.outer(zw, xz) / (sq(zw) * sq(xz)) - np.outer(xw, xz) / (sq(xw) * sq(xz))
            - np.outer(xw, zw) / (sq(zw) * sq(xw)))


def brlog(c, W, X='x'):
    """[(x-w)^i(z-w)^j/((z-w)^2(x-w)^2) + (x-z)^j(x-w)^i/(2(x-z)^2(x-w)^2)]"""
    x, z, w = c.pt[X], c.pt['z'], c.pt[W]
    xz, xw, zw = x - z, x - w, z - w
    return np.outer(xw, zw) / (sq(zw) * sq(xw)) + np.outer(xw, xz) / (2 * sq(xz) * sq(xw))


def kij(c, a, b, M):
    return np.einsum('i,j,ij->', a, b, M)


# ------------------------------------------------------------------ Group II Row I part 1, line 50 (block A and block B)
def g2_1_line50(c, k, Lam, vee):
    import fixed23 as fx
    Lv = np.log(vee / Lam) + np.log(1 - k / vee)
    lk = np.log(vee / (Lam + k))
    # block A: kernel with w (amp), K(x'-w'), colour two(w') x one(w)
    a1 = K(c.pt['xp'] - c.pt['wp']); b1 = K(c.pt['y'] - c.pt['z'])
    kerA = (-1j * G**4 / (16 * PI**5) / k * kij(c, a1, b1, br4(c, 'w')) * logs(c, k, Lam, vee, 'x', 'w')
            + 1j * G**4 / (32 * PI**5) / k * lk * (a1 @ b1) / sq(c.pt['z'] - c.pt['w'])
            - 1j * G**4 / (16 * PI**5) / k * Lv * kij(c, a1, b1, brlog(c, 'w')))
    a2 = K(c.pt['xp'] - c.pt['w'])
    kerB = (1j * G**4 / (16 * PI**5) / k * kij(c, a2, b1, br4(c, 'wp')) * logs(c, k, Lam, vee, 'x', 'wp')
            - 1j * G**4 / (32 * PI**5) / k * lk * (a2 @ b1) / sq(c.pt['z'] - c.pt['wp'])
            + 1j * G**4 / (16 * PI**5) / k * Lv * kij(c, a2, b1, brlog(c, 'wp')))
    colA = np.sum(fx.two_fix(c, 'wp', 'xp', 'y') * fx.Rone_fix(c, 'w', 'z', 'x'))
    colB = np.sum(fx.Rone_fix(c, 'wp', 'z', 'x') * fx.two_fix(c, 'w', 'xp', 'y'))
    return kerA * colA / pi3, kerB * colB / pi3


# ------------------------------------------------------------------ Group II Row I part 2: 2a .. 2f (lines 66-82)
def g2_2_integrated(c, k, Lam, vee, d=2):
    """returns dict of kernels (colour stripped, colour58 common): 2a..2f"""
    gr = us.g2_2_groups62(c, k, 0.5, d)
    A, B = sq(c.pt['x'] - c.pt['z']), sq(c.pt['x'] - c.pt['w'])
    Ap, Bp = sq(c.pt['xp'] - c.pt['z']), sq(c.pt['xp'] - c.pt['wp'])
    zw2, zwp2 = sq(c.pt['z'] - c.pt['w']), sq(c.pt['z'] - c.pt['wp'])
    L = logs(c, k, Lam, vee, 'x', 'w')
    Lp = logs(c, k, Lam, vee, 'xp', 'wp')
    Lv = np.log(vee / Lam) + np.log(1 - k / vee)
    lk = np.log(vee / (Lam + k))
    pre = G**4 / (8 * PI**5) / k / pi3
    out = {}
    out['2a'] = pre / (Bp * A - B * Ap) * gr['pk'] * (np.log(((vee - k) * A + k * B) / (k * Bp + (vee - k) * Ap))
                                                      - np.log((Lam * A + k * B) / (k * Bp + Lam * Ap)))
    out['2b'] = (pre * Lv / (Bp * A * B - B**2 * Ap) * gr['pp'] * (A - B * Ap / Bp)
                 - pre * A / (Bp * A * B - B**2 * Ap) * gr['pp'] * L
                 + pre * Ap / (Bp**2 * A - B * Ap * Bp) * gr['pp'] * Lp)
    out['2c'] = (pre * Bp / Ap / (Bp * A - B * Ap) * gr['kk'] * Lp - pre * B / A / (Bp * A - B * Ap) * gr['kk'] * L)
    out['2d'] = (-pre * lk / (A - B) / (Bp - Ap) * gr['p'] + pre * A / (A - B) / (Bp * A - B * Ap) * gr['p'] * L
                 + pre * Ap / (Bp - Ap) / (Bp * A - B * Ap) * gr['p'] * Lp)
    out['2e'] = (pre * lk / (A - B) / (Bp - Ap) * gr['k'] - pre * B / (A - B) / (Bp * A - B * Ap) * gr['k'] * L
                 - pre * Bp / (Bp - Ap) / (Bp * A - B * Ap) * gr['k'] * Lp)
    pf = G**4 * d / (32 * PI**5) / pi3
    out['2f'] = (pf * lk / k * (Bp * B - A * Ap) / (zw2 * zwp2 * (A - B) * (Ap - Bp))
                 - pf / k * A * B / (zw2 * zwp2 * (A - B)) * (Ap - Bp) / (Bp * A - B * Ap) * L
                 + pf * (1 / vee - 1 / (k + Lam)) / (zw2 * zwp2)
                 + pf / k * Bp * Ap / (zw2 * zwp2) * (A - B) / ((Ap - Bp) * (Bp * A - B * Ap)) * Lp)
    return out


# ------------------------------------------------------------------ Group III Rows I, II after p+ (lines 48, 59)
def g3_integrated(c, k, Lam, vee):
    Lv = np.log(vee / Lam) + np.log(1 - k / vee)
    lk = np.log(vee / (k + Lam))
    a1 = K(c.pt['xp'] - c.pt['wp']); b1 = K(c.pt['y'] - c.pt['z'])
    ker = (-1j * G**4 / (8 * PI**5) / k * kij(c, a1, b1, br4(c, 'w')) * logs(c, k, Lam, vee, 'x', 'w')
           + 1j * G**4 / (16 * PI**5) / k * lk * (a1 @ b1) / sq(c.pt['z'] - c.pt['w'])
           - 1j * G**4 / (8 * PI**5) / k * Lv * kij(c, a1, b1, brlog(c, 'w')))
    return ker / pi3
