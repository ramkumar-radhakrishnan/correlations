"""Three-rho cancellation test for the z-independent colour structure
      C_w = f^{abc} CL^a (U(w) rho(x))^b (U(w) rho(y))^c        (antisymmetric under x <-> y).
It comes only from Group I Row II 2a (symmetric part of O_II = 2 C_w - ...) and from the K_IV family
(G1 IV, G3 I, G3 II, G2 I(1) block A; symmetric part of their colour = -C_w + ..., FS Eq. (8.2)).
Since C_w does not depend on z, its total coefficient is the z-integral of the kernels, antisymmetrized in x <-> y:
   Row II 2a:   -(i g^4/(32 pi^4 k+)) X'.F(w;x,y) x 2
   K_IV family: -(K) x int_z [kernel]_A,   kernel = p8 Ck_iv (exact)  or  p8 (-X'.Y/(2(z-w)^2)) (your list)
Cancellation requires  int_z [Ck_iv]_A = -(pi/2) X'.F(w;x,y).  With your kernel the integral diverges at z -> w."""
import json
import numpy as np
import env  # noqa: F401
from fp import sq

rng = np.random.default_rng(4)


def K(v):
    return v / (v @ v)


def ck_iv(P, z):
    xp, wp, x, y, w = P['xp'], P['wp'], P['x'], P['y'], P['w']
    Xp, Y, Zw, Xw, Xz = K(xp - wp), K(y - z), K(z - w), K(x - w), K(x - z)
    return (Xp @ Zw) * (Y @ Xz) - (Xp @ Xw) * (Y @ Xz) - (Xp @ Xw) * (Y @ Zw)


def user_k(P, z):
    """your K_IV-family kernel in the same units: -X'.Y/(2 (z-w)^2)"""
    return -0.5 * (K(P['xp'] - P['wp']) @ K(P['y'] - z)) / sq(z - P['w'])


def integrate(f, P, centers, nr=500, nt=256, rmin=0.0):
    """int d^2z f(z): partition of unity around the singular points, polar coordinates, r = s/(1-s)"""
    s, ws = np.polynomial.legendre.leggauss(nr)
    s = 0.5 * (s + 1)
    ws = 0.5 * ws
    r = s / (1 - s)
    dr = ws / (1 - s)**2
    th = (np.arange(nt) + 0.5) * 2 * np.pi / nt
    dth = 2 * np.pi / nt
    tot = 0.0
    C = [P[c] for c in centers]
    for a in C:
        for ri, dri in zip(r, dr):
            if ri < rmin:
                continue
            z = a[None, :] + ri * np.stack([np.cos(th), np.sin(th)], 1)
            vals = np.array([f(P, zz) for zz in z])
            wts = np.array([1 / sq(zz - a)**1 / sum(1 / sq(zz - b) for b in C) for zz in z])
            tot += dri * ri * dth * np.sum(vals * wts)
    return tot


def F_dot(P):
    x, y, w = P['x'], P['y'], P['w']
    F = K(x - w) * np.log(sq(w - y) / sq(x - y)) - K(y - w) * np.log(sq(w - x) / sq(x - y))
    return K(P['xp'] - P['wp']) @ F


def swap(P):
    Q = dict(P)
    Q['x'], Q['y'] = P['y'], P['x']
    return Q


if __name__ == '__main__':
    out = []
    for trial in range(3):
        P = {k: rng.normal(size=2) for k in ('xp', 'wp', 'x', 'y', 'w')}
        I1 = integrate(ck_iv, P, ['w', 'x', 'y'])
        I2 = integrate(ck_iv, swap(P), ['w', 'x', 'y'])
        IA = 0.5 * (I1 - I2)
        target = -0.5 * np.pi * F_dot(P)
        # your kernel: cut out a disc |z-w| < eps and watch the growth with log(1/eps)
        Pu = dict(P)
        grow = []
        for eps in (1e-2, 1e-3, 1e-4):
            u1 = integrate(user_k, P, ['w'], rmin=eps)
            u2 = integrate(user_k, swap(P), ['w'], rmin=eps)
            grow.append(float(0.5 * (u1 - u2)))
        rec = dict(int_Ck_iv_A=float(IA), required=float(target), rel=float(abs(IA - target) / abs(target)),
                   your_kernel_with_cutoff=grow, your_kernel_slope_per_log=float((grow[2] - grow[0]) / np.log(100.0)))
        out.append(rec)
        print(rec)
    json.dump(out, open('cw_test.json', 'w'), indent=1)
