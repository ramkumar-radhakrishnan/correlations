#!/usr/bin/env python3
"""
ref3_dilute / v3c: Row I.III, vectorised in s = p+/k+, finite Lambda (lower cut and |s-1|>Lambda), vee -> infinity
(the coefficient of ln vee is returned separately).  PV at s0 = p^2/k^2 by subtraction; the 1/|s-1| behaviour at s=1
subtracted analytically.  R_III(Lambda) = int d^2p/(2pi)^2 h(p;Lambda), Gauss-Legendre in (ln|p|, theta).
Cross-checked against the scalar version (v3_rowIII.Rdens) pointwise.
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v3_rowIII as v3
from v0_common import V0

PF, PAC = v3.P_F, v3.P_AC

def Fvec(s, p, k, l):
    """R density per ds d^2p/(2pi)^2, vectorised in s (numpy array)."""
    v0 = V0(k, l); V = v3.Vv(p, l)
    k2, p2 = k @ k, p @ p
    n = k - p
    vv0 = v0 @ v0
    a_v0V, a_v0k, a_v0p, a_v0n = v0 @ V, v0 @ k, v0 @ p, v0 @ n
    Vk, Vp, Vn = V @ k, V @ p, V @ n
    nk, np_ = n @ k, n @ p
    sq = np.sqrt(s)
    T1 = -1j * sq * (1 + s) / (s - 1) ** 2 / (s * k2 - p2)
    Y = (p2 + k2 * s) / (p2 - k2 * s)
    T2c = -Vp / (k2 * p2 * sq) * 1j * Y                    # coefficient of k_i in V@T2
    d = p[None, :] - s[:, None] * k[None, :]
    G = np.where(s > 1, (s - 1) / (n @ n), -s * (1 - s) / np.einsum('ij,ij->i', d, d))
    T3 = -2j * sq / (1 - s) ** 2 / (k2 * s - p2) * G
    v0w = s * a_v0k - a_v0p; Vw = s * Vk - Vp; nw = s * nk - np_
    v0VT3 = Vw * a_v0n / s - Vn * v0w - nw * a_v0V / (s - 1)
    v0WF = np.sqrt(2) / (2 * np.pi * sq) * (T1 * a_v0V + T2c * a_v0k + T3 * v0VT3)
    v0WA = np.sqrt(2) / (2 * np.pi * sq) * (-Vp / (k2 * p2 * sq)) * a_v0k
    return 2 * (PF * v0WF + PAC * v0WA).real / (np.sqrt(2) * vv0)

GL = {n: np.polynomial.legendre.leggauss(n) for n in (16, 32, 48, 64, 96, 128)}

def gl(f, a, b, n=48):
    x, w = GL[n]
    t = 0.5 * (b - a) * x + 0.5 * (b + a)
    return 0.5 * (b - a) * np.sum(w * f(t))

def gl_log(f, a, b, n=48):
    """int_a^b f(s) ds with s = e^u."""
    x, w = GL[n]
    ua, ub = np.log(a), np.log(b)
    u = 0.5 * (ub - ua) * x + 0.5 * (ub + ua)
    s = np.exp(u)
    return 0.5 * (ub - ua) * np.sum(w * f(s) * s)

def h(p, k, l, Lam, a=0.5):
    k2 = k @ k; s0 = (p @ p) / k2
    if abs(s0 - 1) < 1e-6:
        p = p * (1 + 3e-6); s0 = (p @ p) / k2
    F = lambda s: Fvec(np.atleast_1d(s), p, k, l)
    e = 1e-7
    c = 0.5 * (F(1 + e)[0] + F(1 - e)[0]) * e
    Ainf = F(1e9)[0] * 1e9
    # residue at s0: g(s) = F(s)(s-s0)
    if Lam < s0 and abs(s0 - 1) > 1e-9:
        es = 1e-7 * s0
        r0 = 0.5 * (F(s0 + es)[0] * es - F(s0 - es)[0] * es)
    else:
        r0 = 0.0
    def pv_piece(f, x0, x1, logmap=False):
        """int_{x0}^{x1} f, with PV at s0 if inside (subtracting r0/(s-s0))."""
        inside = x0 < s0 < x1
        if inside:
            fs = lambda s: f(s) - r0 / (s - s0)
            # split at s0 for accuracy
            val = 0.0
            for (y0, y1) in ((x0, s0), (s0, x1)):
                val += (gl_log(fs, y0, y1, 64) if logmap else gl(fs, y0, y1, 64))
            val += r0 * np.log((x1 - s0) / (s0 - x0))
            return val
        return gl_log(f, x0, x1, 64) if logmap else gl(f, x0, x1, 64)
    tot = 0.0
    # [Lam, 1-a]: log map, with extra splits around s0
    edges = sorted(set([Lam, 1 - a] + [x for x in (s0 / 4, s0 / 2, 2 * s0, 4 * s0) if Lam < x < 1 - a]))
    for x0, x1 in zip(edges[:-1], edges[1:]):
        tot += pv_piece(F, x0, x1, logmap=True)
    # around s = 1 (exclusion |s-1|>Lam): subtract c/|s-1|
    fm = lambda s: F(s) - c / np.abs(s - 1)
    for (x0, x1) in ((1 - a, 1 - 1e-12), (1 + 1e-12, 1 + a)):
        # log map in |s-1| for accuracy near 1
        if x0 < 1:
            g = lambda t: fm(1 - t)
            val = 0.0
            if 1 - x1 < s0 and False: pass
            sub_edges = sorted(set([1e-12, a] + [x for x in (abs(1 - s0) / 2, abs(1 - s0), 2 * abs(1 - s0)) if 1e-12 < x < a and s0 < 1]))
            for y0, y1 in zip(sub_edges[:-1], sub_edges[1:]):
                if s0 < 1 and y0 < 1 - s0 < y1:
                    t0 = 1 - s0
                    gs = lambda t: g(t) - (-r0) / (t - t0)
                    val += gl_log(gs, y0, t0, 48) + gl_log(gs, t0, y1, 48) + (-r0) * np.log((y1 - t0) / (t0 - y0))
                else:
                    val += gl_log(g, y0, y1, 48)
            tot += val
        else:
            g = lambda t: fm(1 + t)
            val = 0.0
            sub_edges = sorted(set([1e-12, a] + [x for x in (abs(s0 - 1) / 2, abs(s0 - 1), 2 * abs(s0 - 1)) if 1e-12 < x < a and s0 > 1]))
            for y0, y1 in zip(sub_edges[:-1], sub_edges[1:]):
                if s0 > 1 and y0 < s0 - 1 < y1:
                    t0 = s0 - 1
                    gs = lambda t: g(t) - r0 / (t - t0)
                    val += gl_log(gs, y0, t0, 48) + gl_log(gs, t0, y1, 48) + r0 * np.log((y1 - t0) / (t0 - y0))
                else:
                    val += gl_log(g, y0, y1, 48)
            tot += val
    tot += 2 * c * np.log(a / Lam)
    # [1+a, S] and tail
    S = max(20.0, 8 * s0)
    edges = sorted(set([1 + a, S] + [x for x in (s0 / 2, 2 * s0) if 1 + a < x < S]))
    for x0, x1 in zip(edges[:-1], edges[1:]):
        tot += pv_piece(F, x0, x1, logmap=True)
    tail = lambda s: F(s) - Ainf / s
    tot += gl_log(tail, S, S * 1e8, 96) - Ainf * np.log(S)
    return tot, Ainf

def R_III(k, l, Lam, nr=96, nth=48):
    """int d^2p/(2pi)^2 h(p); returns (value at vee=1, coefficient of ln vee)."""
    lams = np.sqrt(Lam) * np.linalg.norm(k)
    sc = sorted(set([lams, np.linalg.norm(k), np.linalg.norm(l), np.linalg.norm(k + l)]))
    edges = [min(lams, sc[0]) * 1e-3] + sc + [sc[-1] * 1e3]
    th, wth = np.polynomial.legendre.leggauss(nth)
    th = np.pi * (th + 1); wth = np.pi * wth
    xr, wr = np.polynomial.legendre.leggauss(nr // 4)
    tot = np.zeros(2)
    for r0, r1 in zip(edges[:-1], edges[1:]):
        u0, u1 = np.log(r0), np.log(r1)
        for x, w in zip(xr, wr):
            u = 0.5 * (u1 - u0) * x + 0.5 * (u1 + u0); r = np.exp(u)
            for t, wt in zip(th, wth):
                p = np.array([r * np.cos(t), r * np.sin(t)])
                hv = h(p, k, l, Lam)
                tot += np.array(hv) * w * 0.5 * (u1 - u0) * wt * r * r
    return tot / (2 * np.pi) ** 2

if __name__ == '__main__':
    k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7])
    # pointwise check against scalar version
    for p in (np.array([0.4, 0.5]), np.array([0.05, 0.03]), np.array([2.0, -1.0])):
        for s in (0.01, 0.3, 0.97, 1.2, 5.0):
            assert abs(Fvec(np.array([s]), p, k, l)[0] - v3.Rdens(s, p, k, l)) < 1e-10 * (1 + abs(v3.Rdens(s, p, k, l)))
    print("vectorised integrand agrees with v3_rowIII.Rdens")
    import importlib.util, os
    spec = importlib.util.spec_from_file_location('probe3', '/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/ref3/dilute/probe3.py')
    for p in (np.array([0.4, 0.5]), np.array([0.03, 0.04]), np.array([0.6, 0.8]) * 1.001, np.array([2.0, -1.0])):
        print(p, h(p, k, l, 1e-6))
