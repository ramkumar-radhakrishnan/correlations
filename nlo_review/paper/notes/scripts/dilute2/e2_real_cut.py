#!/usr/bin/env python3
"""dilute2 / e2: real emission (Groups II+III, dilute limit) at fixed s = p+/k+, integrated over p in d = 2 - eps.

Fixed q_a (projectile momentum), target distribution phi_T(q) = q^2 exp(-q^2/Q^2), q_b = k + p - q_a.
In units of  LO * (alpha_s N_c/pi) * ds/s :
   j(s) = (1/pi) int d^{2-eps}p  F_s(p),   F_s = k^2 A_gluons/(qa^2 qb^2) w(p) - w(p)/p^2,   w = phi_T(q_b)/phi_T(k-q_a)
        = -(2/eps) c(s) + j_fin(s)                       (MS-bar, mu),
c(s) = angular average of the collinear residue, lim (p-p_c)^2 F_s, p_c = s k (collinear point, x p_c = (1-x) k).
Numerics: Gaussian partition of unity around p_c and around p = 0; DR for the subtraction term:
   (1/pi) int d^{2-eps}p exp(-(p-p_c)^2/rho^2)/(p-p_c)^2 = -2/eps + ln(rho^2/mu^2) - gamma_E.
Claim to be tested:  c(s) -> 1,  j_fin(s) -> ln(s^2 k^2/mu^2)  for s -> 0   (non-uniformity of the Lambda -> 0 limit),
and  I(L) = int_{e^-L}^inf ds/s j_fin(s) = -L^2 + L ln(k^2/mu^2) + const + O(e^-L).
usage: e2_real_cut.py [n]   (n = quadrature multiplier, default 1)"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import d1_qmrk_integrand as d1
EG = 0.5772156649015329

def dot(a, b): return a[..., 0] * b[..., 0] + a[..., 1] * b[..., 1]

def Agl(qa, k1, k2, x, D=4.0):
    """vectorised BSV A_gluons; qa, k1 fixed 2-vectors, k2 array (...,2); qb = k1 + k2 - qa."""
    qa = np.broadcast_to(qa, k2.shape); k1 = np.broadcast_to(k1, k2.shape)
    qb = k1 + k2 - qa; De = k1 + k2
    La = (1 - x) * k1 - x * k2
    sh = dot(La, La) / (x * (1 - x)); De2 = dot(De, De); Sig = sh + De2
    qa2, qb2, k12, k22 = dot(qa, qa), dot(qb, qb), dot(k1, k1), dot(k2, k2)
    th = 2 * dot(qa, k1) - qa2 - k12 / x
    uh = 2 * dot(qa, k2) - qa2 - k22 / (1 - x)
    E = (2 * (2 * x - 1) * qa2 + 4 * dot(La, qa) + (1 - 2 * x) / (x * (1 - x)) * dot(La, La)
         - 2 * x * (1 - x) * ((2 * x - 1) * De2 + 2 * dot(La, De)) * qa2 / (x * (1 - x) * De2 + dot(La, La)))
    r = (1 - x) / x; rr = x / (1 - x)
    A = qa2 * qb2 * (-1 / (th * uh) + qa2 * qb2 / (4 * th * uh * k12 * k22)
                     - 0.25 * (r / (k22 * th) + rr / (k12 * uh)) + 1 / (4 * k12 * k22)
                     + 1 / Sig * (-1 / sh * (2 + (1 / th - 1 / uh) * (r * k12 - rr * k22))
                                  + 0.25 * (Sig / sh + 1) * (r / k22 + rr / k12)
                                  - qb2 / (4 * sh) * (1 / ((1 - x) * th) + 1 / (x * uh))
                                  - qa2 / (4 * sh) * ((1 + rr * k22 / k12) / th + (1 + r * k12 / k22) / uh)))
    d1q = dot(k1 - qa, k1 - qa); d2q = dot(k2 - qa, k2 - qa)
    A = A + (D - 2) / 4 * ((d1q * d2q - k12 * k22) / (th * uh)) ** 2
    A = A - (D - 2) / 4 * 0.25 * ((d2q - rr * k22) / uh + E / sh) * ((d1q - r * k12) / th - E / sh)
    return A, qa2, qb2

class Real:
    def __init__(self, k, qa, Q=1.0):
        self.k, self.qa, self.Q = k, qa, Q
        q0 = k - qa; self.phi0 = (q0 @ q0) * np.exp(-(q0 @ q0) / Q ** 2)
    def F(self, p, s):
        x = 1 / (1 + s)
        A, qa2, qb2 = Agl(self.qa, self.k, p, x)
        w = qb2 * np.exp(-qb2 / self.Q ** 2) / self.phi0
        return (self.k @ self.k) * A / (qa2 * qb2) * w - w / dot(p, p)
    def cres(self, s, n=64):
        """angular average of (p-pc)^2 F at p -> pc (Richardson in the radius)."""
        pc = s * self.k; th = 2 * np.pi * (np.arange(n) + 0.5) / n
        out = []
        for d in (1e-4, 2e-4):
            r = d * min(s, 1.0)
            p = pc + r * np.stack([np.cos(th), np.sin(th)], -1)
            out.append(np.mean(r * r * self.F(p, s)))
        return 2 * out[0] - out[1]

def gl_nodes(edges, n, log=True):
    x, w = np.polynomial.legendre.leggauss(n); X, W = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if log:
            u = 0.5 * (np.log(b) - np.log(a)) * x + 0.5 * (np.log(b) + np.log(a)); X.append(np.exp(u)); W.append(0.5 * (np.log(b) - np.log(a)) * w * np.exp(u))
        else:
            X.append(0.5 * (b - a) * x + 0.5 * (b + a)); W.append(0.5 * (b - a) * w)
    return np.concatenate(X), np.concatenate(W)

def ang_nodes(specials, n):
    brk = {0.0, 2 * np.pi}
    for x in specials:
        x = x % (2 * np.pi); brk.add(x)
        for j in range(0, 22):
            d = 1.2 * 2.0 ** (-j)
            brk.add((x + d) % (2 * np.pi)); brk.add((x - d) % (2 * np.pi))
    return gl_nodes(sorted(brk), n, log=False)

def polar(center, fun, redges, specials, nr, nth):
    r, wr = gl_nodes(redges, nr); th, wt = ang_nodes(specials, nth)
    R, TH = np.meshgrid(r, th, indexing='ij')
    p = center + np.stack([R * np.cos(TH), R * np.sin(TH)], -1)
    return np.einsum('i,j,ij->', wr * r, wt, fun(p))

def j_fin(obj, s, mu=1.0, m=1):
    k, qa = obj.k, obj.qa; pc = s * k; npc = np.linalg.norm(pc)
    rho = 0.4 * min(npc, 1.0); rho0 = 0.4 * min(npc, 1.0)
    c = obj.cres(s)
    chic = lambda p: np.exp(-dot(p - pc, p - pc) / rho ** 2)
    chi0 = lambda p: np.exp(-dot(p, p) / rho0 ** 2)
    ak = np.arctan2(k[1], k[0]); qq = qa - k; aq = np.arctan2(qq[1], qq[0])
    nr, nth = 10 * m, 6 * m
    # (a) around pc
    fa = lambda p: chic(p) * (obj.F(p, s) - c / dot(p - pc, p - pc))
    Ia = polar(pc, fa, list(rho * np.logspace(-9, np.log10(6.5), 22)), [ak, ak + np.pi], nr, nth)
    # (b) around 0
    fb = lambda p: (1 - chic(p)) * chi0(p) * obj.F(p, s)
    Ib = polar(np.zeros(2), fb, list(rho0 * np.logspace(-9, np.log10(6.5), 22)), [ak, ak + np.pi], nr, nth)
    # (c) rest, polar around 0
    fc = lambda p: (1 - chic(p)) * (1 - chi0(p)) * obj.F(p, s)
    scales = sorted([npc, 1.0, np.linalg.norm(qq), np.linalg.norm(qa)])
    ed = {rho0 * 1e-3, 1e4 * max(scales)}
    for sc in scales:
        for j in range(0, 14):
            ed.add(sc * (1 + 2.0 ** (-j))); ed.add(sc * (1 - 2.0 ** (-j - 1)))
        ed.add(sc * 4); ed.add(sc * 16); ed.add(sc / 4); ed.add(sc / 16)
    ed = sorted(e for e in ed if e > 0)
    Ic = polar(np.zeros(2), fc, ed, [ak, ak + np.pi, aq], nr, nth)
    return (Ia + Ib + Ic) / np.pi + c * (np.log(rho ** 2 / mu ** 2) - EG), c

if __name__ == '__main__':
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    k = np.array([1.0, 0.0]); qa = np.array([1.3, 0.7]); obj = Real(k, qa)
    # vectorised A_gluons against d1
    rng = np.random.default_rng(1); mx = 0
    for _ in range(20):
        p = rng.normal(size=2); x = rng.uniform(0.05, 0.95)
        a = Agl(qa, k, p[None, :], x)[0][0]; b = d1.A_gluons(qa, k + p - qa, k, p, x); mx = max(mx, abs(a / b - 1))
    print(f"vectorised A_gluons vs d1: max rel dev {mx:.1e}")
    print("#  ln s        c(s)        j_fin(s)      j_fin - c*ln(s^2)   [mu=|k|=1]")
    us = np.linspace(-20, 8, 57); J = []; C = []
    for u in us:
        s = np.exp(u); jf, c = j_fin(obj, s, m=m); J.append(jf); C.append(c)
        print(f"{u:8.3f}  {c:.9f}  {jf: .9f}  {jf - c * 2 * u: .9f}", flush=True)
    J = np.array(J); C = np.array(C)
    np.savetxt(f'e2_table_m{m}.dat', np.stack([us, C, J], 1))
