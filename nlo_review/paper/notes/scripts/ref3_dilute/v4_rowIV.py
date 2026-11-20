#!/usr/bin/env python3
"""
ref3_dilute / v4: Row I.IV (A_AB2) in the dilute limit, point-like quark, momentum space (k+ = 1, g = 1).
A_AB2 -> 2 int ds/(2pi) int d^2p/(2pi)^2 [ u^{ad}(l) A^dag_j(p) B2^{bd}_{ji}(1,k+l; s,p)
                                         + A^dag_j(p) u^{bc}(l) B2^{ca}_{ji}(1,k; s,p+l)
                                         - A^dag_j(p) B2^{ba}_{ji}[u rho](1,k; s,p) ],
A^dag^b_j(p) = -sqrt2 p_j t^b/(sqrt s p^2) (point source); B2 = 1-rho (beta) + {rho,rho} part (1/4 A A {t,t}).
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
from v0_common import Colour, V0, beta, Amom

C = Colour(3)
f, t = C.f, C.t
AC = np.einsum('aij,bjk->abik', t, t) + np.einsum('bij,ajk->abik', t, t)     # {t^a,t^b}[a,b]
T = {
 'w1': np.einsum('axd,bij,dbc,cjk->axik', f, t, f, t),                     # u^{ad} t^b f^{dbc} t^c
 'w2': np.einsum('axd,bij,bdjk->axik', f, t, AC),                          # u^{ad} t^b {t^b,t^d}
 'x1': np.einsum('bij,bxc,acy,yjk->axik', t, f, f, t),                     # t^b f^{bxc} f^{acy} t^y
 'x2': np.einsum('bij,bxc,cajk->axik', t, f, AC),                          # t^b f^{bxc} {t^c,t^a}
 'd1': -np.einsum('bij,aby,yxd,djk->axik', t, f, f, t),                    # -t^b f^{aby} f^{yxd} t^d
 'd2': -(np.einsum('bij,bxd,dajk->axik', t, f, AC) + np.einsum('bij,axd,bdjk->axik', t, f, AC)),
}
P = {s: C.proj(T[s]) for s in T}

def Adag(s, p):
    return -np.sqrt(2) * p / (np.sqrt(s) * (p @ p))

def Rdens(s, p, k, l, keys=None):
    v0 = V0(k, l)
    A = Adag(s, p)
    W = {}
    W['w1'] = A @ beta(1.0, k + l, s, p)
    W['w2'] = A @ (0.25 * np.outer(Amom(s, p), Amom(1.0, k + l)))
    W['x1'] = A @ beta(1.0, k, s, p + l)
    W['x2'] = A @ (0.25 * np.outer(Amom(s, p + l), Amom(1.0, k)))
    W['d1'] = A @ beta(1.0, k, s, p)
    W['d2'] = A @ (0.25 * np.outer(Amom(s, p), Amom(1.0, k)))
    keys = keys or W.keys()
    pre = 2 / (2 * np.pi)
    return sum(2 * (P[c] * (v0 @ (pre * W[c]))).real for c in keys) / (np.sqrt(2) * (v0 @ v0))

if __name__ == '__main__':
    print("projections:", {c: np.round(P[c], 6) for c in P})
    k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7])
    for s in (0.3, 3.0):
        for r in (1e-2, 1e-3, 1e-4):
            vals = [Rdens(s, r * np.array([np.cos(a), np.sin(a)]), k, l) * r * r for a in np.linspace(0, 2 * np.pi, 7)[:-1]]
            print(f"s={s} |p|={r}: p^2 * R_IV = " + " ".join(f"{v:.5f}" for v in vals))
        for r in (1e2, 1e3):
            vals = [Rdens(s, r * np.array([np.cos(a), np.sin(a)]), k, l) * r * r for a in np.linspace(0, 2 * np.pi, 7)[:-1]]
            print(f"s={s} |p|={r}: p^2 * R_IV = " + " ".join(f"{v:.5f}" for v in vals))
