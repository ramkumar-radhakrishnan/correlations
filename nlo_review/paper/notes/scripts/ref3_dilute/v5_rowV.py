#!/usr/bin/env python3
"""
ref3_dilute / v5: Row I.V (A_CB2) in the dilute limit, point-like quark projectile, momentum space.

A_CB2 -> 2 int dp+/(2pi) int_x,z C1^dag [ (u(x)+u(z)) B2 - delta B2 ]   (O(u); C1 real)
       = 2/(2pi) int_0^1 dxi int d^2kappa/(2pi)^2 ( -i sqrt(xi(1-xi))/sqrt2 ) kappa^m/kappa^2 Br_{m;jki} X_{jk} * colour
three cases: kick on the (1-xi) gluon (z), on the xi gluon (x), rotated source (delta); 1-rho and {rho,rho} parts of B2.
UV (kappa -> inf): only the triple-gluon term of B2 contributes; its contraction with Br, averaged over directions,
is d-independent and equals  W_UV = (2/(2pi)) (1/sqrt2) Phat(xi) K_i/K^2 / kappa^2,
Phat = 1/(xi(1-xi)) - 2 + xi(1-xi)  (K = q for the z,x cases, K = k for delta).
It is subtracted with kappa^2/(kappa^2+M^2), M^2 = xi(1-xi) K^2, and integrated in DR (MS-bar):
   int d^dkappa/(2pi)^d 1/(kappa^2+M^2) = (1/4pi) [2/eps - ln(M^2/mu^2)].
Output for fixed (k,l,xi): finite part (MS-bar, given mu) and the coefficient of (1/4pi)(2/eps), both as
contributions to R (ratio to LO, units g^2).
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
from v0_common import Colour, V0, beta, Amom
from scipy.integrate import quad

C = Colour(3)
f, t = C.f, C.t
AC = np.einsum('eij,djk->edik', t, t) + np.einsum('dij,ejk->edik', t, t)   # {t^e,t^d}
T = {
 'z1': np.einsum('adc,cxe,edy,yij->axij', f, f, f, t),
 'z2': np.einsum('adc,cxe,edij->axij', f, f, AC),
 'x1': np.einsum('adc,dxb,cby,yij->axij', f, f, f, t),
 'x2': np.einsum('adc,dxb,cbij->axij', f, f, AC),
 'd1': -np.einsum('abc,cby,yxd,dij->axij', f, f, f, t),
 'd2': -(np.einsum('abc,cxd,dbij->axij', f, f, AC) + np.einsum('abc,bxd,cdij->axij', f, f, AC)),
}
P = {s: C.proj(T[s]) for s in T}

def Ycontract(X, kap, xi):
    return kap * np.trace(X) - (kap @ X) / xi - (X @ kap) / (1 - xi)
    # kap@X : sum_j kap_j X[j,i] ; X@kap : sum_k X[i,k] kap_k

def pref(xi):
    return 2 / (2 * np.pi) * (-1j * np.sqrt(xi * (1 - xi)) / np.sqrt(2))

def W_cases(xi, kap, k, l):
    """returns dict case -> W vector (kinematic, before colour), integrand of int d^2kappa/(2pi)^2"""
    out = {}
    pr = pref(xi) / (kap @ kap)
    # z: a = kap + l + (1-xi) k, b = xi k - kap
    a = kap + l + (1 - xi) * k; b = xi * k - kap
    out['z1'] = pr * Ycontract(beta(1 - xi, a, xi, b), kap, xi)
    out['z2'] = pr * Ycontract(0.25 * np.outer(Amom(xi, b), Amom(1 - xi, a)), kap, xi)
    # x: a = (1-xi)k + kap, b = xi k - kap + l
    a = (1 - xi) * k + kap; b = xi * k - kap + l
    out['x1'] = pr * Ycontract(beta(1 - xi, a, xi, b), kap, xi)
    out['x2'] = pr * Ycontract(0.25 * np.outer(Amom(xi, b), Amom(1 - xi, a)), kap, xi)
    # delta: a = (1-xi)k + kap, b = xi k - kap
    a = (1 - xi) * k + kap; b = xi * k - kap
    out['d1'] = pr * Ycontract(beta(1 - xi, a, xi, b), kap, xi)
    out['d2'] = pr * Ycontract(0.25 * np.outer(Amom(xi, b), Amom(1 - xi, a)), kap, xi)
    return out

def Phat(xi):
    return 1 / (xi * (1 - xi)) - 2 + xi * (1 - xi)

def W_UV_unavg(xi, kap, K):
    """leading large-kappa term (case independent apart from K), unaveraged, d=2."""
    k2 = kap @ kap
    X = -1j * np.sqrt(xi * (1 - xi)) / ((K @ K) * k2) * (-(K @ kap) * np.eye(2) + np.outer(K, kap) / (1 - xi) + np.outer(kap, K) / xi)
    return pref(xi) / k2 * Ycontract(X, kap, xi)

def R_of(Ws, v0, keys):
    return sum(2 * (P[s] * (v0 @ Ws[s])).real for s in keys) / (np.sqrt(2) * (v0 @ v0))

def rowV_xi(xi, k, l, mu=None, eps=1e-7, verbose=False):
    if mu is None: mu = np.linalg.norm(k)
    q = k + l; v0 = V0(k, l)
    Kz = q; Kd = k
    Mz2 = xi * (1 - xi) * (q @ q); Md2 = xi * (1 - xi) * (k @ k)
    colUVq = P['z1'] + P['x1']; colUVk = P['d1']
    def integrand(r, th):
        kap = np.array([r * np.cos(th), r * np.sin(th)])
        Ws = W_cases(xi, kap, k, l)
        val = R_of(Ws, v0, Ws.keys())
        # subtraction
        wq = W_UV_unavg(xi, kap, Kz) * r * r / (r * r + Mz2)
        wk = W_UV_unavg(xi, kap, Kd) * r * r / (r * r + Md2)
        sub = 2 * ((colUVq * (v0 @ wq)) + (colUVk * (v0 @ wk))).real / (np.sqrt(2) * (v0 @ v0))
        return (val - sub) * r / (2 * np.pi) ** 2
    inner = lambda r: quad(lambda th: integrand(r, th), 0, 2 * np.pi, limit=200, epsabs=eps, epsrel=1e-8)[0]
    # radial: split at the scales
    scales = sorted(set([np.linalg.norm(xi * k), np.linalg.norm((1 - xi) * k), np.linalg.norm(l + (1 - xi) * k),
                         np.linalg.norm(xi * k + l), np.linalg.norm(q), 1e-9]))
    pts = [0.0] + scales + [10 * max(scales), 100 * max(scales)]
    tot = 0.0
    for r0, r1 in zip(pts[:-1], pts[1:]):
        tot += quad(inner, r0, r1, limit=200, epsabs=eps, epsrel=1e-7)[0]
    tot += quad(inner, pts[-1], np.inf, limit=200, epsabs=eps, epsrel=1e-7)[0]
    # analytic DR of the subtraction (averaged W_UV = (2/(2pi))(1/sqrt2) Phat K/K^2 / kappa^2)
    def uvR(col, K, M2):
        w = (2 / (2 * np.pi)) / np.sqrt(2) * Phat(xi) * K / (K @ K)
        return 2 * (col * (v0 @ w)).real / (np.sqrt(2) * (v0 @ v0))
    fin_sub = (uvR(colUVq, Kz, Mz2) * (-np.log(Mz2 / mu ** 2)) + uvR(colUVk, Kd, Md2) * (-np.log(Md2 / mu ** 2))) / (4 * np.pi)
    pole = (uvR(colUVq, Kz, Mz2) + uvR(colUVk, Kd, Md2))       # coefficient of (1/4pi)(2/eps)
    return tot + fin_sub, pole, tot

if __name__ == '__main__':
    print("colour projections:", {s: np.round(P[s], 6) for s in P})
    k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7])
    for xi in [float(a) for a in sys.argv[1:]] or [0.3]:
        fin, pole, tot = rowV_xi(xi, k, l)
        print(f"{xi:.6f} fin {fin:.10f} pole {pole:.10f} numeric {tot:.10f}", flush=True)
