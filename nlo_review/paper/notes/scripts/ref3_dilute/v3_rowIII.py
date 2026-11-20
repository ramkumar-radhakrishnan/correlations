#!/usr/bin/env python3
"""
ref3_dilute / v3: Row I.III (A_B3) in the dilute limit, point-like quark, momentum space (k+ = 1, g = 1).
A_B3^a_i = int ds/(2pi) int d^2p/(2pi)^2 [B3^{ba}_{ji}(k,p)]^dag (sqrt2/sqrt s) u^{bd}(l) t^d V_j(p,l),
V(p,l) = p/p^2 - (p+l)/(p+l)^2,  s = p+/k+, B3 from Eq. (B3mom) (rho -> t, dagger: i -> -i).
Colour structures: F = f^{abc} t^c u^{bd} t^d,  AC = {t^a,t^b} u^{bd} t^d.
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
from v0_common import Colour, V0

C = Colour(3)
f, t = C.f, C.t
T_F = np.einsum('abc,cij,bxd,djk->axik', f, t, f, t)
AC = np.einsum('aij,bjk->abik', t, t) + np.einsum('bij,ajk->abik', t, t)
T_AC = np.einsum('abij,bxd,djk->axik', AC, f, t)
P_F, P_AC = C.proj(T_F), C.proj(T_AC)

def Vv(p, l):
    pl = p + l
    return p / (p @ p) - pl / (pl @ pl)

def K_F(s, p, k):
    """kinematic matrix [j,i] multiplying colour F (dagger already taken); T1 + T2(f part) + T3"""
    k2, p2 = k @ k, p @ p
    I = np.eye(2)
    T1 = -1j * np.sqrt(s) * (1 + s) / (s - 1) ** 2 * I / (s * k2 - p2)
    Y = (p2 + k2 * s) / (p2 - k2 * s)
    T2 = -np.outer(p, k) / (k2 * p2 * np.sqrt(s)) * 1j * Y
    n = k - p; w = k * s - p
    br = np.outer(w, n) / s - np.outer(n, w) - (n @ w) * I / (s - 1)
    if s > 1:
        G = (s - 1) / (n @ n)
    else:
        d = p - s * k
        G = -s * (1 - s) / (d @ d)
    T3 = -2j * np.sqrt(s) / (1 - s) ** 2 / (k2 * s - p2) * br * G
    return T1 + T2 + T3

def K_AC(s, p, k):
    return -np.outer(p, k) / ((k @ k) * (p @ p) * np.sqrt(s))

def Rdens(s, p, k, l):
    """contribution to R per ds d^2p/(2pi)^2 (integrand)."""
    v0 = V0(k, l); V = Vv(p, l)
    WF = np.sqrt(2) / (2 * np.pi * np.sqrt(s)) * (V @ K_F(s, p, k))      # sum_j V_j K[j,i]
    WA = np.sqrt(2) / (2 * np.pi * np.sqrt(s)) * (V @ K_AC(s, p, k))
    return 2 * (P_F * (v0 @ WF) + P_AC * (v0 @ WA)).real / (np.sqrt(2) * (v0 @ v0))

if __name__ == '__main__':
    print("projections P_F, P_AC:", P_F, P_AC)
    k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7])
    p = np.array([0.4, 0.5])
    print("s0 = p^2/k^2 =", p @ p)
    for s in (1 - 1e-2, 1 - 1e-3, 1 - 1e-4, 1 + 1e-4, 1 + 1e-3, 1 + 1e-2):
        print(f"s={s}: R*(s-1) = {Rdens(s,p,k,l)*(s-1):.6e}   R*(s-1)^2 = {Rdens(s,p,k,l)*(s-1)**2:.3e}")
    for s in (1e-4, 1e-6, 1e4, 1e6):
        print(f"s={s}: s*R = {s*Rdens(s,p,k,l):.8e}")
    s0 = p @ p
    for d in (1e-3, 1e-5):
        print(f"near s0: (s-s0)R at s0+-{d}: {(d)*Rdens(s0+d,p,k,l):.6e} {(-d)*Rdens(s0-d,p,k,l):.6e}")
