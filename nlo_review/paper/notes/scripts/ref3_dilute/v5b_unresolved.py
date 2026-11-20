#!/usr/bin/env python3
"""
ref3_dilute / v5b: the unresolved part of Row I.V versus the self-energy part of Row I.II (point quark, dilute).
Unresolved I.V: the target colour rotation acts on the pair at its centre of mass, i.e. the kick l is shared
according to the momentum fractions and the relative momentum kappa is unchanged:
   cases z+x (colour P_z1+P_x1):  a = (1-xi) q + kappa, b = xi q - kappa   (pair before the collision has total q)
   case delta (colour P_d1):      a = (1-xi) k + kappa, b = xi k - kappa
Self-energy of Row I.II (app_wavefunctions, Psi^1_{g rho}):  c_SE(K) = (alpha_s N_c/2pi) int dxi Phat(xi)[-2/eps + ln(K^2 xi(1-xi)/mu^2)].
The subtraction term of v5 integrates (DR, MS-bar) to exactly -R_SE(xi).  Hence
     R_unres(xi) + R_SE(xi) = int d^2kappa/(2pi)^2 [unres - S]   (2D, finite)
and the claim "unresolved I.V cancels the self-energy identically" is the statement that this integral vanishes.
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
from v0_common import V0, beta, Amom
import v5_rowV as v5
from scipy.integrate import quad

def resid(xi, k, l, eps=1e-8):
    q = k + l; v0 = V0(k, l)
    Mq2 = xi * (1 - xi) * (q @ q); Mk2 = xi * (1 - xi) * (k @ k)
    cq = v5.P['z1'] + v5.P['x1']; ck = v5.P['d1']
    def integrand(r, th):
        kap = np.array([r * np.cos(th), r * np.sin(th)])
        pr = v5.pref(xi) / (kap @ kap)
        Wq = pr * v5.Ycontract(beta(1 - xi, (1 - xi) * q + kap, xi, xi * q - kap), kap, xi)
        Wk = pr * v5.Ycontract(beta(1 - xi, (1 - xi) * k + kap, xi, xi * k - kap), kap, xi)
        Sq = v5.W_UV_unavg(xi, kap, q) * r * r / (r * r + Mq2)
        Sk = v5.W_UV_unavg(xi, kap, k) * r * r / (r * r + Mk2)
        val = 2 * (cq * (v0 @ (Wq - Sq)) + ck * (v0 @ (Wk - Sk))).real / (np.sqrt(2) * (v0 @ v0))
        return val * r / (2 * np.pi) ** 2
    inner = lambda r: quad(lambda th: integrand(r, th), 0, 2 * np.pi, limit=200, epsabs=eps, epsrel=1e-9)[0]
    sc = sorted(set([np.linalg.norm(xi * k), np.linalg.norm((1 - xi) * k), np.linalg.norm(xi * q), np.linalg.norm((1 - xi) * q)]))
    pts = [0.0] + sc + [10 * sc[-1], 100 * sc[-1]]
    tot = sum(quad(inner, a, b, limit=200, epsabs=eps, epsrel=1e-8)[0] for a, b in zip(pts[:-1], pts[1:]))
    tot += quad(inner, pts[-1], np.inf, limit=200, epsabs=eps, epsrel=1e-8)[0]
    # separate the two K-structures to see whether each vanishes
    return tot

if __name__ == '__main__':
    k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7])
    for xi in (0.05, 0.3, 0.5, 0.8):
        r = resid(xi, k, l)
        print(f"xi={xi}: int[unres - S] = {r:.3e}   (compare |R_SE(xi)| scale: Phat*N_c/(4pi^2) = {v5.Phat(xi)*3/(4*np.pi**2):.3e})", flush=True)
