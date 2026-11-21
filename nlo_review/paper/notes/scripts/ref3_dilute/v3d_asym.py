#!/usr/bin/env python3
"""ref3_dilute / v3d: asymptotic coefficients of Row I.III that generate ln^2(k+/Lambda):
 soft region  Lambda k^2 << p^2 << k^2:      p^2 h ~ alpha(th) L + beta(th) ln p^2 + gamma(th)
 forward region Lambda << |p-k| << |k|:      |p-k|^2 h ~ a(th) ln(|p-k|^2/Lambda^2) + c(th)
 angular averages (Gauss-Legendre with refinement at +-k, -l); L^2 coefficient = (abar - bbar/2)/(4 pi) + afwd/(2 pi)."""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v3c_fast as f
k = np.array([float(a) for a in sys.argv[1:3]]); l = np.array([float(a) for a in sys.argv[3:5]])
th, w = f.angular_nodes([np.arctan2(k[1], k[0]), np.arctan2(k[1], k[0]) + np.pi, np.arctan2(-l[1], -l[0])], 6)
w = w / (2 * np.pi)
def avg(fun):
    return np.sum([wi * fun(np.array([np.cos(t), np.sin(t)])) for t, wi in zip(th, w)])
# soft: p = r e, two r and two Lambda
r1, r2 = 1e-3, 1e-4; L1, L2 = 1e-10, 1e-12
S = {(r, La): avg(lambda e: f.h(r * e, k, l, La)[0] * r * r) for r in (r1, r2) for La in (L1, L2)}
alpha = (S[(r1, L2)] - S[(r1, L1)]) / np.log(L1 / L2)
beta = (S[(r1, L1)] - S[(r2, L1)]) / np.log(r1 ** 2 / r2 ** 2)
alpha2 = (S[(r2, L2)] - S[(r2, L1)]) / np.log(L1 / L2); beta2 = (S[(r1, L2)] - S[(r2, L2)]) / np.log(r1 ** 2 / r2 ** 2)
print(f"soft: alpha={alpha:.6f} ({alpha2:.6f})  beta={beta:.6f} ({beta2:.6f})", flush=True)
# forward: p = k + d e
d1, d2 = 1e-4, 1e-5; La = 1e-9
thk, wk = f.angular_nodes([np.arctan2(k[1], k[0]), np.arctan2(-k[1], -k[0]), np.arctan2(-(l + k)[1], -(l + k)[0])], 6)
wk = wk / (2 * np.pi)
def avgk(fun):
    return np.sum([wi * fun(np.array([np.cos(t), np.sin(t)])) for t, wi in zip(thk, wk)])
Fd = {d: avgk(lambda e: f.h(k + d * e, k, l, La)[0] * d * d) for d in (d1, d2)}
a = (Fd[d1] - Fd[d2]) / np.log(d1 ** 2 / d2 ** 2)
cfw = Fd[d1] - a * np.log(d1 ** 2 / La ** 2)
FdL = avgk(lambda e: f.h(k + d1 * e, k, l, 1e-10)[0] * d1 * d1)
print(f"forward: a={a:.6f}  c={cfw:.6f}; check Lambda-slope: {(FdL - Fd[d1]) / (2*np.log(1e-9/1e-10)):.6f}", flush=True)
L2coef = (alpha - beta / 2) / (4 * np.pi) + a / (2 * np.pi)
print(f"L^2 coefficient of R_III = {(alpha - beta/2)/(4*np.pi):.6f} (soft) + {a/(2*np.pi):.6f} (forward) = {L2coef:.6f}")
