#!/usr/bin/env python3
"""dilute2 / e1: soft limit of the paper's dilute two-gluon cross section (Groups II+III) at fixed s = p+/k+ -> 0.

Claim: for p+ = s k+, |p| = O(s|k|) -> 0 (second gluon soft, at a finite rapidity distance from the observed one),
   A_gluons/(qa^2 qb^2)  ->  (1/k^2) S(p; p_c),   S = (p^2 + p_c^2 - p.p_c) / (p^2 (p-p_c)^2),   p_c = s k,
i.e. S - 1/p^2 = p.p_c/(p^2 (p-p_c)^2)  [= (1/p^2) (1/2)(1 + cosh(Dy)/(cosh(Dy)-cos(phi))) - 1/p^2, e^{Dy} = |p_c|/|p|].
Checked (i) with the paper's amplitude (d1.cgc_sq, Eq. (A2gamp) at O(u)), (ii) with BSV's A_gluons, also for D != 4.
Then: (1/pi) int d^2p [S - 1/p^2 - chi-subtraction] against the DR result  -2/eps + ln(p_c^2/mu^2)."""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import d1_qmrk_integrand as d1
from scipy.integrate import dblquad

def S(p, pc):
    return (p @ p + pc @ pc - p @ pc) / ((p @ p) * ((p - pc) @ (p - pc)))

k = np.array([1.0, 0.0]); qa = np.array([1.3, 0.7]); kp = 1.0
T = d1.colour_tensors(3); CA = 9 * 8
print("soft limit: ratio  k^2 A/(qa^2 qb^2) / S   [paper amplitude | BSV D=4 | BSV D=5 | BSV D=3.4]")
for v in (np.array([0.4, 0.3]), np.array([1.3, -0.2]), np.array([-0.7, 0.5]), np.array([1.05, 0.1])):
    for s in (1e-2, 1e-3, 1e-4, 1e-5):
        p = s * v; pp = s * kp; x = kp / (kp + pp)
        qb = k + p - qa; l = -qb
        amp = d1.cgc_sq(kp, k, pp, p, l, T) * kp * pp * (qa @ qa) ** 2 / CA      # = A_gluons by Eq. (d:identity)
        r = [amp] + [d1.A_gluons(qa, qb, k, p, x, D=D) for D in (4.0, 5.0, 3.4)]
        r = [(k @ k) * a / ((qa @ qa) * (qb @ qb)) / S(p, s * k) for a in r]
        print(f"  p/s={v}, s={s:.0e}: " + "  ".join(f"{a:.8f}" for a in r))

# DR integral of S - 1/p^2:  (1/pi) int d^{2-eps}p [S - 1/p^2] = -2/eps + ln(pc^2/mu^2)   (MS-bar).
# 2D check with a cut-off subtraction: (1/pi) int d^2p [S - 1/p^2 - theta(|p-pc|<lam)/(p-pc)^2] = ln(pc^2/lam^2).
pc = np.array([1.0, 0.0])
for lam in (0.5, 0.2):
    def f(r, th, cx, cy, sub):
        p = np.array([cx + r * np.cos(th), cy + r * np.sin(th)])
        val = S(p, pc) - 1 / (p @ p)
        if sub and r < lam: val -= 1 / r ** 2
        return val * r
    # disc around pc (radius lam) in polar coordinates around pc; rest in polar coordinates around pc as well
    I1 = dblquad(lambda th, r: f(r, th, 1.0, 0.0, True), 0, lam, 0, 2 * np.pi, epsabs=1e-10, epsrel=1e-10)[0]
    I2 = dblquad(lambda th, r: f(r, th, 1.0, 0.0, False), lam, 1.0, 0, 2 * np.pi, epsabs=1e-10, epsrel=1e-10)[0]
    I3 = dblquad(lambda th, u: f(1 / u, th, 1.0, 0.0, False) / u ** 2, 0, 1.0, 0, 2 * np.pi, epsabs=1e-10, epsrel=1e-10)[0]
    print(f"lam={lam}: (1/pi) int [S - 1/p^2 - theta/(p-pc)^2] = {(I1 + I2 + I3) / np.pi:.6f}   expected ln(pc^2/lam^2) = {np.log(1 / lam ** 2):.6f}")
