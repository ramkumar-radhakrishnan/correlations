#!/usr/bin/env python3
"""
ref3_dilute / v3b: Row I.III integrated (point quark, dilute), Lambda -> 0, vee -> infinity with logs kept:
   R_III = R_T2 + R_13,
 T2 part ({t,t} and f parts of the second line of B3), integrated over s analytically (PV):
   R_T2 = -(3/2pi) (V0.k)/(k^2 V0^2) * J,   J = int d^2p/(2pi)^2 g(p) ln|p^2/(Lambda k^2) - 1|,  g = (p.V)/p^2
        = -(3/2pi) (V0.k)/(k^2 V0^2) * { (1/4pi)[ (ln(lam^2/k^2)+L)^2/2 - pi^2/3 ] + J1 + L J0 }
   J1 = int [g - theta(p<lam)/p^2] ln(p^2/k^2),  J0 = int [g - theta(p<lam)/p^2]       (d^2p/(2pi)^2), L = ln(k+/Lambda)
 T1+T3 part (f only), s integral numerically at fixed p (PV at s0=p^2/k^2, |s-1|>Lambda):
   h13(p) = Int(p) + 2 c(p) L + A_inf(p) ln(vee/k+),  R_13 = int d^2p/(2pi)^2 h13.
 The on-shell (i pi delta) parts of B3 do not contribute to R (checked: purely imaginary relative to LO).
Output: R_III = a2 L^2 + a1 L + a_vee ln(vee) + a0.
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v3_rowIII as v3
from v0_common import V0
from scipy.integrate import quad

def make(k, l):
    v0 = V0(k, l); k2 = k @ k
    def F13(s, p):
        p2 = p @ p
        if s == 1.0 or s * k2 == p2:
            s = s * (1 + 1e-12) + 1e-14
        V = v3.Vv(p, l)
        Y = (p2 + k2 * s) / (p2 - k2 * s)
        T2 = -np.outer(p, k) / (k2 * p2 * np.sqrt(s)) * 1j * Y
        W = np.sqrt(2) / (2 * np.pi * np.sqrt(s)) * (V @ (v3.K_F(s, p, k) - T2))
        return 2 * (v3.P_F * (v0 @ W)).real / (np.sqrt(2) * (v0 @ v0))
    return F13, v0

def h13(p, F13, k, a=0.5):
    k2 = k @ k; s0 = (p @ p) / k2
    if abs(s0 - 1) < 1e-7:                      # measure-zero circle |p| = |k|: nudge
        p = p * (1 + 2e-7); s0 = (p @ p) / k2
    e = 1e-7
    c = 0.5 * (F13(1 + e, p) * e + F13(1 - e, p) * e)          # F ~ c/|s-1| on both sides
    S = max(20.0, 4 * s0)
    Ainf = F13(1e9, p) * 1e9
    def seg(f, x0, x1):
        if x0 < s0 < x1 and min(s0 - x0, x1 - s0) > 1e-9:
            g = lambda s: f(s) * (s - s0)
            return quad(g, x0, x1, weight='cauchy', wvar=s0, limit=200)[0]
        return quad(f, x0, x1, limit=200, points=[s0] if x0 < s0 < x1 else None)[0]
    f = lambda s: F13(s, p)
    def fm(s):
        if abs(s - 1) < 1e-13:
            s = 1 + 1e-13
        return F13(s, p) - c / abs(s - 1)
    tot = seg(f, 1e-13, 1 - a) + seg(fm, 1 - a, 1.0) + seg(fm, 1.0, 1 + a) + seg(f, 1 + a, S)
    tot += quad(lambda s: F13(s, p) - Ainf / s, S, np.inf, limit=200)[0]
    tot += 2 * c * np.log(a) - Ainf * np.log(S)
    return tot, c, Ainf

def p_integral(fun, k, l, eps=1e-6, ncomp=3):
    """int d^2p/(2pi)^2 of a vector-valued fun(p), polar coordinates, vector adaptive quadrature (quad_vec)."""
    from scipy.integrate import quad_vec
    sc = sorted(set([np.linalg.norm(k), np.linalg.norm(l), np.linalg.norm(k + l), 0.5 * np.linalg.norm(k)]))
    pts = [0.0] + sc + [3 * sc[-1], 10 * sc[-1], 40 * sc[-1]]
    inner = lambda r: quad_vec(lambda th: fun(np.array([r * np.cos(th), r * np.sin(th)])), 0, 2 * np.pi,
                               epsabs=eps, epsrel=1e-5, limit=200)[0] * r
    acc = np.zeros(ncomp)
    for r0, r1 in zip(pts[:-1], pts[1:]):
        acc += quad_vec(inner, r0, r1, epsabs=eps, epsrel=1e-5, limit=200)[0]
    acc += quad_vec(inner, pts[-1], np.inf, epsabs=eps, epsrel=1e-5, limit=200)[0]
    return acc / (2 * np.pi) ** 2

def rowIII(k, l):
    F13, v0 = make(k, l)
    k2 = k @ k
    def fun13(p):
        return np.array(h13(p, F13, k))
    I13 = p_integral(fun13, k, l)           # [Int, c, Ainf]
    lam = 1.0
    def g(p):
        V = v3.Vv(p, l); p2 = p @ p
        gg = (p @ V) / p2 - (1.0 / p2 if p2 < lam ** 2 else 0.0)
        return np.array([gg * np.log(p2 / k2), gg])
    J = p_integral(g, k, l, ncomp=2)
    pref = -(3 / (2 * np.pi)) * (v0 @ k) / (k2 * (v0 @ v0))
    Lk = np.log(lam ** 2 / k2)
    a2 = pref / (8 * np.pi)
    a1 = pref * (Lk / (4 * np.pi) + J[1]) + 2 * I13[1]
    avee = I13[2]
    a0 = pref * (Lk ** 2 / (8 * np.pi) - np.pi / 12 + J[0]) + I13[0]
    return dict(a2=a2, a1=a1, avee=avee, a0=a0, J=J, I13=I13)

if __name__ == '__main__':
    k = np.array([float(sys.argv[1]), float(sys.argv[2])]) if len(sys.argv) > 4 else np.array([1.0, 0.0])
    l = np.array([float(sys.argv[3]), float(sys.argv[4])]) if len(sys.argv) > 4 else np.array([0.3, 0.7])
    r = rowIII(k, l)
    print(f"k={k} l={l}: R_III = {r['a2']:.8f} L^2 + {r['a1']:.8f} L + {r['avee']:.8f} ln(vee) + {r['a0']:.8f}")
    print("   J1,J0 =", r['J'], " I13(Int,c,Ainf) =", r['I13'])
