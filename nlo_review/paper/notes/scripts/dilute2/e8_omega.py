#!/usr/bin/env python3
"""dilute2 / e8: moments of the Regge-region integrand of the virtual rows (hard cut scheme),
   tau3(p; k, l) = lim_{s->0} s [F_III + F_IV + F_V(resolved)](s, p)     (p = momentum of the soft loop gluon after the target),
   J_0 = int d^2p/(2pi)^2 [tau3 - c_0 th(|p|<lam)/p^2 - c_inf th(|p|>lam)/p^2],
   J_1 = the same with an extra factor ln(p^2/lam^2)   (tails c_0, c_inf: universal numbers, measured).
J_0 must reproduce the sum of the L coefficients of Rows I.III, I.IV, I.V (up to a universal constant): check of the method.
J_1 is the kinematics-dependent part of  int tau ln p^2, i.e. of the change of the constant C when the lower cut on the loop gluon
is put on its light-cone energy p^- (p+ > p^2/M) instead of p+ > Lambda, as far as Rows I.III-I.V are concerned.
usage: e8_omega.py name [nr nth]   (corner name of e0_points; computed in the frame k=(1,0), lam = |k| = 1)"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import e7c_region_b as eb
from e0_points import PTS
import e2_real_cut as e2

def tau3(p, k, l, s=1e-9):
    return sum(eb.Hb(p / s, k, l, s)) / (p @ p)

def moments(k, l, nr=6, nth=4, c0=None, cinf=None):
    q = k + l
    if c0 is None:
        c0 = np.mean([1e-8 * tau3(1e-4 * np.array([np.cos(t), np.sin(t)]), k, l, 1e-12) for t in np.linspace(0.1, 6.2, 12)])
        cinf = np.mean([1e8 * tau3(1e4 * np.array([np.cos(t), np.sin(t)]), k, l, 1e-9) for t in np.linspace(0.1, 6.2, 12)])
    sing = [k, -k, -l, -k - l, q, -q]
    scales = sorted(set(round(float(np.linalg.norm(v)), 12) for v in sing))
    ed = {1e-5, 1e5}
    for sc in scales:
        for j in range(0, 12):
            ed.add(sc * (1 + 2.0 ** (-j))); ed.add(sc * (1 - 2.0 ** (-j - 1)))
        ed.add(sc * 4); ed.add(sc * 16); ed.add(sc / 4); ed.add(sc / 16)
    ed.add(1.0)
    ed = sorted(e for e in ed if e > 0)
    r, wr = e2.gl_nodes(ed, nr); th, wt = e2.ang_nodes([np.arctan2(v[1], v[0]) for v in sing], nth)
    J0 = J1 = 0.0
    for ri, wi in zip(r, wr):
        sub = (c0 if ri < 1 else cinf) / ri ** 2
        acc = 0.0
        for t, w in zip(th, wt):
            acc += w * (tau3(ri * np.array([np.cos(t), np.sin(t)]), k, l) - sub)
        J0 += wi * ri * acc; J1 += wi * ri * acc * np.log(ri ** 2)
    return J0 / (2 * np.pi) ** 2, J1 / (2 * np.pi) ** 2, c0, cinf

if __name__ == '__main__':
    name = sys.argv[1]; nr = int(sys.argv[2]) if len(sys.argv) > 2 else 6; nth = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    P = PTS[name]; k = np.array([1.0, 0.0]); l = P['l']
    J0, J1, c0, cinf = moments(k, l, nr, nth)
    print(f"{name} nr={nr} nth={nth}: c_0 = {c0:.6f}  c_inf = {cinf:.6f}  J_0 = {J0:.6f}  J_1 = {J1:.6f}", flush=True)
