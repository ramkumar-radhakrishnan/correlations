#!/usr/bin/env python3
"""dilute2 / e7c: the soft region p = s v, s -> 0 (soft loop gluon at the rapidity of the observed gluon) of the virtual
rows: H_b(v) = lim s p^2 [F_III + F_IV + F_V(resolved)](s, p = s v).  Is it the same function of v at different (q_a,q_b),
i.e. universal (eikonal), so that the hard cut s > Lambda differs from DR by the universal soft-collinear term only?
Compared with the real soft factor:  -4 pi abar [p^2 S(p;p_c) - 1],  S = (p^2+p_c^2-p.p_c)/(p^2 (p-p_c)^2), p_c = s k."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v3c_fast as v3, v4_rowIV as v4, v5_rowV as v5
from v0_common import V0
from e0_points import PTS
abar = 3 / (4 * np.pi ** 2)
def Hb(v, k, l, s):
    p = s * v; k2 = k @ k; q = k + l; v0 = V0(k, l)
    h3 = s * (p @ p) * v3.Fvec(np.array([s]), p, k, l)[0]
    h4 = s * (p @ p) * v4.Rdens(s, p, k, l)
    def v5dens(xi, kap):
        Ws = v5.W_cases(xi, kap, k, l); val = v5.R_of(Ws, v0, Ws.keys())
        Mz2 = xi * (1 - xi) * (q @ q); Md2 = xi * (1 - xi) * k2; rr = kap @ kap
        wq = v5.W_UV_unavg(xi, kap, q) * rr / (rr + Mz2); wk = v5.W_UV_unavg(xi, kap, k) * rr / (rr + Md2)
        sub = 2 * (((v5.P['z1'] + v5.P['x1']) * (v0 @ wq)) + (v5.P['d1'] * (v0 @ wk))).real / (np.sqrt(2) * (v0 @ v0))
        return val - sub
    # Row V: the soft daughter has transverse momentum b = xi k - kappa (cases z, d), i.e. kappa = xi k - p
    h5 = s * (p @ p) * (v5dens(s, s * k - p) + v5dens(1 - s, -(s * k - p)))
    return h3, h4, h5
k = np.array([1.0, 0.0])
vs = [np.array(x) for x in ((0.05, 0.02), (0.3, 0.4), (0.9, 0.3), (1.0, 0.05), (1.3, -0.4), (-0.6, 0.5), (2.5, 1.0), (-4.0, 3.0), (20.0, 10.0), (200.0, -100.0))]
print("columns: v ; for each point: H_III, H_IV, H_V, sum ;  last: -4 pi abar [p^2 S - 1]")
for v in vs:
    line = f"v=({v[0]:7.2f},{v[1]:7.2f}) "
    for name in ('a1b1', 'a2b2', 'a1b2'):
        l = PTS[name]['l']
        for s in (1e-6,):
            h = Hb(v, k, l, s); line += f"| {name}: {h[0]: .5f} {h[1]: .5f} {h[2]: .5f} sum {sum(h): .5f} "
    S = (v @ v + 1 - v[0]) / ((v @ v) * ((v - k) @ (v - k)))
    line += f"| eik {-4 * np.pi * abar * ((v @ v) * S - 1): .5f}"
    print(line)
