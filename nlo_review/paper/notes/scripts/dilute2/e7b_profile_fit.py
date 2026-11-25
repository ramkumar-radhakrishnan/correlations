#!/usr/bin/env python3
"""dilute2 / e7b: is the soft-region profile of the resolved Row I.V of the single-scale form
   <H_V>(u) = H_inf + sigma [ P_k psi(u) - P_q psi(u k^2/q_a^2) ],   psi(u) = 1/(1+u),  sigma = 8 pi abar ?
(then its u-dependence is cancelled scale by scale by the commutator part of Row I.II, weights P_k, -P_q, and the hard cut
s > Lambda differs from dimensional regularisation in this region by a pure constant only)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e7_region_a as e7
from e0_points import PTS
abar = 3 / (4 * np.pi ** 2); sigma = 8 * np.pi * abar
for name in sys.argv[1:]:
    P = PTS[name]; k = np.array([1.0, 0.0]); l = P['l']; q = k + l
    v0 = k / (k @ k) - q / (q @ q); vv = v0 @ v0; Pk = (v0 @ k) / vv; Pq = (v0 @ q) / (q @ q) / vv; rq = q @ q
    phis = 2 * np.pi * (np.arange(32) + 0.5) / 32
    print(f"# {name}: P_k={Pk:.5f} P_q={Pq:.5f} q_a^2/k^2={rq:.5f};  ln u, <H_V>, model, difference (two values of rho)")
    out = []
    for t in (-8, -5, -3, -2, -1, -0.5, 0.0, 0.5, 1, 2, 3, 5, 8, 12):
        u = np.exp(t); row = [t]
        for rho in (1e-4, 1e-5):
            hv = np.mean([e7.H_rows(u, ph, k, l, rho)[2] for ph in phis]); row.append(hv)
        out.append(row)
    out = np.array(out); Hinf = out[-1, 2]
    for t, h1, h2 in out:
        u = np.exp(t); model = Hinf + sigma * (Pk / (1 + u) - Pq / (1 + u / rq))
        print(f"{t:6.2f}  {h1: .6f}  {h2: .6f}  model {model: .6f}  diff {h2 - model: .6f}")
