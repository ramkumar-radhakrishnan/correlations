#!/usr/bin/env python3
"""dilute2 / e2b: (i) convergence of j_fin(s) (quadrature multiplier m=1,2), (ii) fit of the small-s behaviour
j_fin = a ln s + b (claim a = 2, b = ln(k^2/mu^2) = 0), (iii) I(L) = int_{e^-L}^inf ds/s j_fin(s) against -L^2 + const."""
import numpy as np, sys
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/dilute2')
import e2_real_cut as e2
from scipy.interpolate import CubicSpline
k = np.array([1.0, 0.0]); qa = np.array([1.3, 0.7]); obj = e2.Real(k, qa)
print("(i) convergence: u = ln s, j_fin for m = 1, 2, 3")
for u in (-12.0, -10.0, -8.0, -5.0, -2.0, 0.0, 1.0, 3.0):
    print(f"   u={u:6.1f}: " + "  ".join(f"{e2.j_fin(obj, np.exp(u), m=m)[0]: .9f}" for m in (1, 2, 3)), flush=True)
d = np.loadtxt('e2_table_m2.dat'); u, c, J = d[:, 0], d[:, 1], d[:, 2]
sel = (u >= -12.01) & (u <= -7.99)
a, b = np.polyfit(u[sel], J[sel], 1)
print(f"(ii) fit j_fin = a u + b on u in [-12,-8]: a = {a:.6f} (claim 2), b = {b:.6f} (claim 0 up to O(s ln^2 s))")
sel2 = (u >= -12.01) & (u <= -9.99); a, b = np.polyfit(u[sel2], J[sel2], 1)
print(f"     same on u in [-12,-10]: a = {a:.6f}, b = {b:.6f};   collinear residue c(s) at u=-12,-10,-8: {c[u==-12][0]:.7f} {c[u==-10][0]:.7f} {c[u==-8][0]:.7f} (claim -> 1)")
# remainder r(u) = j_fin - theta(u<0) 2u ; I(L) + L^2 = int_{-L}^inf du r(u)
ok = u >= -12.01
r = J[ok] - np.where(u[ok] < 0, 2 * u[ok], 0.0); uu = u[ok]
cs1 = CubicSpline(uu[uu <= 0], r[uu <= 0]); cs2 = CubicSpline(uu[uu >= 0], r[uu >= 0])
tail = r[-1] / 0.8   # r ~ exp(-0.8 u) at large u (fitted below)
sl = np.polyfit(uu[-6:], np.log(-r[-6:]), 1)[0]; tail = r[-1] / (-sl)
print(f"(iii) large-s tail: r ~ exp({sl:.3f} u);  tail beyond u={uu[-1]}: {tail:.6f}")
for L in (2, 3, 4, 5, 6, 8, 10, 12):
    val = cs1.integrate(-L, 0) + cs2.integrate(0, uu[-1]) + tail
    print(f"   L={L:3d}:  I(L) = {val - L * L: .6f}   I(L) + L^2 = {val: .6f}")
