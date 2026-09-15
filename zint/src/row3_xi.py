"""The UV coefficient as a function of xi -- is it P_gg?"""
import numpy as np, sys
sys.path.insert(0,'/home/user/correlations/zint/src')
from row3_uv import integrand
x = np.array([0.3, -0.2]); y = np.array([-0.6, 0.9]); Vv = np.array([0.4, 0.7])
th = np.linspace(0, 2*np.pi, 4096, endpoint=False)
u = np.stack([np.cos(th), np.sin(th)], -1)
rho = 1e-5
def coef(xi):
    return np.mean([integrand(x, y, x+uu*rho, Vv, xi)[0] for uu in u])*rho**2
Cuv = lambda xi: xi*(1-xi) + xi/(1-xi) + (1-xi)/xi
print("A. the xi-structure.  The bracket carries 1, 1/xib^2, 1/xi, 1/xi^2, 1/xib, 1/(xi xib);")
print("   multiplied by the overall xi*xib these become")
print("      xi xib ,  xi/xib ,  xib ,  xib/xi ,  xi ,  1")
print("   i.e. exactly the ingredients of C_UV = xi xib + xi/xib + xib/xi, plus xi + xib + 1.")
print()
print("B. measured UV coefficient vs C_UV(xi)")
print(f"   {'xi':>7} {'C(xi) measured':>18} {'C_UV(xi)':>14} {'ratio':>12}")
for xi in (0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95):
    c = coef(xi); r = c/Cuv(xi)
    print(f"   {xi:7.2f} {c:18.8f} {Cuv(xi):14.6f} {r:12.6f}")
print()
print("C. endpoint behaviour (this is what decides whether it is DGLAP)")
for xi in (1e-2, 1e-3, 1e-4):
    print(f"   xi={xi:.0e}:  C = {coef(xi):14.4f}   C_UV = {Cuv(xi):14.4f}   ratio {coef(xi)/Cuv(xi):10.6f}")
for xb in (1e-2, 1e-3, 1e-4):
    xi = 1-xb
    print(f"   xib={xb:.0e}: C = {coef(xi):14.4f}   C_UV = {Cuv(xi):14.4f}   ratio {coef(xi)/Cuv(xi):10.6f}")
