#!/usr/bin/env python3
"""
ref3_dilute / v9: assemble Group I (point quark, dilute) at one (k,l):
   R_I = R_IIcomm + R_Vres + R_III + R_IV   (I.II self-energy + unresolved I.V = 0 exactly, v5b)
in units g^2 = 1 (multiply by 4 pi alpha_s), MS-bar, mu = |k|, L = ln(k+/Lambda), V = ln(vee/k+).
Inputs: run outputs r3_*.out (Row I.III Lambda scan), r4.out (Row I.IV s-grid), r5_end.out + r5b.out (resolved I.V).
"""
import sys, numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
D = sys.argv[1]; r3file = sys.argv[2]
k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7]); q = k + l; Nc = 3
v0 = k / (k @ k) - q / (q @ q); vv = v0 @ v0
mu2 = k @ k
out = {}
# --- Row I.II, commutator part (analytic) ---
coefL = -(Nc / (4 * np.pi ** 2)) * ((v0 @ k) / (k @ k) * np.log((k @ k) / mu2) - (v0 @ q) / (q @ q) * np.log((q @ q) / mu2)) / vv
out['IIcomm'] = dict(L2=Nc / (4 * np.pi ** 2), L=coefL, V=0.0, c=Nc / 24, pole_L=Nc / (4 * np.pi ** 2), pole_V=0.0)
# --- Row I.IV ---
d = np.loadtxt(f'{D}/r4.out'); u = d[:40, 0]; sG = d[:40, 1]
g0, ginf = d[40, 1], d[43, 1]
csL = CubicSpline(u[u < 0], sG[u < 0] - g0); csR = CubicSpline(u[u > 0], sG[u > 0] - ginf)
uu = np.linspace(u[0], 0, 20001); vals = np.interp(uu, u, sG) - g0
# spline on each side, extended to 0 by the spline
IL = CubicSpline(u, sG - g0).integrate(-12, 0) if False else None
cs = CubicSpline(u, sG)
I_IV = cs.integrate(-12, 0) - 12 * g0 + cs.integrate(0, 12) - 12 * ginf
out['IV'] = dict(L2=0.0, L=g0, V=ginf, c=I_IV, pole_L=d[0, 2], pole_V=d[0, 2])
# --- resolved Row I.V ---
e = np.loadtxt(f'{D}/r5_end.out'); xe, te = e[:, 0], e[:, 1]
xs_all = np.concatenate([xe, [1e-3, 1e-2]]); ts_all = np.concatenate([te, [123.7017019657, 7.7711490621]])
f = lambda x, A, B, C, Dd: A + B * np.log(x) + C * x * np.log(x) + Dd * x
popt, _ = curve_fit(f, xs_all, xs_all * ts_all, p0=[-0.01, -0.02, 0, 0])
A, B = popt[0], popt[1]
b = np.loadtxt(f'{D}/r5b.out'); xb, tb = b[:, 0], b[:, 1]
xg, wg = np.polynomial.legendre.leggauss(len(xb)); tg = 0.5 * (xg + 1); jac = 0.5 * tg * 0.5 * 2  # xi = t^2/2, dxi = t dt, dt = dx/2
sub = tb - (A + B * np.log(xb)) / xb
Ifin = np.sum(wg * 0.5 * tg * sub)                    # int_0^{1/2} dxi [..] = int_0^1 dt t [..]
out['Vres'] = dict(L2=-B, L=2 * A, V=0.0, c=2 * (Ifin - A * np.log(2) + B * 0.5 * np.log(2) ** 2), pole_L=0.0, pole_V=0.0,
                   A=A, B=B, fitpars=popt)
# --- Row I.III ---
lines = [ln for ln in open(f'{D}/{r3file}') if ln.startswith('Lam=')]
Lam = np.array([float(ln.split()[0][4:]) for ln in lines]); Rv = np.array([float(ln.split('=')[2].split()[0]) for ln in lines])
cvee = np.mean([float(ln.split('=')[3].split()[0]) for ln in lines])
Ls = np.log(1 / Lam)
M = np.vstack([Ls ** 2, Ls, np.ones_like(Ls)]).T
c3, *_ = np.linalg.lstsq(M, Rv, rcond=None)
out['III'] = dict(L2=c3[0], L=c3[1], V=cvee, c=c3[2], pole_L=0.0, pole_V=0.0, resid=np.max(np.abs(M @ c3 - Rv)))
tot = {key: sum(out[r][key] for r in ('IIcomm', 'IV', 'Vres', 'III')) for key in ('L2', 'L', 'V', 'c', 'pole_L', 'pole_V')}
for r in out:
    print(r, {kk: (round(vv_, 6) if np.isscalar(vv_) else vv_) for kk, vv_ in out[r].items()})
print("TOTAL Group I (units g^2): R = %.6f L^2 + %.6f L + %.6f ln(vee/k+) + %.6f  + (2/eps)[%.6f L + %.6f ln vee]" %
      (tot['L2'], tot['L'], tot['V'], tot['c'], tot['pole_L'], tot['pole_V']))
