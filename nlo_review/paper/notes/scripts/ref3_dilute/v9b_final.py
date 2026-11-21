#!/usr/bin/env python3
"""ref3_dilute / v9b: final assembly of Group I (point quark, dilute) at two (k,l) points, MS-bar mu = |k| = 1, units g^2 = 1.
R_I = cV * ln(vee/k+) + cL * L + c2 * L^2 + (2/eps)(pV * ln vee + pL * L) + C,  L = ln(k+/Lambda).
Inputs (copied run outputs): v4b_sgrid_pt*.out, v5_endpoints_pt*.out, v5_nodes_pt*.out, h2 scans, A_inf integrals."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
import sys
S = sys.argv[1]
Nc = 3; c8 = Nc / (8 * np.pi ** 2)
pts = {
 'pt1': dict(l=np.array([0.3, 0.7]), III=[(1e-5, -2.46802238), (1e-6, -3.57719450), (1e-7, -4.89472033), (1e-8, -6.42073748)], Ainf=-0.04492649),
 'pt2': dict(l=np.array([-0.5, 0.4]), III=[(1e-5, 5.60561281), (1e-6, 8.19570618), (1e-7, 11.27714490), (1e-8, 14.84983799)], Ainf=-0.00743644),
}
res = {}
for name, P in pts.items():
    k = np.array([1.0, 0.0]); l = P['l']; q = k + l
    v0 = k / (k @ k) - q / (q @ q); vv = v0 @ v0; Pk = (v0 @ k) / vv; Pq = (v0 @ q) / (q @ q) / vv
    out = {}
    out['I.II comm'] = dict(c2=Nc / (4 * np.pi ** 2), cL=(Nc / (4 * np.pi ** 2)) * Pq * np.log(q @ q), cV=0, C=Nc / 24, pL=Nc / (4 * np.pi ** 2), pV=0)
    d = np.loadtxt(f'{S}/v4b_sgrid_{name}.out'); u, sG = d[:20, 0], d[:20, 1]; g0, gi = d[20, 1], d[21, 1]
    cs = CubicSpline(u, sG); I_IV = cs.integrate(-12, 0) - 12 * g0 + cs.integrate(0, 12) - 12 * gi
    out['I.IV'] = dict(c2=0, cL=g0, cV=gi, C=I_IV, pL=d[0, 2], pV=d[0, 2])
    e = np.loadtxt(f'{S}/v5_endpoints_{name}.out'); xe, te = e[:, 0], e[:, 1]
    f = lambda x, A, B, C1, D1: A + B * np.log(x) + C1 * x * np.log(x) + D1 * x
    (A, B, _, _), _ = curve_fit(f, xe, xe * te, p0=[-0.01, -0.02, 0, 0])
    b = np.loadtxt(f'{S}/v5_nodes_{name}.out'); xb, tb = b[:, 0], b[:, 1]
    xg, wg = np.polynomial.legendre.leggauss(len(xb)); tg = 0.5 * (xg + 1)
    Ifin = np.sum(wg * 0.5 * tg * (tb - (A + B * np.log(xb)) / xb))
    out['I.V resolved'] = dict(c2=-B, cL=2 * A, cV=0, C=2 * (Ifin - A * np.log(2) + 0.5 * B * np.log(2) ** 2), pL=0, pV=0)
    La = np.array([x for x, _ in P['III']]); R = np.array([y for _, y in P['III']]); Ls = np.log(1 / La)
    M = np.vstack([Ls ** 2, Ls, np.ones_like(Ls)]).T; (a2, a1, a0), *_ = np.linalg.lstsq(M, R, rcond=None)
    out['I.III'] = dict(c2=a2, cL=a1, cV=P['Ainf'], C=a0, pL=0, pV=0)
    tot = {kk: sum(o[kk] for o in out.values()) for kk in ('c2', 'cL', 'cV', 'C', 'pL', 'pV')}
    print(f"== {name}: k={k}, l={l}, q={q}; P_k={Pk:.5f}, P_q={Pq:.5f}")
    for r, o in out.items():
        print(f"   {r:13s}: L^2 {o['c2']: .6f}  L {o['cL']: .6f}  V {o['cV']: .6f}  C {o['C']: .6f}  poles (2/eps)[{o['pL']:.6f} L + {o['pV']:.6f} V]")
    print(f"   TOTAL        : L^2 {tot['c2']: .6f}  L {tot['cL']: .6f}  V {tot['cV']: .6f}  C {tot['C']: .6f}  poles (2/eps)[{tot['pL']:.6f} L + {tot['pV']:.6f} V]")
    print(f"   expectations : L^2 Nc/(4pi^2)={Nc/(4*np.pi**2):.6f};  L: 2omega(q_b) finite={-Nc/(4*np.pi**2)*np.log(l@l):.6f};  V: 2omega(q_a) finite={-Nc/(4*np.pi**2)*np.log(q@q):.6f};"
          f"  pole of 2omega: {Nc/(4*np.pi**2):.6f}")
    print(f"   I.III L^2 fit {a2:.6f} vs -(Nc/8pi^2)P_q = {-c8*Pq:.6f};  resolved I.V L^2 {-B:.6f} vs +(Nc/8pi^2)P_q = {c8*Pq:.6f}")
    # Fadin-Lipatov virtual kernel (BSV eq. kernelv), finite part at mu = |Delta| = |k| = 1, n_f = 0, units g^2: gbar^2 = Nc/(16 pi^2)
    qa = q; qb = -l; D2 = k @ k; a2_, b2_ = qa @ qa, qb @ qb; ab = qa @ qb; dqq = (qa - qb) @ (qa - qb)
    def X(a2_, b2_):
        Lab = np.log(a2_ / b2_); den = a2_ - b2_
        br = D2 / den * (1 - D2 * (a2_ + b2_ - 4 * ab) / (3 * den ** 2)) * Lab - D2 * dqq / (6 * a2_ * b2_) + D2 ** 2 * (a2_ + b2_) * (a2_ + b2_ - 4 * ab) / (6 * a2_ * b2_ * den ** 2)
        return Nc / (16 * np.pi ** 2) * (np.pi ** 2 + 3 * D2 / den * Lab - Lab ** 2 + br)
    if abs(a2_ - b2_) < 1e-6:
        xs = [X(a2_ * (1 + h), b2_) for h in (1e-3, 2e-3, -1e-3, -2e-3)]
        Xfl = (4 * (xs[0] + xs[2]) / 2 - (xs[1] + xs[3]) / 2) / 3          # Richardson limit h -> 0
    else:
        Xfl = X(a2_, b2_)
    print(f"   FL virtual kernel K^(v)/K^B, finite part at mu=|k| (for reference, NOT regulator-matched): {Xfl:.6f};  CGC Group I C: {tot['C']:.6f}")
    res[name] = tot
