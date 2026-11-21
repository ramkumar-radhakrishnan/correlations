#!/usr/bin/env python3
"""
ref3_dilute / v4b: Row I.IV integrated (point quark, dilute), DR (MS-bar) in p, cut-offs in s:
  G(s) = int d^dp/(2pi)^d F_IV(s,p)
       = int d^2p/(2pi)^2 [F - c_IR(s) th(p<lam)/p^2 - c_UV(s) th(p>lam)/p^2]
         + (c_UV - c_IR)/(4pi) [2/eps - ln(lam^2/mu^2)],        F ~ c_IR/p^2 (p->0), c_UV/p^2 (p->inf)
  R_IV = int_Lambda^vee ds G(s).  c_IR, c_UV are proportional to 1/s (checked), so the pole multiplies ln(vee/Lambda).
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v4_rowIV as v4

def cIRUV(s, k, l):
    th = np.linspace(0, 2 * np.pi, 33)[:-1]
    cir = np.mean([v4.Rdens(s, 1e-7 * np.array([np.cos(t), np.sin(t)]), k, l) * 1e-14 for t in th])
    cuv = np.mean([v4.Rdens(s, 1e7 * np.array([np.cos(t), np.sin(t)]), k, l) * 1e14 for t in th])
    return cir, cuv

def G(s, k, l, lam=1.0, mu=1.0, nr=24, nth=48):
    cir, cuv = cIRUV(s, k, l)
    import v3c_fast
    # point singularities (integrable 1/|p-p0|) at p0 = -l, -k-l, -k: geometric angular refinement towards them
    sp = [np.arctan2(-l[1], -l[0]), np.arctan2(-(k + l)[1], -(k + l)[0]), np.arctan2(-k[1], -k[0])]
    th, wth = v3c_fast.angular_nodes(sp, max(4, nth // 6))
    xr, wr = np.polynomial.legendre.leggauss(nr)
    base = [np.linalg.norm(k), np.linalg.norm(l), np.linalg.norm(k + l)]
    sc = sorted(set(base + [lam] + [np.sqrt(s) * b for b in base] + [np.sqrt(s) * b * f_ for b in base for f_ in (0.3, 3.0)]))
    edges = [1e-6 * sc[0]] + sc + [1e6 * sc[-1]]
    tot = 0.0
    for r0, r1 in zip(edges[:-1], edges[1:]):
        u0, u1 = np.log(r0), np.log(r1)
        for x, w in zip(xr, wr):
            u = 0.5 * (u1 - u0) * x + 0.5 * (u1 + u0); r = np.exp(u)
            sub = (cir if r < lam else cuv) / (r * r)
            acc = 0.0
            for t, wt in zip(th, wth):
                acc += wt * (v4.Rdens(s, np.array([r * np.cos(t), r * np.sin(t)]), k, l) - sub)
            tot += acc * w * 0.5 * (u1 - u0) * r * r
    num = tot / (2 * np.pi) ** 2
    # tails beyond edges: |F - c/p^2| decays faster; neglected (checked by moving edges)
    fin = num - (cuv - cir) / (4 * np.pi) * np.log(lam ** 2 / mu ** 2)
    pole = (cuv - cir) / (4 * np.pi)            # coefficient of 2/eps
    return fin, pole, cir, cuv

if __name__ == '__main__':
    k = np.array([1.0, 0.0]); l = np.array([0.3, 0.7])
    for s in (1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0, 1e3, 1e4):
        fin, pole, cir, cuv = G(s, k, l)
        print(f"s={s:g}: G_fin={fin:.8f}  s*G_fin={s*fin:.8f}  pole={pole:.6e}  s*c_IR={s*cir:.6f} s*c_UV={s*cuv:.6f}", flush=True)
