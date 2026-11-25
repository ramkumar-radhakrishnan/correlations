#!/usr/bin/env python3
"""dilute2 / e5: refinement checks of the cheap rows at one point: Row I.IV on a dense grid in u = ln s (step 0.5, and
higher transverse resolution at a few u), Row I.V resolved with 24 instead of 12 nodes and two more end-point values.
usage: e5_refine_IV_V.py lx ly tag"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v4b_integrate as v4, v5_rowV as v5
k = np.array([1.0, 0.0]); l = np.array([float(sys.argv[1]), float(sys.argv[2])]); tag = sys.argv[3]
with open(f'e5_V24_{tag}.out', 'w') as f:
    x, w = np.polynomial.legendre.leggauss(24); t = 0.5 * (x + 1)
    for xi in list(0.5 * t ** 2) + [3e-6, 1e-3]:
        fin, pole, tot = v5.rowV_xi(xi, k, l, eps=1e-9); f.write(f"{xi:.14e} {tot:.12e}\n"); f.flush()
with open(f'e5_IVdense_{tag}.out', 'w') as f:
    for u in np.arange(-12, 12.01, 0.5):
        s = np.exp(u); fin, pole, cir, cuv = v4.G(s, k, l, nr=16, nth=48); f.write(f"{u:.6f} {s * fin:.12e}\n"); f.flush()
    for u in (-3.0, 0.0, 2.0):
        s = np.exp(u); fin, pole, cir, cuv = v4.G(s, k, l, nr=32, nth=96); f.write(f"# hires {u:.6f} {s * fin:.12e}\n"); f.flush()
