#!/usr/bin/env python3
"""ref3_dilute / v5c: resolved remainder of Row I.V (= R_V + R_SE, the subtracted 2D integral of v5) at fixed xi
as the target momentum l -> 0 (pair unresolved): does it vanish (no LO-type constant) or tend to a constant?"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v5_rowV as v5
k = np.array([1.0, 0.0]); e = np.array([0.6, 0.8])
for xi in (0.3,):
    for L in (0.7, 0.2, 0.05, 0.0125):
        fin, pole, tot = v5.rowV_xi(xi, k, L * e)
        print(f"xi={xi} |l|={L}: resolved remainder = {tot:.8f}", flush=True)
