#!/usr/bin/env python3
"""dilute2 / e9: complete target-side Regge integrand of the virtual rows as a function of the momentum P of the soft line:
   tau_0(P)   = lim_{s->0} s [F_III + F_IV + F_V(res)](s, p = P)                 (loop gluon soft),
   tau_fwd(P) = lim_{t->0} t [F_III(1+t, k + sg1*P) + F_III(1-t, k + sg2*P)]      (third line of Row I.III soft, |p+ - k+| = t k+).
Question: is tau_0 + tau_fwd independent of k at fixed (P, l)?  (e7d: tau_0 alone is not.)"""
import sys, os, io, contextlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
with contextlib.redirect_stdout(io.StringIO()):
    import e8_omega as e8
import v3c_fast as v3
def tau_fwd(P, k, l, sg1, sg2, t=1e-7):
    a = t * v3.Fvec(np.array([1 + t]), k + sg1 * P, k, l)[0]; b = t * v3.Fvec(np.array([1 - t]), k + sg2 * P, k, l)[0]
    return a, b
l = np.array([0.3, 0.7]); ks = [np.array(x) for x in ((1.0, 0.0), (0.6, 0.5), (-0.8, 1.7), (2.5, -0.3))]
for P in [np.array(x) for x in ((0.4, 0.5), (-0.3, -0.2), (1.5, -0.7), (3.0, 4.0))]:
    print(f"P=({P[0]:5.2f},{P[1]:5.2f})")
    for k in ks:
        t0 = e8.tau3(P, k, l)
        f = {(s1, s2): tau_fwd(P, k, l, s1, s2) for s1 in (1, -1) for s2 in (1, -1)}
        print(f"   k=({k[0]:5.2f},{k[1]:5.2f}): tau_0 = {t0: .6f} | fwd s>1: +P {f[(1,1)][0]: .6f}  -P {f[(-1,1)][0]: .6f} | fwd s<1: +P {f[(1,1)][1]: .6f}  -P {f[(1,-1)][1]: .6f}"
              f" | sums: {t0 + f[(1,1)][0] + f[(1,1)][1]: .6f} {t0 + f[(1,1)][0] + f[(1,-1)][1]: .6f} {t0 + f[(-1,1)][0] + f[(1,1)][1]: .6f} {t0 + f[(-1,1)][0] + f[(1,-1)][1]: .6f}")
