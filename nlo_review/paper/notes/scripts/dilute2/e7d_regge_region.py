#!/usr/bin/env python3
"""dilute2 / e7d: Regge region of the virtual rows, s = p+/k+ -> 0 at fixed p:
   tau(p; k, l) = lim s [F_III + F_IV + F_V(resolved)](s,p)    (ratio to LO, per d^2p/(2pi)^2 ds/s).
Regge factorisation of the target side (and the conversion of the cut p+ > Lambda into a k-independent redefinition of
phi_T) requires tau to depend on (p, l) only, not on k.  Checked by changing k at fixed l."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e7c_region_b as eb
def tau(p, k, l, s):
    return sum(eb.Hb(p / s, k, l, s)) / (p @ p)      # Hb multiplies by s p^2 with p = s v
if __name__ == '__main__':
    l = np.array([0.3, 0.7])
    ks = [np.array(x) for x in ((1.0, 0.0), (0.6, 0.5), (-0.8, 1.7), (2.5, -0.3))]
    print("l =", l, ";  columns: p ; tau for k =", [tuple(k) for k in ks])
    for p in [np.array(x) for x in ((0.4, 0.5), (-0.3, -0.2), (1.5, -0.7), (0.05, 0.02), (-0.31, -0.69), (3.0, 4.0), (30.0, 10.0))]:
        print(f"p=({p[0]:6.2f},{p[1]:6.2f}) " + "  ".join(f"{tau(p, k, l, 1e-9): .8f}" for k in ks) + f"   | p^2 tau: " + " ".join(f"{(p@p)*tau(p, k, l, 1e-9): .5f}" for k in ks[:2]))
