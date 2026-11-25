#!/usr/bin/env python3
"""dilute2 / e6: grid check of Row I.III (v3c_fast.R_III_pu) at one point and one Lambda with a different radial resolution
and a different partition-of-unity radius rho (a Lambda-independent quadrature error would not be seen in the fit in L).
usage: e6_rowIII_grid_check.py lx ly Lam nr nper rho"""
import sys, time, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import v3c_fast as f
k = np.array([1.0, 0.0]); l = np.array([float(sys.argv[1]), float(sys.argv[2])]); Lam = float(sys.argv[3])
nr, nper, rho = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
t = time.time(); R = f.R_III_pu(k, l, Lam, nr=nr, nper=nper, rho=rho)
print(f"l={l} Lam={Lam:.0e} nr={nr} nper={nper} rho={rho}: R_III(vee=1)={R[0]:.8f}  coeff ln vee={R[1]:.8f}  ({time.time()-t:.0f}s)")
