import sys; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v5_rowV as v5
k=np.array([float(a) for a in sys.argv[1:3]]); l=np.array([float(a) for a in sys.argv[3:5]])
x,w=np.polynomial.legendre.leggauss(12)
t=0.5*(x+1); xis=0.5*t**2     # xi = t^2/2 on [0,1/2]
for xi in xis:
    fin,pole,tot=v5.rowV_xi(xi,k,l,eps=1e-9)
    print(f"{xi:.14e} {tot:.12e}",flush=True)
