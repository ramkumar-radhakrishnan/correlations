import sys; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v4b_integrate as v
k=np.array([float(a) for a in sys.argv[1:3]]); l=np.array([float(a) for a in sys.argv[3:5]])
x,w=np.polynomial.legendre.leggauss(20)
for u in list(12*x)+[-20.0,20.0]:
    s=np.exp(u); fin,pole,cir,cuv=v.G(s,k,l,nr=16,nth=48)
    print(f"{u:.12f} {s*fin:.12e} {s*pole:.12e}",flush=True)
