import sys; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v4b_integrate as v
k=np.array([float(a) for a in sys.argv[1:3]]); l=np.array([float(a) for a in sys.argv[3:5]])
# u = ln s nodes on [-12, 12] (GL 40) plus far points for the tails
x,w=np.polynomial.legendre.leggauss(40)
us=12*x
for u in list(us)+[-18.0,-15.0,15.0,18.0]:
    s=np.exp(u); fin,pole,cir,cuv=v.G(s,k,l)
    print(f"{u:.12f} {s*fin:.12e} {s*pole:.12e}",flush=True)
