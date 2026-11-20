import sys; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v5_rowV as v5
k=np.array([float(a) for a in sys.argv[1:3]]); l=np.array([float(a) for a in sys.argv[3:5]])
xs=[float(a) for a in sys.argv[5:]]
for xi in xs:
    fin,pole,tot=v5.rowV_xi(xi,k,l,eps=1e-9)
    print(f"{xi:.8e} {tot:.12e}",flush=True)
