import sys,time; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v3c_fast as f
k=np.array([float(a) for a in sys.argv[1:3]]); l=np.array([float(a) for a in sys.argv[3:5]])
nr=int(sys.argv[5]); nper=int(sys.argv[6]); lams=[float(a) for a in sys.argv[7:]]
for Lam in lams:
    t=time.time(); R=f.R_III_pu(k,l,Lam,nr=nr,nper=nper)
    print(f"Lam={Lam:.0e}  R_III(vee=1)={R[0]:.8f}  coeff ln vee={R[1]:.8f}  ({time.time()-t:.0f}s)",flush=True)
