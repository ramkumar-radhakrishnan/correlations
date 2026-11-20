import sys,time; sys.path.insert(0,'/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import numpy as np, v3c_fast as f
k=np.array([float(a) for a in sys.argv[1:3]]); l=np.array([float(a) for a in sys.argv[3:5]])
nr=int(sys.argv[5]) if len(sys.argv)>5 else 96
Ls=[]; Rs=[]
for Lam in (1e-3,1e-4,1e-5,1e-6,1e-7):
    t=time.time(); R=f.R_III(k,l,Lam,nr=nr,nth=int(sys.argv[6]) if len(sys.argv)>6 else 8)
    Ls.append(np.log(1/Lam)); Rs.append(R[0])
    print(f"Lam={Lam:.0e}  R_III(vee=1)={R[0]:.8f}  coeff ln vee={R[1]:.8f}  ({time.time()-t:.0f}s)",flush=True)
A=np.vstack([np.array(Ls)**2,Ls,np.ones(len(Ls))]).T
c,res,_,_=np.linalg.lstsq(A,np.array(Rs),rcond=None)
print("fit R_III = %.8f L^2 + %.8f L + %.8f ; max resid %.2e"%(c[0],c[1],c[2],np.max(np.abs(A@c-Rs))))
