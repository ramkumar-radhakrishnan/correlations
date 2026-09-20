"""Row 20: what the N_c collapse needs -- consistency of the two U(z) index orders."""
import numpy as np, itertools
from scipy.linalg import expm
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3): f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(9)
def Uadj(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Ux,Uxp,Uz=Uadj(rng.normal(size=8)),Uadj(rng.normal(size=8)),Uadj(rng.normal(size=8)*1.2)
print("="*78); print("A.  THE z-SIDE DECIDES WHETHER YOU GET N_c AT ALL"); print("="*78)
# C22^{e'e} = sum f^{ac'd'} [bra U(z)]^{bc'} U^{e'd'}(x') f^{acd} [ket U(z)]^{bc} U^{de}(x)
def C22(braUz,ketUz,braUxp=None):
    bx = Uxp if braUxp is None else braUxp
    return np.einsum('axy,azw,bx,fy,bz,we->fe',f,f,braUz,bx,ketUz,Ux)
cases=[("both U(z)          (as your six parts have it)",Uz,Uz),
       ("both U-dagger(z)                              ",Uz.T,Uz.T),
       ("MIXED: bra U(z), ket U-dagger(z)  (top line)  ",Uz,Uz.T)]
for nm,b_,k_ in cases:
    C=C22(b_,k_); ref=Nc*(Uxp@Ux)
    print("   %s :  ||C22 - N_c U(x')U(x)|| = %10.3e   %s"%(nm,np.abs(C-ref).max(),
          "COLLAPSES" if np.abs(C-ref).max()<1e-10 else "DOES NOT COLLAPSE"))
print()
print("   sum_b U^{bc'}U^{bc} = delta^{c'c}  and  sum_b U^{c'b}U^{cb} = delta^{c'c},")
print("   but  sum_b U^{bc'}U^{cb} = [U U]^{cc'}  is NOT a delta.  So the collapse needs")
print("   the SAME index order on both sides -- either both U(z) or both U^dag(z).")
print("   Your top line writes  U^{dag bc}(z)  in the ket and  U^{bc'}(z)  in the bra;")
print("   your expanded parts silently use U^{bc}(z) for both.  Make the top line match.")
print()
print("="*78); print("B.  THE x-SIDE DECIDES WHICH N_c STRUCTURE YOU GET"); print("="*78)
Cas=np.einsum('axy,azw,bx,fy,bz,we->fe',f,f,Uz,Uxp,Uz,Ux)       # bra U^{e'd'}(x')
Ccj=np.einsum('axy,azw,bx,yf,bz,we->fe',f,f,Uz,Uxp,Uz,Ux)       # bra U^{d'e'}(x')
print("   bra U^{e'd'}(x')  (as written) -> N_c [U(x')U(x)]      : %.2e"%np.abs(Cas-Nc*(Uxp@Ux)).max())
print("   bra U^{d'e'}(x')  (conjugate)  -> N_c [U^dag(x')U(x)]  : %.2e"%np.abs(Ccj-Nc*(Uxp.T@Ux)).max())
print("   difference between the two structures: %.4f"%np.abs(Nc*(Uxp@Ux)-Nc*(Uxp.T@Ux)).max())
print("   The second is the adjoint DIPOLE and is what you want physically.")
