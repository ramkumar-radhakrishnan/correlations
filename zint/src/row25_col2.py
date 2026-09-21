"""Row 25: the colour structure with uniform index orders."""
import numpy as np, itertools
from scipy.linalg import expm
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3): f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(23)
def U(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Uy,Ux,Uxp,Uz,Uw=[U(rng.normal(size=8)) for _ in range(5)]
print("="*80); print("ADOPTING THE UNIFORM ORDERING (f index first everywhere)"); print("="*80)
print("   bra = f^{c'd'a} U^{bd'}(z) [ U^{c'e'}(y') - U^{c'e'}(x') ] rho^{e'}(x')")
print("   ket = f^{cbd} U^{ac}(w) rho^d(x)  -  f^{acd} U^{bc}(z) U^{de}(x) rho^e(x)")
print("   (your bra term 1 has U^{e'c'}(y') and U^{d'b}(z) -- both transposed relative to term 2)")
print()
# the four products, indices [e' (or d), e (or d)] -- keep them as separate objects
b1k1=np.einsum('xya,by,xe,cbd,ac->ed',f,Uz,Uy,f,Uw)
b1k2=-np.einsum('xya,by,xe,acd,bc,dg->eg',f,Uz,Uy,f,Uz,Ux)
b2k1=-np.einsum('xya,by,xe,cbd,ac->ed',f,Uz,Uxp,f,Uw)
b2k2=np.einsum('xya,by,xe,acd,bc,dg->eg',f,Uz,Uxp,f,Uz,Ux)
print("   %-12s %12s   %s"%("product","||.||_max","collapses?"))
for nm,M in (("bra1 x ket1",b1k1),("bra1 x ket2",b1k2),("bra2 x ket1",b2k1),("bra2 x ket2",b2k2)):
    cands={"N_c U^dag(x')U(x)":Nc*(Uxp.T@Ux),"-N_c U^dag(x')U(x)":-Nc*(Uxp.T@Ux),
           "N_c I":Nc*np.eye(8),"-N_c I":-Nc*np.eye(8)}
    best=min(cands.items(),key=lambda kv: np.abs(M-kv[1]).max())
    ok = np.abs(M-best[1]).max()<1e-10
    print("   %-12s %12.4f   %s"%(nm,np.abs(M).max(),
          ("YES -> %s"%best[0]) if ok else "no (closest %s, residual %.3f)"%(best[0],np.abs(M-best[1]).max())))
print()
print("   Only bra2 x ket2 collapses, and it gives  -N_c [U^dag(x')U(x)]^{e'e} :  residual %.2e"
      %np.abs(b2k2+Nc*(Uxp.T@Ux)).max())
print()
print("   WHY the others do not: the collapse needs sum_b U^{b.}U^{b.} = delta, which requires")
print("   BOTH sides to carry a Wilson line at the SAME point on the SAME free index.  Here")
print("     bra: a sits on f, b sits on U(z)      -- always")
print("     ket term 1: a sits on U(w), b sits on f      <- no U(z) at all, so no delta in b")
print("     ket term 2: a sits on f,    b sits on U(z)   <- matches the bra, and only then does it fire")
print("   So unlike the two-C row, this interference does NOT reduce to a single N_c structure.")
print("   Three of the four terms keep their full Wilson-line content and must be carried along.")
