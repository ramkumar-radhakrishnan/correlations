"""Row 22: the Wilson-line reduction of the two-C row."""
import numpy as np, itertools
from scipy.linalg import expm
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3): f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(13)
def Uadj(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Uy,Uyp,Ux,Uxp,Uz=[Uadj(rng.normal(size=8)) for _ in range(5)]
print("="*80); print("1.  THE REDUCTION, WITH CONSISTENT INDEX ORDERS"); print("="*80)
print("   ket = f^{cda} U^{db}(z) [ U^{ec}(y) - U^{ec}(x) ] rho^e(x)")
print("   bra = f^{c'd'a} U^{d'b}(z) [ U^{e'c'}(y') - U^{e'c'}(x') ] rho^{e'}(x')")
print("     sum_b U^{d'b}(z)U^{db}(z) = delta^{d'd} ;  sum_{a,d} f^{c'da}f^{cda} = N_c delta^{c'c}")
Dy,Dyp=Uy-Ux,Uyp-Uxp
good=np.einsum('xya,zwa,yb,wb,ex,fz->fe',f,f,Uz,Uz,Dy,Dyp)   # [e',e]
ref =Nc*np.einsum('fc,ec->fe',Dyp,Dy)
print("   ||result - N_c [U(y')-U(x')][U(y)-U(x)]|| = %.2e    CORRECT"%np.abs(good-ref).max())

print(); print("="*80); print("2.  BUT YOUR INDEX ORDERS ARE NOT CONSISTENT"); print("="*80)
print("   you wrote   ket : f^{cda} U^{ec}(y) U^{db}(z)  -  f^{cda} U^{bd}(z) U^{ce}(x)")
print("                                    ^^^^^^^^                ^^^^^^^^")
print("               bra : f^{c'd'a} U^{e'c'}(y') U^{d'b}(z)  -  f^{c'd'a} U^{bd'}(z) U^{c'e'}(x')")
print("   two problems:")
print("     (i)  U^{db}(z) in term 1 vs U^{bd}(z) in term 2  -- the z lines are transposed")
print("          relative to each other, so sum_b does NOT give a delta in the cross terms:")
print("            sum_b U^{d'b} U^{bd} = [U(z)U(z)]^{d'd}  (not delta)      ||.-1|| = %.3f"
      %np.abs(Uz@Uz-np.eye(8)).max())
print("     (ii) U^{ec}(y) vs U^{ce}(x) -- transposed, so the DIFFERENCE structure breaks:")
print("            [U(y) - U(x)^T]  does not vanish at y = x :  ||U(x)-U(x)^T|| = %.3f"%np.abs(Ux-Ux.T).max())
print()
bad=( np.einsum('xya,zwa,ex,yb,fz,wb->fe',f,f,Uy,Uz,Uyp,Uz)       # (1)(1'): U^{ec}(y)U^{db}(z) x ...
     -np.einsum('xya,zwa,ex,yb,fz,bw->fe',f,f,Uy,Uz,Uxp,Uz)
     -np.einsum('xya,zwa,xe,by,fz,wb->fe',f,f,Ux,Uz,Uyp,Uz)
     +np.einsum('xya,zwa,xe,by,fz,bw->fe',f,f,Ux,Uz,Uxp,Uz) )
print("   literally as you wrote it :  ||result - N_c[U(y')-U(x')][U(y)-U(x)]|| = %.4f"%np.abs(bad-ref).max())
print("   (i.e. the quoted N_c form is right, but only after the orderings are made uniform)")
print()
print("   also: your final line reads  [U^{ec}(y) - U^{ce}(x)] -- the second should be U^{ec}(x),")
print("   otherwise the bracket does not vanish at y = x and the collinear taming is lost.")
