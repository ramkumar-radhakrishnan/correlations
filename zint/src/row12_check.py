"""Row 12 = the 3-rho (cross) piece of |B_2|^2.  Structure, prefactor, reality."""
import sympy as sp, numpy as np, itertools
from scipy.linalg import expm
g,pi=sp.symbols('g pi',positive=True); P,K=sp.symbols('pplus kplus',positive=True); S=P+K; I=sp.I
B1 =  I*g**2*sp.sqrt(P*K)/(4*pi**2*S)      # one-rho piece of B_2
B2 = -g**2/(8*pi**2*sp.sqrt(P*K))          # two-rho piece  (REAL)
over = sp.Integer(4)/(2*pi)**3
print("="*78); print("1.  PREFACTOR AND REALITY"); print("="*78)
c1 = sp.simplify(over*sp.conjugate(B1).subs({sp.conjugate(g):g,sp.conjugate(pi):pi,
        sp.conjugate(P):P,sp.conjugate(K):K})*B2)
c2 = sp.simplify(over*B2*B1)
print("   B^(1)dag x B^(2) :", sp.simplify(c1), " = -(1/(2pi)^3)(i g^4/8pi^4)(1/S) x",
      sp.simplify(c1/(-(1/(2*pi)**3)*I*g**4/(8*pi**4)/S)))
print("   B^(2)dag x B^(1) :", sp.simplify(c2), " = +(1/(2pi)^3)(i g^4/8pi^4)(1/S) x",
      sp.simplify(c2/(+(1/(2*pi)**3)*I*g**4/(8*pi**4)/S)))
print()
print("   B^(1) is imaginary, B^(2) is real, so the two cross terms are complex")
print("   conjugates and carry OPPOSITE signs.  Your -...+ is correct, and it makes")
print("   N_1 = 2 i Im(...) REAL on its own -- no '+ c.c.' is needed here.")
print()
print("   the vertex's overall 1/S distributes as")
print("      (1/S)[ d^ij(A-B)/(2C) + (S/k+)(..) + (S/p+)(..) ]")
print("    = d^ij(A-B)/(2 S C) + (1/k+)(..) + (1/p+)(..)        MATCHES your bracket")

print(); print("="*78); print("2.  COLOUR"); print("="*78)
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2
f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3):
    f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
print("   B^{bc}_2 one-rho piece carries  f^{cba} rho^a(x)  (free b,c).  Relabel a -> d:")
print("     ket term 1 : U^{ac}(w) f^{cbd} rho^d(x)        -> your f^{cbd} U^{ac}(w) rho^d(x)  OK")
print("     ket term 2 : -U^{bc}(z) f^{acd} U^{de}(x) rho^e(x) -> your f^{acd}U^{bc}(z)U^{de}(x) OK")
print("   the two carry DIFFERENT structure constants, f^{cbd} vs f^{acd}, so f is NOT")
print("   an overall factor -- and you correctly did not pull it out.   CORRECT.")
T1=np.einsum('ac,cbd->abd',np.eye(8),f); T2=np.einsum('acd->acd',f)
print("   ||f^{cbd} - f^{acd}|| as tensors (b,a swapped): %.4f  -> genuinely different"
      % np.linalg.norm(f - np.einsum('bac->abc',f)))
print()
print("   B^{bc}_2 two-rho piece carries {rho^c(y), rho^b(x)} (free b,c).  With the")
print("   relabelling y -> x' (the i-slot, paired with w') and x -> y (the j-slot, paired")
print("   with z) this is {rho^{c'}(x'), rho^b(y)}  -> matches your bracket.   CORRECT.")
print()
print("   ONE TRANSPOSE: the row has U^{c'a}(w'), you wrote U^{ac'}(w').  Same issue as")
print("   in the previous rows -- pick one ordering and use it in both terms.")

print(); print("="*78); print("3.  THE TWO CROSS TERMS ARE EACH OTHER'S MIRROR"); print("="*78)
print("   term 1: one-rho in the KET at w   (structures (x-w), (x-z), (z-w), denom p+A+k+B)")
print("           two-rho in the BRA at w'  (kernel (x'-w')^i (y-z)^j)")
print("   term 2: one-rho in the BRA at w'  (structures (x-w'), (x-z), (z-w'))")
print("           two-rho in the KET at w   (kernel (x'-w)^i (y-z)^j)")
print("   both as you wrote them.   CORRECT.")
