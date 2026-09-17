"""Row 10: the colour structures of N1 and N2."""
import numpy as np, itertools
from scipy.linalg import expm
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2
f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3):
    f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(13)
def U(al):
    H=sum(al[a]*T[a] for a in range(8)); V=expm(1j*H)
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Uw,Ux,Uy,Uxp,Uwp = [U(rng.normal(size=8)) for _ in range(5)]

print("="*76); print("1.  N1's BRACKETS  (one rho each) -- these are fine"); print("="*76)
print("   bra :  -Abar^dag + A^dag U^{b'a}(w')  ->  [ U^{ab'}(x') - U^{b'a}(w') ] rho^{b'}(x')")
print("   you wrote                                 [ U^{ab'}(x') - U^{ab'}(w') ] rho^{b'}(x')")
print("   ||U(w') - U(w')^T|| = %.4f   -> the two differ; pick ONE ordering." % np.abs(Uwp-Uwp.T).max())
print("   ket :  U^{ab}(w) A^(3)b - Abar^(3)a      ->  [ U^{ab}(w) - U^{ab}(x) ] rho^b(x)")
print("   you wrote the same                                                        MATCHES")

print(); print("="*76); print("2.  N2's KET BRACKET IS INDEX-INCONSISTENT"); print("="*76)
print("   A^(3)a (two-rho piece) carries  f^{abc} {rho^b(x), rho^c(y)} , free index a.")
print()
print("   term 1 of the ket is  U^{ab}(w) A^(3)b  =  U^{ab}(w) f^{bcd} {rho^c(x), rho^d(y)}")
print("   you wrote             f^{abc} U^{ab}(w) {rho^b(x), rho^c(y)}")
print("        -> the index b appears THREE times (in f, in U, in rho).  Ill-formed.")
print()
print("   term 2 of the ket is  f^{abc} U^{bd}(x) U^{ce}(y) {rho^d(x), rho^e(y)}   -- well formed,")
print("        and that IS what you wrote.  So f^{abc} can be pulled out of term 2 but NOT term 1.")
print()
T1 = np.einsum('ab,bde->ade', Uw, f)                 # correct term 1
T2 = np.einsum('abc,bd,ce->ade', f, Ux, Uy)          # term 2
print("   are the two colour tensors proportional?  ||T1|| = %.4f , ||T2|| = %.4f ,"
      % (np.linalg.norm(T1), np.linalg.norm(T2)))
print("   ||T1 - T2|| = %.4f   -> genuinely different tensors; f is NOT an overall factor."
      % np.linalg.norm(T1-T2))
same = np.einsum('abc,bd,ce->ade', f, Ux, Ux)
pred = np.einsum('ax,xde->ade', Ux, f)
print("   (only if x = y does term 2 collapse:  f^{abc}U^{bd}(x)U^{ce}(x) = U^{aa'}(x) f^{a'de} ,")
print("    max deviation %.2e )" % np.abs(same-pred).max())

print(); print("="*76); print("3.  OTHER THINGS TO FIX"); print("="*76)
print("   (a) N2's measure reads  int_{w,w'} ... int_{x,y}  but x' is in the integrand")
print("       (both in (x'-w')^i/(x'-w')^2 and in rho^{b'}(x')).  Should be int_{x,y,x'}.")
print("   (b) the symbol k is used for TWO different momenta in N2: the observed one")
print("       (in e^{-ik(w'-w)} and d^3k) and the loop one (in int d^2k and e^{ik(w-y)}).")
print("       Rename the loop momentum, or the w dependence will be mis-combined.")
print("   (c) the row is  [A^(1)dag] x [A^(3)] only -- one ordering.  Its value is complex")
print("       (G is attached to (x-w) but not to (x'-w')), so the physical result needs")
print("       + c.c.  i.e. 2 Re, or the companion row A^(3)dag x A^(1).")
print("   (d) minor: the definitions write A^a_{2i} and A^(3)a_{2i} while the row uses")
print("       A^(1)_i and A^(3)_i -- same object, fix the subscript.")
