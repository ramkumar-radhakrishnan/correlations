"""Is the new colour structure the old one with U(z) -> 1, or a unitarity contraction?"""
import numpy as np
rng = np.random.default_rng(5)
N2 = 8
def rand_adj():
    A = rng.normal(size=(N2,N2)); Q,_ = np.linalg.qr(A); return Q
Uy, Uw, Ux, Uyp, Uxp, Uz = [rand_adj() for _ in range(6)]
rx, ry, rxp = [rng.normal(size=N2) for _ in range(3)]
I = np.eye(N2)
# totally antisymmetric f (any antisymmetric tensor suffices for a linear identity check)
F = rng.normal(size=(N2,N2,N2))
F = (F - F.transpose(1,0,2) + F.transpose(1,2,0) - F.transpose(2,1,0)
     + F.transpose(2,0,1) - F.transpose(0,2,1))/6

def old_total(Uz):
    # Bra^{c' d' b}
    bra = (np.einsum('ec,db,e->cdb', Uyp, Uz, rxp)
           - np.einsum('bd,ce,e->cdb', Uz, Uxp, rxp))
    # Ket^{a b}
    ket = (np.einsum('bc,ce,ad,e,d->ab', Uz, Ux, Uy, rx, ry)
           - np.einsum('ac,de,db,e,c->ab', Uw, Ux, Uz, rx, ry))
    return np.einsum('CDA,CDB,AB->', F, bra, ket)          # f^{c'd'a} Bra Ket, summed over b

def new_total():
    bra = np.einsum('ec,e->c', Uyp, rxp) - np.einsum('ce,e->c', Uxp, rxp)      # ^{c'}
    ket = (np.einsum('ce,e->c', Ux, rx)[None,:]
           * (np.einsum('ad,d->a', Uy, ry) - np.einsum('ad,d->a', Uw, ry))[:,None])  # ^{a c}
    return np.einsum('CBA,C,AB->', F, bra, ket)            # f^{c' c a} Bra^{c'} Ket^{c a}

print("A. is the new colour structure the OLD one with U(z) -> 1 ?")
print(f"   old at U(z) = 1 : {old_total(I):+.12f}")
print(f"   new             : {new_total():+.12f}")
print(f"   difference      : {abs(old_total(I)-new_total()):.2e}    <-- 0 means YES")
print(f"   old at a generic U(z) : {old_total(Uz):+.12f}   (quite different)")
print()

print("B. could it instead follow from UNITARITY of U(z)?  Check the four products.")
print("   (adjoint U is real orthogonal: U^{ab}U^{cb} = U^{ba}U^{bc} = delta^{ac})")
pairs = [("bra1 x ket1 : U^{d'b}(z) U^{bc}(z)", np.einsum('db,bc->dc', Uz, Uz)),
         ("bra1 x ket2 : U^{d'b}(z) U^{db}(z)", np.einsum('eb,db->ed', Uz, Uz)),
         ("bra2 x ket1 : U^{bd'}(z) U^{bc}(z)", np.einsum('bd,bc->dc', Uz, Uz)),
         ("bra2 x ket2 : U^{bd'}(z) U^{db}(z)", np.einsum('bd,eb->de', Uz, Uz))]
for nm, M in pairs:
    print(f"   {nm:38s}  max|M - delta| = {np.abs(M-I).max():.3f}"
          f"   {'-> delta (unitarity)' if np.abs(M-I).max()<1e-10 else '-> U(z)^2, NOT delta'}")
print()
print("   so only the CROSS terms collapse by unitarity; the diagonal ones give U(z)^2.")
print("   The quoted expression is therefore NOT a unitarity rearrangement --")
print("   it is the U(z) -> 1 limit.")
