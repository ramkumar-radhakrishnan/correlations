"""The new row differs from the previous one ONLY in the ket colour bracket.
Check what both brackets become at large |z|, where U(z) -> 1 (adjoint identity)."""
import numpy as np
rng = np.random.default_rng(11)
N2 = 8                                  # N_c^2 - 1 for SU(3)
def rand_adj():                         # adjoint Wilson line: real orthogonal
    A = rng.normal(size=(N2,N2)); Q,_ = np.linalg.qr(A); return Q
Uy, Uw, Ux, Uyp, Uxp = [rand_adj() for _ in range(5)]
rho_x, rho_y, rho_xp = [rng.normal(size=N2) for _ in range(3)]
I = np.eye(N2)

print("FIRST (bra) bracket -- unchanged between the two rows:")
print("   [ U^{e'c'}(y') U^{d'b}(z) rho^{e'}(x')  -  U^{bd'}(z) U^{c'e'}(x') rho^{e'}(x') ]")
def bra(Uz):
    t1 = np.einsum('ec,db,e->cdb', Uyp, Uz, rho_xp)     # indices (c', d', b)
    t2 = np.einsum('bd,ce,e->cdb', Uz, Uxp, rho_xp)
    return t1 - t2
lim = bra(I); pred = np.einsum('ec,e->c', Uyp, rho_xp)[:,None,None]*np.eye(N2)[None,:,:] \
                    - np.einsum('ce,e->c', Uxp, rho_xp)[:,None,None]*np.eye(N2)[None,:,:]
print(f"   at U(z)=1 :  max|bracket - d^{{bd'}}[U^{{e'c'}}(y') - U^{{c'e'}}(x')] rho^{{e'}}(x')| = "
      f"{np.abs(lim-pred).max():.2e}")
print(f"   its magnitude: max|.| = {np.abs(lim).max():.4f}   -> z-INDEPENDENT and NONZERO")
print()

print("SECOND (ket) bracket -- PREVIOUS row:")
print("   [ U^{bc}(z) U^{ad}(y) U^{ce}(x) rho^d(y) rho^e(x)  -  U^{ad}(y) rho^d(y) rho^b(x) ]")
def ket_old(Uz):
    t1 = np.einsum('bc,ad,ce,d,e->ab', Uz, Uy, Ux, rho_y, rho_x)
    t2 = np.einsum('ad,d,b->ab', Uy, rho_y, rho_x)
    return t1 - t2
lim_old = ket_old(I)
pred_old = np.einsum('ad,d->a', Uy, rho_y)[:,None]*(np.einsum('be,e->b', Ux, rho_x) - rho_x)[None,:]
print(f"   at U(z)=1 :  max|bracket - U^{{ad}}(y)rho^d(y)[U^{{be}}(x)rho^e(x) - rho^b(x)]| = "
      f"{np.abs(lim_old-pred_old).max():.2e}")
print(f"   magnitude {np.abs(lim_old).max():.4f}  -> z-independent, nonzero (a difference in x)")
print()

print("SECOND (ket) bracket -- THIS row:")
print("   [ U^{bc}(z) U^{ce}(x) U^{ad}(y) rho^e(x) rho^d(y)")
print("     - U^{ac}(w) U^{de}(x) U^{db}(z) rho^e(x) rho^c(y) ]")
def ket_new(Uz):
    t1 = np.einsum('bc,ce,ad,e,d->ab', Uz, Ux, Uy, rho_x, rho_y)
    t2 = np.einsum('ac,de,db,e,c->ab', Uw, Ux, Uz, rho_x, rho_y)
    return t1 - t2
lim_new = ket_new(I)
# at U(z)=1 : t1 -> U^{be}(x)rho^e(x) U^{ad}(y)rho^d(y) ;  t2 -> U^{ac}(w)rho^c(y) U^{be}(x)rho^e(x)
pred_new = np.einsum('be,e->b', Ux, rho_x)[None,:]*(np.einsum('ad,d->a', Uy, rho_y)
                                                    - np.einsum('ac,c->a', Uw, rho_y))[:,None]
print(f"   at U(z)=1 :  max|bracket - U^{{be}}(x)rho^e(x)[U^{{ad}}(y) - U^{{ad}}(w)]rho^d(y)| = "
      f"{np.abs(lim_new-pred_new).max():.2e}")
print(f"   magnitude {np.abs(lim_new).max():.4f}  -> z-INDEPENDENT and NONZERO")
print()
print("   NOTE the new structure is a DIFFERENCE of Wilson lines at y and w,")
print("   so it vanishes when w -> y.  Check:")
for eps in (1e-1, 1e-2, 1e-3):
    # U(w) -> U(y) continuously
    Uw_near = Uy @ np.linalg.qr(np.eye(N2) + eps*rng.normal(size=(N2,N2)))[0]
    p = np.einsum('be,e->b', Ux, rho_x)[None,:]*(np.einsum('ad,d->a', Uy, rho_y)
                                                 - np.einsum('ac,c->a', Uw_near, rho_y))[:,None]
    print(f"      |U(w)-U(y)| ~ {eps:.0e} :  max|limit bracket| = {np.abs(p).max():.4e}")
print()
print("CONCLUSION: in BOTH rows the large-|z| colour factor is z-independent and nonzero,")
print("so the coefficient of the large-|z| logarithm survives.  The new colour structure")
print("softens the w -> y coincidence, but does nothing for the z integral.")
