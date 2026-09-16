"""Row 8 colour algebra: does the z-dependence really cancel, and what is the residual factor?"""
import numpy as np, itertools
from scipy.linalg import expm

Nc = 3
# SU(3) structure constants
def build_f():
    lam = np.zeros((8,3,3), dtype=complex)
    lam[0][0,1]=lam[0][1,0]=1
    lam[1][0,1]=-1j; lam[1][1,0]=1j
    lam[2][0,0]=1; lam[2][1,1]=-1
    lam[3][0,2]=lam[3][2,0]=1
    lam[4][0,2]=-1j; lam[4][2,0]=1j
    lam[5][1,2]=lam[5][2,1]=1
    lam[6][1,2]=-1j; lam[6][2,1]=1j
    lam[7]=np.diag([1,1,-2])/np.sqrt(3)
    T = lam/2
    f = np.zeros((8,8,8))
    for a,b,c in itertools.product(range(8),repeat=3):
        f[a,b,c] = (-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
    return T, f
T, f = build_f()

def adjoint(alpha):
    """random adjoint Wilson line U^{ab} = 2 tr(t^a V t^b V^dag), real orthogonal."""
    H = np.zeros((3,3), dtype=complex)
    for a in range(8): H += alpha[a]*T[a]
    V = expm(1j*H)
    U = np.zeros((8,8))
    for a in range(8):
        for b in range(8):
            U[a,b] = (2*np.trace(T[a]@V@T[b]@V.conj().T)).real
    return U

rng = np.random.default_rng(7)
Uz  = adjoint(rng.normal(size=8))
Uy  = adjoint(rng.normal(size=8))
Uyp = adjoint(rng.normal(size=8))
Ux  = adjoint(rng.normal(size=8))
Uxp = adjoint(rng.normal(size=8))

print("orthogonality  ||U^T U - 1|| =", np.abs(Uz.T@Uz-np.eye(8)).max())

# ---------- bra bracket, indices exactly as written ----------
# [ f^{c'd'a} U^{e'c'}(y') U^{d'b}(z)  -  f^{c'd'a} U^{bd'}(z) U^{c'e'}(x') ]   (free: a,b,e')
BRA = (np.einsum('xya,ex,yb->abe', f, Uyp, Uz)
     - np.einsum('xya,by,xe->abe', f, Uz, Uxp))
# ---------- ket bracket, indices exactly as written ----------
# [ f^{cda} U^{ec}(y) U^{db}(z)  -  f^{cda} U^{bd}(z) U^{ce}(x) ]
KET = (np.einsum('xya,ex,yb->abe', f, Uy, Uz)
     - np.einsum('xya,by,xe->abe', f, Uz, Ux))

full_asis = np.einsum('abe,abE->eE', BRA, KET)   # sum a,b ; free e',e

# ---------- ket bracket with the transposes made uniform ----------
KET2 = (np.einsum('xya,ex,yb->abe', f, Uy, Uz)
      - np.einsum('xya,yb,ex->abe', f, Uz, Ux))
full_fix = np.einsum('abe,abE->eE', BRA, KET2)

# ---------- bra bracket also made uniform ----------
BRA2 = (np.einsum('xya,ex,yb->abe', f, Uyp, Uz)
      - np.einsum('xya,yb,ex->abe', f, Uz, Uxp))
full_fix2 = np.einsum('abe,abE->eE', BRA2, KET2)

# ---------- candidate reduced forms ----------
red_noNc = np.einsum('Ec,ec->Ee', Uyp-Uxp, Uy-Ux)          # [U(y')-U(x')][U(y)-U(x)]
red_Nc   = Nc*red_noNc

print()
print("as written (mixed transposes) vs Nc[U(y')-U(x')][U(y)-U(x)] : maxdiff %.3e" %
      np.abs(full_asis-red_Nc).max())
print("ket transposes made uniform                                : maxdiff %.3e" %
      np.abs(full_fix-red_Nc).max())
print("both brackets uniform                                      : maxdiff %.3e" %
      np.abs(full_fix2-red_Nc).max())
print("both uniform, compared WITHOUT the Nc                      : maxdiff %.3e" %
      np.abs(full_fix2-red_noNc).max())
print()
print("ratio  full_uniform / reduced   (elementwise, first 3x3):")
print(np.round((full_fix2/red_noNc)[:3,:3], 8))
print()
print("f^{c'da} f^{cda} = Nc delta^{c'c} :", np.allclose(np.einsum('xya,zya->xz',f,f), Nc*np.eye(8)))
print("U^{d'b}(z) U^{db}(z) = delta^{d'd}:", np.allclose(Uz@Uz.T, np.eye(8)))
print("U^{bd'}(z) U^{bd}(z) = delta^{d'd}:", np.allclose(Uz.T@Uz, np.eye(8)))
print("U^{d'b}(z) U^{bd}(z) = (UU)^{d'd} NOT delta:", np.abs(Uz@Uz-np.eye(8)).max())
