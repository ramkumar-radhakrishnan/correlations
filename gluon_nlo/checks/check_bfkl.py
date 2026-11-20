"""BFKL block: (1.1) == (1.2) ?  (1.1) with a colour-singlet projectile == (1.3) ?
(1.3) == exact one-gluon (dipole) evolution of the LO cross section ?"""
import numpy as np, texeval as te
from lattice import Lattice

B = te.blocks('../Final_terms_original.tex')
m11 = te.parse(te.body(B[0])); m12 = te.parse(te.body(B[1]))
print('monomials (1.1):', len(m11), ' (1.2):', len(m12))
INT = ['x', "x'", 'z', 'w', "w'"]
PH = [('w', 1), ("w'", -1)]

def tr_F(lat, n):
    """(1.3) evaluated directly from its trace form.  F(a,b) = sum_{w,w'} e^{-ik(w-w')}
    K(a-w).K(b-w') Tr[(U_a-U_w)^T (U_b-U_w')]   (general a,b), returned as matrix."""
    U, KV, ph = lat.U, lat.KV, lat.ph
    S = lat.S
    # A_i(a)^{ce} = sum_w ph(w) KV_i(a,w) (U_a - U_w)^{ce}
    A = np.einsum('w,awi,ace->aice', ph, KV, U) - np.einsum('w,awi,wce->aice', ph, KV, U)
    return np.einsum('aice,bice->ab', A, A.conj())

def eq13(lat, n):
    """Author's (1.3) bracket, typed from the trace form (8 lines)."""
    U, KV, ph, KK = lat.U, lat.KV, lat.ph, lat.KK
    # generic trace object  T(p,q;w,w') = Tr[(U_p - U_w)^T (U_q - U_w')]
    S = lat.S; N = lat.N
    D = np.einsum('pce,qce->pq', U, U)          # Tr[U_p^T U_q]
    # build explicitly: Tr[(U_p-U_w)^T(U_q-U_w')] = D[p,q]-D[p,w']-D[w,q]+D[w,w']
    def T4():
        return (D[:, None, :, None] - D[:, None, None, :] - D[None, :, :, None] + D[None, :, None, :])
    T = T4()                                       # T[p,w,q,w']
    php, phm = ph, ph.conj()
    tot = 0
    # line1: K(x-w).K(y-w') K(x,y,z) 2Nc Tr[(Ux-Uw)(Uy-Uw')]
    pass
    G = np.einsum('xwi,yvi,w,v,xwyv->xy', KV, KV, php, phm, T, optimize=True)   # F(x,y)
    Kxyz = KK; Kxxz = np.einsum('xxz->xz', KK); 
    tot += np.einsum('xy,xyz,xy->', n, Kxyz, G) * 2 * N
    tot += -N * np.einsum('xy,xz,xy->', n, Kxxz, G)
    tot += -N * np.einsum('xy,yz,xy->', n, Kxxz, G)
    # line4: K(w,w',z) K(x,y,z) 2Nc Tr[(Uz-Uw)(Uz-Uw')]
    # T[z,w,z,w'] needs diagonal p=q=z:
    Tzz = np.einsum('zwzv->zwv', T)
    Fzz = np.einsum('wvz,w,v,zwv->z', KK, php, phm, Tzz, optimize=True)
    tot += 2 * N * np.einsum('xy,xyz,z->', n, Kxyz, Fzz)
    # line5: K(x-w').K(w-z) K(x,y,z) 2Nc Tr[(Uz-Uw)(Ux-Uw')]
    H5 = np.einsum('xvi,wzi,w,v,zwxv->zx', KV, KV, php, phm, T, optimize=True)
    tot += 2 * N * np.einsum('xy,xyz,zx->', n, Kxyz, H5)
    tot += -N * np.einsum('xy,yz,zx->', n, Kxxz, H5)
    # line7: K(x-w).K(w'-z) K(x,y,z) 2Nc Tr[(Ux-Uw)(Uz-Uw')]
    H7 = np.einsum('xwi,vzi,w,v,xwzv->xz', KV, KV, php, phm, T, optimize=True)
    tot += 2 * N * np.einsum('xy,xyz,xz->', n, Kxyz, H7)
    tot += -N * np.einsum('xy,yz,xz->', n, Kxxz, H7)
    return tot / (lat.NA)

def bfkl_expected(lat, n):
    """Independent prediction: exact one-step evolution of an arbitrary colour-singlet
    projectile, written through dipoles (finite N_c, see report):
       sum_{x,y} n(x,y) [...]  ==  -(N_c/2) sum_{x,y} n(x,y) sum_z M_{xyz}... ;
    implemented as the BFKL kernel acting on n, contracted with F."""
    N = lat.N; KK = lat.KK; S = lat.S
    F = tr_F(lat, n)                      # F[a,b]
    Kxxz = np.einsum('xxz->xz', KK)
    # decompose n = sum over 'dipoles': n(x,y) with zero row sums.  For each ordered pair
    # (x,y), x!=y, a dipole d_{xy} = e_xx+e_yy-e_xy-e_yx; n = -sum_{x<y} n(x,y) d_{xy}
    tot = 0
    for x in range(S):
        for y in range(x + 1, S):
            c = -n[x, y]
            if c == 0: continue
            M = Kxxz[x] + Kxxz[y] - 2 * KK[x, y]          # dipole kernel over z
            O = lambda a, b: (F[a, a] + F[b, b] - F[a, b] - F[b, a])
            Oxz = F[x, x] + np.diag(F) - F[x, :] - F[:, x]
            Ozy = np.diag(F) + F[y, y] - F[:, y] - F[y, :]
            Oxy = O(x, y)
            # one-gluon real+virtual for a singlet dipole (derived in the report):
            # 2 * C_F*N/2 * M * [O_xz + O_zy - O_xy] / (N^2-1) per unit C_F normalisation
            CF = (N * N - 1) / (2 * N)
            tot += c * N * CF * np.sum(M * (Oxz + Ozy - Oxy)) / CF
    return tot / lat.NA

for N in (2, 3):
    lat = Lattice(L=5, N=N, seed=7)
    S, NA = lat.S, lat.NA
    rng = np.random.default_rng(3)
    P = rng.normal(size=(S, NA, S, NA)) + 1j * rng.normal(size=(S, NA, S, NA))
    ev = te.Evaluator(lat, P2=P)
    a, b = ev.total(m11, INT, PH), ev.total(m12, INT, PH)
    print(f'N={N} generic ordered <rho rho>:  (1.1)={a:.6e}  (1.2)={b:.6e}  rel.diff={abs(a-b)/abs(a):.1e}')
    # colour-singlet projectile, symmetric n with zero row sums
    n = rng.normal(size=(S, S)); n = n + n.T; n -= n.mean(1)[:, None] + n.mean(0)[None, :] - n.mean()
    Pg = np.einsum('xy,bc->xbyc', n, np.eye(NA)) / NA
    ev = te.Evaluator(lat, P2=Pg)
    a = ev.total(m11, INT, PH)
    c = eq13(lat, n)
    e = bfkl_expected(lat, n)
    print(f'      singlet: (1.1)={a:.6e}  (1.3)={c:.6e}  BFKL-expected={e:.6e}  diffs {abs(a-c)/abs(a):.1e} {abs(c-e)/abs(c):.1e}')
    # n with nonzero row sums (non-singlet total charge) -- does (1.3)==(1.1) still hold?
    n2 = rng.normal(size=(S, S)); n2 = n2 + n2.T
    Pg = np.einsum('xy,bc->xbyc', n2, np.eye(NA)) / NA
    a = te.Evaluator(lat, P2=Pg).total(m11, INT, PH); c = eq13(lat, n2)
    print(f'      generic symmetric n: (1.1)={a:.6e} (1.3)={c:.6e} diff {abs(a-c)/abs(a):.1e}')
