"""JIMWLK block: (2.1)+(2.2) == (2.3) ?  singlet (2.4) ?  and all == -H_JIMWLK O_LO ?"""
import numpy as np, texeval as te
from lattice import Lattice

B = te.blocks('../Final_terms_original.tex')
m21 = te.parse(te.body(B[3])); m22 = te.parse(te.body(B[4])); m23 = te.parse(te.body(B[5]))
m24 = te.parse(te.singlet_body(B[6]))
print('monomials: (2.1)', len(m21), '(2.2)', len(m22), '(2.3)', len(m23), '(2.4)', len(m24))
INT = ['x', "x'", 'z', 'w', "w'"]; INT4 = ['x', 'y', 'z', 'w', "w'"]
PH = [('w', 1), ("w'", -1)]

def lo_G(lat, P):
    """O_LO = sum_{p,q} G[p,q,e',e] (U_p^T U_q)_{e'e}  from
       sum e^{-ik(w-w')} K(x-w).K(x'-w') [U(w')-U(x')]^{ce'} [U(w)-U(x)]^{ce} <rho^{e'}(x') rho^{e}(x)>
       (the author's (IV.29)/(4.2) ordering; amplitude at w, conjugate at w')."""
    S = lat.S; KV, ph = lat.KV, lat.ph
    # amplitude weights a[x,w] = ph(w) KV(x-w)  (vector), conj weights b[x',w'] = ph*(w') KV(x'-w')
    G = np.zeros((S, S, lat.NA, lat.NA), complex)
    # (p,q) = (w',w):  + sum P[x',e',x,e] KV(x'-w').KV(x-w) ph*(w') ph(w)
    core = np.einsum('vi,wi->vw', np.ones((1,1)), np.ones((1,1))) if False else None
    A = np.einsum('xwi,w->xwi', KV, ph)          # amplitude kernel (x -> w)
    C = np.einsum('xwi,w->xwi', KV, ph.conj())   # conjugate kernel (x' -> w')
    # term [U(w')][U(w)]
    G += np.einsum('uvi,xwi,uExe->vwEe', C, A, P)            # p=w', q=w
    G -= np.einsum('uvi,xwi,uExe->uwEe', C, A, P)            # p=x', q=w
    G -= np.einsum('uvi,xwi,uExe->vxEe', C, A, P)            # p=w', q=x
    G += np.einsum('uvi,xwi,uExe->uxEe', C, A, P)            # p=x', q=x
    return G

def lo_value(lat, P):
    G = lo_G(lat, P)
    return np.einsum('pqEe,pcE,qce->', G, lat.U, lat.U)

def c_matrix(lat):
    """c[z,e,f] defined by U_z^dagger T^e U_z = c^{ef} T^f   (U real orthogonal)."""
    U, T = lat.U, lat.T
    lhs = np.einsum('zba,ebc,zcd->zead', U, T, U)           # (U^T T^e U)_{ad}
    # project on T^f: tr(T^f T^g) = N delta
    return np.einsum('zead,fda->zef', lhs, T).real / lat.N

def minus_H_on_LO(lat, P):
    """-Htilde O_LO,  H = (alpha_s/2pi^2) Htilde,
    Htilde = sum_{xyz} K_xyz [J_L^e(x)J_L^e(y)+J_R^e(x)J_R^e(y) - 2 c^{ef}(z) J_L^e(x) J_R^f(y)],
    J_L^e U = T^e U, J_R^e U = U T^e (Hermitian generators -> H >= 0, d/dY = -H).
    Exact because O_LO is bilinear in U."""
    G = lo_G(lat, P); U, T, KK, N = lat.U, lat.T, lat.KK, lat.N
    c = c_matrix(lat)
    Kd = np.einsum('ppz->pz', KK)                           # K_ppz
    X = np.einsum('pcE,qce->pqEe', U, U)                     # (U_p^T U_q)
    M = Kd[:, None, :] + Kd[None, :, :] - 2 * KK              # M_pqz
    S = lat.S
    out = 0
    # C_A (M + K_pp + K_qq) X
    w1 = (M + Kd[:, None, :] + Kd[None, :, :]).sum(-1)
    out += N * np.einsum('pq,pqEe,pqEe->', w1, G, X)
    # -2 K_pq T^e X T^e
    TXT = np.einsum('gEa,pqab,gbe->pqEe', T, X, T)
    out += -2 * np.einsum('pq,pqEe,pqEe->', KK.sum(-1), G, TXT)
    # -2 c^{ef}(z) [(K_pp - K_pq) T^f U_p^T T^e U_q + (K_qq - K_pq) U_p^T T^e U_q T^f]
    A1 = np.einsum('fEa,pca,ecd,qde->pqfeEe'.replace('fe','fg'), T, U, T, U) if False else None
    # T^f U_p^T T^e U_q : indices (E,e)
    Y1 = np.einsum('fEa,pba,gbd,qde->pqgfEe', T, U, T, U)     # g=e (generator idx), f
    Y2 = np.einsum('pbE,gbd,qda,fae->pqgfEe', U, T, U, T)
    w_pp = (Kd[:, None, :] - KK)                              # [p,q,z]
    w_qq = (Kd[None, :, :] - KK)
    out += -2 * np.einsum('pqz,zgf,pqEe,pqgfEe->', w_pp, c, G, Y1, optimize=True)
    out += -2 * np.einsum('pqz,zgf,pqEe,pqgfEe->', w_qq, c, G, Y2, optimize=True)
    return -out

for N in (2, 3):
    lat = Lattice(L=5, N=N, seed=11)
    S, NA = lat.S, lat.NA
    c = c_matrix(lat)
    # gauge-consistency of the cross term: Theta(a,b,a) = N D(a,b)
    a, b = 3, 17
    Th = np.einsum('ef,ca,ecd,db,fba->', c[a], lat.U[a], lat.T, lat.U[b], lat.T)
    D = np.einsum('ce,ce->', lat.U[a], lat.U[b])
    print(f'N={N}: gauge check Theta(a,b,a)/(N D(a,b)) = {Th/(N*D):.6f}')
    rng = np.random.default_rng(5)
    P = rng.normal(size=(S, NA, S, NA)) + 1j * rng.normal(size=(S, NA, S, NA))
    ev = te.Evaluator(lat, P2=P)
    v21, v22, v23 = (ev.total(m, INT, PH) for m in (m21, m22, m23))
    ex = minus_H_on_LO(lat, P)
    print(f'  generic <rho rho>: (2.1)+(2.2)={v21+v22:.6e}  (2.3)={v23:.6e}  -H O_LO={ex:.6e}')
    print(f'     rel diffs: (2.1+2.2 vs 2.3) {abs(v21+v22-v23)/abs(v23):.1e}   (2.3 vs -H O) {abs(v23-ex)/abs(ex):.1e}')
    n = rng.normal(size=(S, S)); n = n + n.T; n -= n.mean(1)[:, None] + n.mean(0)[None, :] - n.mean()
    Pg = np.einsum('xy,bc->xbyc', n, np.eye(NA)) / NA
    evg = te.Evaluator(lat, P2=Pg)
    v23g = evg.total(m23, INT, PH)
    v24 = evg.total(m24, INT4, PH) / NA
    exg = minus_H_on_LO(lat, Pg)
    print(f'  singlet: (2.3)={v23g:.6e}  (2.4)={v24:.6e}  -H O_LO={exg:.6e}')
    print(f'     rel diffs: (2.3 vs 2.4) {abs(v23g-v24)/abs(exg):.1e}  (2.4 vs -H O) {abs(v24-exg)/abs(exg):.1e}')
