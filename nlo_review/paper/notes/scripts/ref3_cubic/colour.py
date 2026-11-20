# Self-contained SU(N) colour algebra and the thirteen cubic colour tensors of the finite terms,
# transcribed from Appendix app:fin_sym of the paper draft (app_finite.tex, eqs. app_sym_*).
# Points: x, y, xp (= x'), yp (= y'), z, w, wp (= w').
import itertools
import numpy as np


def gellmann(N):
    ts = []
    for j in range(N):
        for k in range(j + 1, N):
            m = np.zeros((N, N), complex); m[j, k] = m[k, j] = 0.5; ts.append(m)
            m = np.zeros((N, N), complex); m[j, k] = -0.5j; m[k, j] = 0.5j; ts.append(m)
    for l in range(1, N):
        m = np.zeros((N, N), complex)
        for q in range(l): m[q, q] = 1
        m[l, l] = -l; ts.append(m / np.sqrt(2 * l * (l + 1)))
    return np.array(ts)


class SUN:
    """Generators t^a (tr t^a t^b = delta/2), f^{abc}, d^{abc}, adjoint (T^a)_{bc} = -i f^{abc}."""
    def __init__(s, N):
        s.N = N; s.t = gellmann(N); s.n = len(s.t)
        t = s.t
        AB = np.einsum('aij,bjk->abik', t, t)
        s.f = np.real(-2j * np.einsum('abij,cji->abc', AB - AB.transpose(1, 0, 2, 3), t))
        s.d = np.real(2 * np.einsum('abij,cji->abc', AB + AB.transpose(1, 0, 2, 3), t))
        s.T = -1j * s.f
        # eta_a: (t^a)^T = eta_a t^a (Gell-Mann basis: symmetric +1, antisymmetric -1)
        s.eta = np.array([1.0 if np.allclose(m.T, m) else -1.0 for m in t])
        assert all(np.allclose(t[a].T, s.eta[a] * t[a]) for a in range(s.n))
        s.CF = (N * N - 1) / (2 * N)

    def randV(s, rng, eps=None):
        if eps is None:
            A = rng.normal(size=(s.N, s.N)) + 1j * rng.normal(size=(s.N, s.N))
            Q, R = np.linalg.qr(A); Q = Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))
            return Q / np.linalg.det(Q) ** (1.0 / s.N)
        al = rng.normal(size=s.n)                      # weak field V = exp(i eps alpha.t)
        H = np.einsum('a,aij->ij', al, s.t)
        w, v = np.linalg.eigh(H)
        return v @ np.diag(np.exp(1j * eps * w)) @ v.conj().T

    def adj(s, V):
        return np.real(2 * np.einsum('aij,jk,bkl,li->ab', s.t, V, s.t, V.conj().T))


POINTS = ['x', 'y', 'xp', 'yp', 'z', 'w', 'wp']


def structures(G, U):
    """The 13 cubic tensors C^{abc}: name -> (tensor, type, rho points (p1,p2,p3)).
    type: 'ord' = rho(p1)rho(p2)rho(p3); 'A23' = rho(p1){rho(p2),rho(p3)}; 'A12' = {rho(p1),rho(p2)}rho(p3)."""
    f = G.f; e = lambda *a: np.einsum(*a, optimize=True)
    Lx = U['xp'] - U['wp']                                    # [U^{ga}(x') - U^{ga}(w')]
    out = {}
    # Row I.II, N_2b  (eq. app_sym_IIb)
    C = e('ga,hbc,gh->abc', Lx, f, U['w']) - e('ga,gde,db,ec->abc', Lx, f, U['x'], U['y'])
    out['I.II.N2b'] = (C, 'A23', ('xp', 'x', 'y'))
    # Row I.III, C_D (u=z), C_A (u=w)  (eq. app_sym_IIII)
    for u, lab in (('z', 'I.III(z)'), ('w', 'I.III(w)')):
        C = e('ga,ghd,db,hc->abc', Lx, f, U['x'], U[u] - U['y'])
        out[lab] = (C, 'ord', ('xp', 'x', 'y'))
    # Row I.V, first part (eq. app_sym_IV1)
    C = e('ga,gde,dc,eb->abc', Lx, f, U['x'], U['z']) - e('ga,gde,eb,dc->abc', Lx, f, U['y'], U['yp'])
    out['I.V.N1'] = (C, 'A23', ('xp', 'y', 'yp'))
    # Row II.I, N_1a (eq. app_sym_III1a)
    C = (e('ga,kbc,gk->abc', U['wp'], f, U['w']) - e('ga,gkl,kb,lc->abc', U['wp'], f, U['z'], U['x'])
         - e('ga,mh,mb,khc,gk->abc', U['xp'], U['z'], U['y'], f, U['w'])
         + e('ga,mh,mb,gkl,kh,lc->abc', U['xp'], U['z'], U['y'], f, U['z'], U['x']))
    out['II.I.N1a'] = (C, 'A12', ('xp', 'y', 'x'))
    # Row II.I, N_1b (eq. app_sym_III1b)
    C = (e('gb,kca,gk->abc', U['w'], f, U['wp']) - e('gb,gkl,kc,la->abc', U['w'], f, U['z'], U['x'])
         - e('gb,mh,mc,kha,gk->abc', U['xp'], U['z'], U['y'], f, U['wp'])
         + e('gb,mh,mc,gkl,kh,la->abc', U['xp'], U['z'], U['y'], f, U['z'], U['x']))
    out['II.I.N1b'] = (C, 'A23', ('x', 'xp', 'y'))
    # Row III.II (eq. app_sym_IIIII)
    C = (e('ma,gb,mh,khc,gk->abc', U['y'], U['xp'], U['z'], f, U['w'])
         - e('ma,gb,mh,gkl,kh,lc->abc', U['y'], U['xp'], U['z'], f, U['z'], U['x'])
         - e('gb,kac,gk->abc', U['xp'], f, U['w']) + e('gb,gkl,ka,lc->abc', U['xp'], f, U['z'], U['x']))
    out['III.II'] = (C, 'ord', ('y', 'xp', 'x'))
    # Rows III.III (A23), III.V (ord), III.VI (ord), v = y' or z
    for v, vs in (('yp', "y'"), ('z', 'z')):
        Lv = U[v] - U['xp']
        C = e('mnk,ma,nc,kb->abc', f, Lv, U['z'], U['w']) - e('mnk,ma,kb,nc->abc', f, Lv, U['y'], U['x'])
        out['III.III.1(%s)' % vs] = (C, 'A23', ('xp', 'y', 'x'))
        C = e('mkl,ma,kb,lc->abc', f, Lv, U['x'], U['y'] - U['w'])
        out['III.V(%s)' % vs] = (C, 'ord', ('xp', 'x', 'y'))
        C = e('mnk,ma,kb,nc->abc', f, Lv, U['y'], U['x'] - U['z'])
        out['III.VI(%s)' % vs] = (C, 'ord', ('xp', 'y', 'x'))
    return out


ORDER = ['I.II.N2b', 'I.III(z)', 'I.III(w)', 'I.V.N1', 'II.I.N1a', 'II.I.N1b', 'III.II',
         "III.III.1(y')", 'III.III.1(z)', "III.V(y')", 'III.V(z)', "III.VI(y')", 'III.VI(z)']


def tdelta(typ):
    return 1 if typ == 'ord' else 2


# ---- invariants ----
def Dinv(U, a, b): return np.einsum('ij,ij->', U[a], U[b])
def Finv(G, U, a, b, c): return np.einsum('ijk,lmn,il,jm,kn->', G.f, G.f, U[a], U[b], U[c])
def Ginv(G, U, a, b, c): return np.einsum('ijk,lmn,il,jm,kn->', G.f, G.d, U[a], U[b], U[c])
def Qdinv(G, U, a, b, c, e):
    A = U[a].T @ U[b]; B = U[c].T @ U[e]
    return np.einsum('hij,hkl,jk,li->', G.f, G.d, A, B)
def Qinv(G, U, a, b, c, e):
    A = U[a].T @ U[b]; B = U[c].T @ U[e]
    return np.einsum('eij,jk,ekl,li->', G.T, A, G.T, B).real


def delta_table(G, U):
    """Delta = sum over orderings C d, as printed in fix_dabc (table) / paper eq. Deltas (times the factors)."""
    g = lambda a, b, c: Ginv(G, U, a, b, c)
    q = lambda a, b, c, e: Qdinv(G, U, a, b, c, e)
    return {
        'I.II.N2b': 2 * (g('x', 'y', 'wp') - g('x', 'y', 'xp')),
        'I.III(z)': g('x', 'y', 'xp') - g('x', 'y', 'wp') + g('x', 'xp', 'z') + g('x', 'z', 'wp'),
        'I.III(w)': g('x', 'y', 'xp') - g('x', 'y', 'wp') + g('x', 'xp', 'w') + g('x', 'w', 'wp'),
        'I.V.N1': -2 * (g('x', 'xp', 'z') + g('x', 'z', 'wp') + g('y', 'xp', 'yp') + g('y', 'yp', 'wp')),
        'II.I.N1a': 2 * (g('x', 'z', 'wp') - g('x', 'y', 'xp') - q('z', 'y', 'xp', 'w')),
        'II.I.N1b': 2 * (g('x', 'z', 'w') - g('x', 'y', 'xp') - q('z', 'y', 'xp', 'wp')),
        'III.II': g('x', 'y', 'xp') + g('x', 'xp', 'z') + q('z', 'y', 'xp', 'w'),
        "III.III.1(y')": 2 * (g('x', 'y', 'xp') - g('x', 'y', 'yp') - g('xp', 'z', 'w') + g('yp', 'z', 'w')),
        'III.III.1(z)': 2 * (g('x', 'y', 'xp') - g('x', 'y', 'z') - g('xp', 'z', 'w')),
        "III.V(y')": -g('x', 'y', 'xp') + g('x', 'y', 'yp') - g('x', 'xp', 'w') + g('x', 'yp', 'w'),
        'III.V(z)': -g('x', 'y', 'xp') + g('x', 'y', 'z') - g('x', 'xp', 'w') + g('x', 'z', 'w'),
        "III.VI(y')": -g('x', 'y', 'xp') + g('x', 'y', 'yp') + g('y', 'xp', 'z') - g('y', 'yp', 'z'),
        'III.VI(z)': -g('x', 'y', 'xp') + g('x', 'y', 'z') + g('y', 'xp', 'z'),
    }
