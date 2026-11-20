# Colour framework: SU(N) generators, f, d, random adjoint Wilson lines, traces.
import numpy as np
def gellmann(N):
    ts = []
    for j in range(N):
        for k in range(j+1, N):
            m = np.zeros((N, N), complex); m[j, k] = m[k, j] = 1; ts.append(m/2)
            m = np.zeros((N, N), complex); m[j, k] = -1j; m[k, j] = 1j; ts.append(m/2)
    for l in range(1, N):
        m = np.zeros((N, N), complex)
        for q in range(l): m[q, q] = 1
        m[l, l] = -l; ts.append(m/np.sqrt(2*l*(l+1)))
    return np.array(ts)
class SUN:
    def __init__(s, N):
        s.N = N; t = s.t = gellmann(N); s.n = len(t)
        AB = np.einsum('aij,bjk->abik', t, t)
        s.f = np.real(-2j*np.einsum('abij,cji->abc', AB-AB.transpose(1, 0, 2, 3), t))
        s.d = np.real(2*np.einsum('abij,cji->abc', AB+AB.transpose(1, 0, 2, 3), t))
        s.T = -1j*s.f          # (T^e)_{bc} = -i f^{ebc}
        s.D = s.d              # (D^e)_{bc} = d^{ebc}
    def randV(s, rng):
        A = rng.normal(size=(s.N, s.N))+1j*rng.normal(size=(s.N, s.N))
        Q, R = np.linalg.qr(A); Q = Q@np.diag(np.diag(R)/np.abs(np.diag(R)))
        return Q/np.linalg.det(Q)**(1.0/s.N)
    def adj(s, V):
        return np.real(2*np.einsum('aij,jk,bkl,li->ab', s.t, V, s.t, V.conj().T))
# traces in the notation of the explicit file (adjoint U real orthogonal, U^dagger = U^T)
def Dtr(S, U, a, b): return np.trace(U[a].T@U[b])
def Ftr(S, U, a, b, c): return np.einsum('eij,jk,kl,elm,mn,ni->', S.T, U[a].T, U[b], S.T, U[b].T, U[c], optimize=True).real
def Qtr(S, U, a, b, c, d): return np.einsum('eij,jk,kl,elm,mn,ni->', S.T, U[a].T, U[b], S.T, U[c].T, U[d], optimize=True).real
def Trd(S, U, a, b, c, e):   # Tr{T^h U^dag(a) U(b) D^h U^dag(c) U(e)}  (purely imaginary)
    return np.einsum('hij,jk,kl,hlm,mn,ni->', S.T, U[a].T, U[b], S.D, U[c].T, U[e], optimize=True)
def Gdef(S, U, a, b, c): return np.einsum('ijk,lmn,il,jm,kn->', S.f, S.d, U[a], U[b], U[c], optimize=True)
def Fdef(S, U, a, b, c): return np.einsum('ijk,lmn,il,jm,kn->', S.f, S.f, U[a], U[b], U[c], optimize=True)
