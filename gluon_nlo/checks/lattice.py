"""Periodic transverse lattice, random adjoint Wilson lines, kernels.

All kernels are regularised by setting 1/0 -> 0 at coincident points, so any
algebraic identity that holds pointwise (or after relabelling dummy variables)
in the continuum holds exactly on the lattice.
"""
import numpy as np
import itertools

def su_n_generators(N):
    """Hermitian fundamental generators t^a with tr(t^a t^b)=delta/2."""
    ts = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex); m[i, j] = m[j, i] = 0.5; ts.append(m)
            m = np.zeros((N, N), complex); m[i, j] = -0.5j; m[j, i] = 0.5j; ts.append(m)
    for k in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(k):
            m[i, i] = 1
        m[k, k] = -k
        ts.append(m / np.sqrt(2 * k * (k + 1)))
    return np.array(ts)

def structure_constants(t):
    """f^{abc} = -2 i tr([t^a,t^b] t^c)."""
    n = len(t)
    f = np.zeros((n, n, n))
    for a, b, c in itertools.product(range(n), repeat=3):
        f[a, b, c] = (-2j * np.trace((t[a] @ t[b] - t[b] @ t[a]) @ t[c])).real
    return f

def random_su_n(N, rng, scale=1.0):
    t = su_n_generators(N)
    h = np.einsum('a,aij->ij', rng.normal(size=len(t)) * scale, t)
    w, v = np.linalg.eigh(h)
    return v @ np.diag(np.exp(1j * w)) @ v.conj().T

def adjoint(V, t):
    """U^{ab} = 2 tr(t^a V t^b V^dagger): real orthogonal."""
    return 2 * np.einsum('aij,jk,bkl,li->ab', t, V, t, V.conj().T).real

class Lattice:
    def __init__(self, L=5, N=3, seed=1, k=(1, 2), scale=1.5):
        assert L % 2 == 1, "odd L so that the minimal image is unique"
        self.L, self.N = L, N
        rng = np.random.default_rng(seed)
        self.rng = rng
        self.S = L * L
        self.pos = np.array([(i, j) for i in range(L) for j in range(L)])
        d = self.pos[None, :, :] - self.pos[:, None, :]          # d[a,b] = b - a
        d = (d + L // 2) % L - L // 2                             # minimal image
        self.d = d
        r2 = (d ** 2).sum(-1).astype(float)
        inv = np.where(r2 > 0, 1.0 / np.where(r2 > 0, r2, 1), 0.0)
        # KV[A,B,i] = (A-B)^i/(A-B)^2  ->  (A-B) = -(d[A,B])
        self.KV = -d * inv[..., None]
        self.INV2 = inv                                           # 1/(A-B)^2, symmetric
        # K(a,b,c) = KV(a-c).KV(b-c)
        self.KK = np.einsum('aci,bci->abc', self.KV, self.KV)
        self.t = su_n_generators(N)
        self.f = structure_constants(self.t)
        self.NA = len(self.t)
        self.T = -1j * self.f                                     # adjoint generators (T^a)_{bc}
        Vs = [random_su_n(N, rng, scale) for _ in range(self.S)]
        self.V = np.array(Vs)
        self.U = np.array([adjoint(V, self.t) for V in Vs])
        kk = 2 * np.pi * np.array(k) / L
        self.ph = np.exp(-1j * self.pos @ kk)                      # e^{-i k.w}

    def set_U_identity(self):
        self.U = np.array([np.eye(self.NA) for _ in range(self.S)])

    def set_U_constant(self, rng=None):
        rng = rng or self.rng
        V = random_su_n(self.N, rng, 1.5)
        U0 = adjoint(V, self.t)
        self.U = np.array([U0 for _ in range(self.S)])
