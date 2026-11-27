r"""Numerical check of the 4D (k,p) Fourier transforms used for the one-gluon LCWF.

FT[F](a,b) = (2pi)^-2 \int d^2k d^2p e^{i k.a + i p.b} F(k,p), with a = X' = x-z', b = -X = z-x.
F = (u.K)^m (v.K)^l / (Q1 Q2), K=(k,p). Schwinger parameters 1/(Q1 Q2) = \int ds dt e^{-sQ1-tQ2};
the Gaussian (k,p) integral is done exactly and (ln s, ln t) by the trapezoidal rule.
Euclidean regime: p^2 - xi k^2 -> p^2 + eta k^2 (xi = -eta).
"""
import numpy as np

def gauss_ft(Q1, Q2, u, v, a, b, L=16, n=900):
    ls = np.linspace(-L, L, n); h = ls[1]-ls[0]
    S, T = np.meshgrid(np.exp(ls), np.exp(ls), indexing='ij')
    A = S[..., None, None]*Q1 + T[..., None, None]*Q2
    det = A[..., 0, 0]*A[..., 1, 1]-A[..., 0, 1]*A[..., 1, 0]
    B = np.empty_like(A)
    B[..., 0, 0] = A[..., 1, 1]; B[..., 1, 1] = A[..., 0, 0]
    B[..., 0, 1] = -A[..., 0, 1]; B[..., 1, 0] = -A[..., 1, 0]
    B /= det[..., None, None]
    J = np.stack([a, b])
    BJ = np.einsum('...ij,jc->...ic', B, J)
    E = np.einsum('ic,...ic->...', J, BJ)/4
    uBv = np.einsum('i,...ij,j->...', u, B, v)
    uBJ = np.einsum('i,...ic->...c', u, BJ); vBJ = np.einsum('i,...ic->...c', v, BJ)
    out = uBv[..., None, None]*np.eye(2)/2 - uBJ[..., :, None]*vBJ[..., None, :]/4
    w = np.pi**2/det*np.exp(-E)*S*T
    return np.einsum('st,stml->ml', w, out)*h*h/(2*np.pi)**2

rng = np.random.default_rng(1)
a = rng.normal(size=2); b = rng.normal(size=2)
T0 = gauss_ft(np.diag([1., 0.]), np.diag([0., 1.]), np.array([1., 0.]), np.array([0., 1.]), a, b)
print("sanity k^m p^l/(k^2p^2):", np.round(T0/(-np.outer(a, b)/(a@a*b@b)), 5).ravel())

eta = 0.7
for trial in range(3):
    Xp = rng.normal(size=2); X = rng.normal(size=2)
    a, b = Xp, -X
    Q2 = np.array([[eta, 0.], [0., 1.]])                       # p^2 + eta k^2
    # case 1: (p-k)^m (p+eta k)^l / ((p-k)^2 (p^2+eta k^2))
    T1 = gauss_ft(np.array([[1., -1.], [-1., 1.]]), Q2, np.array([-1., 1.]), np.array([eta, 1.]), a, b)
    S1 = np.outer(Xp+eta*X, Xp-X)/((Xp@Xp+eta*X@X)*((Xp-X)@(Xp-X)))
    # case 2: (k-p)^m (p+eta k)^l / ((p+eta k)^2 (p^2+eta k^2))
    T2 = gauss_ft(np.array([[eta**2, eta], [eta, 1.]]), Q2, np.array([1., -1.]), np.array([eta, 1.]), a, b)
    S2 = np.outer(Xp+eta*X, Xp-X)/((Xp@Xp+eta*X@X)*((Xp+eta*X)@(Xp+eta*X)))
    print("trial", trial, " case1 FT/S1:", np.round(T1/S1, 4).ravel(), " case2 FT/S2:", np.round(T2/S2, 4).ravel())
