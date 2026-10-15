"""First-principles evaluation of the Group II / III integrands of Calculation_main.tex.

Everything is evaluated pointwise: transverse points are random 2-vectors, the charges rho^a(x) are
random classical 8-vectors (charges at different points commute), the adjoint Wilson lines are random
SU(3) adjoint matrices (general, U^T != U, or symmetric U = U^T when SYM is set).

Coefficients (Calculation_main.tex):
  A^a_i(P,-W) = i g/(sqrt2 pi sqrt P) int_X K^i(X-W) rho^a(X)                       (A1coeff, unitarity)
  Abar: rho -> U rho
  B2^{bc}_{ji}(k,-W;p,-Z): gluon (c,i,k) at W, gluon (b,j,p) at Z
      one-rho: i g^2 f^{cba} rho^a(X) sqrt(pk) / (4 pi^2 D (p+k)) * Phi^{ij},  D = p (X-Z)^2 + k (X-W)^2
      two-rho: - g^2 {rho^c(Y), rho^b(X)} K^i(Y-W) K^j(X-Z) / (8 pi^2 sqrt(pk))
  C^{dac}_{1jik}(P, y-W; p, W-Z) = - g f^{cda}/(2 pi sqrt(2P)) sqrt(pk)/P delta^2(y-W+(p/P)(W-Z))
                                   (W-Z)^m/(W-Z)^2 [d_km d_ji - (P/p) d_jm d_ki - (P/k) d_im d_jk]
Amplitude pieces (2/sqrt(2pi) times):
  frakB = U^{ac}(W) B2^{bc}_{ji} - U^{dag bc}(Z) B2bar^{ca}_{ji}
  X1'   = U^{dag bc}(Z) Abar^c_j(p,-Z) Abar^a_i(k,-W) - U^{ac}(W) Abar^d_j(p,-Z) U^{db}(Z) A^c_i(k,-W)
          (main.tex: U^{dag db}(Z); the unobserved gluon is created by abar^dag = U a^dag, which gives U^{db})
  X2'   = U^{dag bc}(Z) Abar^a_i(k,-W) Abar^c_j(p,-Z) - Abar^a_i(k,-W) A^b_j(p,-Z)
  C'    = int_y [ C^{dac} U^{dag bd}(Z) U^{ce}(y) A^e_k(P,-y) - U^{dag bd}(Z) C^{dac} Abar^c_k(P,-y) ]   (frakC = -C')
Cross section (Eq. CS2): 4/(2pi)^3 int dp/2pi int e^{-ik(w'-w)} { B^dag B + 1/4 X1'^dag X1' + 1/4 X2'^dag X2' + C'^dag C' }
  Group III: 4/(2pi)^3 int dp/2pi int e^{-ik(w'-w)} { 1/2 X1'^dag B + 1/2 X2'^dag B - C'^dag B + 1/4 X2'^dag X1'
             - 1/2 C'^dag X1' - 1/2 C'^dag X2' } + c.c.
All arrays are indexed [a, b, i, j]: a,i = measured gluon (k+, at W), b,j = unobserved gluon (p+, at Z)."""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'three_rho'))
import numerics as nm  # noqa: E402

TS, F = nm.structure_constants(3)
NC = 3
G = 1.3                     # coupling (kept != 1 to test the powers of g)
PI = np.pi
I2 = np.eye(2)


class Cfg:
    def __init__(self, seed, sym=False, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')):
        rng = np.random.default_rng(seed)
        self.rng = rng
        self.pt = {n: rng.normal(size=2) * 0.8 for n in names}
        self.U = {n: nm.random_adjoint(3, TS, rng, sym) for n in names}
        self.rho = {n: rng.normal(size=8) for n in names}
        self.kv = rng.normal(size=2)          # transverse momentum of the measured gluon
        self.sym = sym

    def add(self, name, vec):
        """a new point (e.g. w' fixed by a delta function) with its own Wilson line"""
        self.pt[name] = np.asarray(vec, float)
        self.U[name] = nm.random_adjoint(3, TS, self.rng, self.sym)
        self.rho[name] = self.rng.normal(size=8)


def sq(v):
    return float(v @ v)


def K(v):
    return v / sq(v)


def alpha(P):
    return G / (np.sqrt(2) * PI * np.sqrt(P))


# ------------------------------------------------------------------ coefficients
def A(c, P, W, X):
    """A^a_i(P,-W) with charge at X: [a, i]"""
    return 1j * alpha(P) * np.outer(c.rho[X], K(c.pt[X] - c.pt[W]))


def Abar(c, P, W, X):
    return 1j * alpha(P) * np.outer(c.U[X] @ c.rho[X], K(c.pt[X] - c.pt[W]))


def Phi(c, k, W, p, Z, X):
    """Phi^{ij} of B2 (i: gluon k at W, j: gluon p at Z)"""
    x, z, w = c.pt[X], c.pt[Z], c.pt[W]
    xz, xw, zw = x - z, x - w, z - w
    P = p + k
    return (I2 * (sq(xz) - sq(xw)) / (2 * sq(zw))
            + P / k * (np.outer(zw, xz) / sq(zw) - np.outer(xw, xz) / (2 * sq(xw)))
            + P / p * (np.outer(xw, zw) / sq(zw) + np.outer(xw, xz) / (2 * sq(xz))))


def B2one(c, k, W, p, Z, X, bar=False):
    """B2^{bc}_{ji}(k,-W;p,-Z), one-rho part, charge at X: [b, c, j, i]"""
    x, z, w = c.pt[X], c.pt[Z], c.pt[W]
    D = p * sq(x - z) + k * sq(x - w)
    r = c.U[X] @ c.rho[X] if bar else c.rho[X]
    col = np.einsum('cba,a->bc', F, r)
    coef = 1j * G**2 * np.sqrt(p * k) / (4 * PI**2 * D * (p + k))
    return coef * np.einsum('bc,ij->bcji', col, Phi(c, k, W, p, Z, X))


def B2two(c, k, W, p, Z, Y, X, bar=False):
    """two-rho part: charge Y with gluon k at W, charge X with gluon p at Z ({rho,rho} = 2 rho rho)"""
    rY = c.U[Y] @ c.rho[Y] if bar else c.rho[Y]
    rX = c.U[X] @ c.rho[X] if bar else c.rho[X]
    coef = -G**2 / (8 * PI**2 * np.sqrt(p * k)) * 2
    return coef * np.einsum('c,b,i,j->bcji', rY, rX, K(c.pt[Y] - c.pt[W]), K(c.pt[X] - c.pt[Z]))


def Cstrip(c, P, p, W, Z):
    """C^{dac}_{1jik}(P, y-W; p, W-Z) without its delta function: [d, a, c, j, i, k]"""
    k = P - p
    wz = c.pt[W] - c.pt[Z]
    v = wz / sq(wz)
    br = (np.einsum('km,ji->jikm', I2, I2) - P / p * np.einsum('jm,ki->jikm', I2, I2)
          - P / k * np.einsum('im,jk->jikm', I2, I2))
    vec = np.einsum('jikm,m->jik', br, v)
    pre = -G / (2 * PI * np.sqrt(2 * P)) * np.sqrt(p * k) / P
    fcda = np.einsum('cda->dac', F)
    return pre * np.einsum('dac,jik->dacjik', fcda, vec)


# ------------------------------------------------------------------ amplitude pieces [a, b, i, j]
def frakB_one(c, k, p, W, Z, X):
    B = B2one(c, k, W, p, Z, X)
    Bb = B2one(c, k, W, p, Z, X, bar=True)
    t1 = np.einsum('ac,bcji->abij', c.U[W], B)
    t2 = np.einsum('cb,caji->abij', c.U[Z], Bb)     # U^{dag bc}(Z) = U^{cb}(Z)
    return t1 - t2


def frakB_two(c, k, p, W, Z, Y, X):
    B = B2two(c, k, W, p, Z, Y, X)
    Bb = B2two(c, k, W, p, Z, Y, X, bar=True)
    return np.einsum('ac,bcji->abij', c.U[W], B) - np.einsum('cb,caji->abij', c.U[Z], Bb)


def X1p(c, k, p, W, Z, Yk, Xp):
    """Yk: charge of the measured gluon, Xp: charge of the unobserved gluon"""
    t1 = np.einsum('cb,cj,ai->abij', c.U[Z], Abar(c, p, Z, Xp), Abar(c, k, W, Yk))
    # unobserved gluon created after the target by abar^dag = U a^dag: factor U^{db}(Z) (main.tex writes U^{dag db})
    t2 = np.einsum('ac,dj,db,ci->abij', c.U[W], Abar(c, p, Z, Xp), c.U[Z], A(c, k, W, Yk))
    return t1 - t2


def X2p(c, k, p, W, Z, Yk, Xp):
    t1 = np.einsum('cb,ai,cj->abij', c.U[Z], Abar(c, k, W, Yk), Abar(c, p, Z, Xp))
    t2 = np.einsum('ai,bj->abij', Abar(c, k, W, Yk), A(c, p, Z, Xp))
    return t1 - t2


def Cp_strip(c, k, p, W, Z, Ypar, X):
    """C' without the delta function delta^2(Ypar - W + (p/P)(W - Z)); parent at Ypar, its charge at X"""
    P = p + k
    Cs = Cstrip(c, P, p, W, Z)
    Kk = K(c.pt[X] - c.pt[Ypar])
    r1 = c.U[Ypar] @ c.rho[X]
    r2 = c.U[X] @ c.rho[X]
    t = np.einsum('dacjik,db,c,k->abij', Cs, c.U[Z], r1 - r2, Kk)   # U^{dag bd}(Z) = U^{db}(Z)
    return 1j * alpha(P) * t


def dot(Y, Z):
    """sum_{a,b,i,j} [Y]^dagger Z  (charges classical, U real)"""
    return np.sum(np.conj(Y) * Z)


PREF = 4 / (2 * PI)**3          # times int dp/2pi (kept outside)
