# Common library for the referee-3 finiteness checks (contact kernels, conditional convergence).
# SU(3) adjoint Wilson lines (smooth fields), singlet invariants D, F, Q, and all kernels of the
# (2rho) contact equations of sec_finite.tex, transcribed from sec_finite.tex / app_finite.tex.
import numpy as np
from scipy.special import sici, exp1
from scipy.linalg import expm

GE = np.euler_gamma
KVEC = np.array([1.0, 0.0])          # measured momentum k (|k| = 1 along the first axis)
EPERP = np.array([0.0, 1.0])         # unit vector perpendicular to k

# ---------------------------------------------------------------- colour
def gellmann(N):
    ts = []
    for j in range(N):
        for k in range(j + 1, N):
            m = np.zeros((N, N), complex); m[j, k] = m[k, j] = 1; ts.append(m / 2)
            m = np.zeros((N, N), complex); m[j, k] = -1j; m[k, j] = 1j; ts.append(m / 2)
    for l in range(1, N):
        m = np.zeros((N, N), complex)
        for q in range(l): m[q, q] = 1
        m[l, l] = -l; ts.append(m / np.sqrt(2 * l * (l + 1)))
    return np.array(ts)

NC = 3
TF = gellmann(NC)
_AB = np.einsum('aij,bjk->abik', TF, TF)
FABC = np.real(-2j * np.einsum('abij,cji->abc', _AB - _AB.transpose(1, 0, 2, 3), TF))
DABC = np.real(2 * np.einsum('abij,cji->abc', _AB + _AB.transpose(1, 0, 2, 3), TF))
TADJ = -1j * FABC                      # (T^e)_{bc} = -i f^{ebc}
NA = NC * NC - 1


class SmoothField:
    """Smooth adjoint Wilson-line field U(x) = exp(A^a(x) f^a), A^a = sum of plane waves, times an
    optional envelope (finite target: U -> 1 at large distance)."""
    def __init__(self, seed=1, nmodes=4, amp=1.2, qscale=1.0, envelope=None):
        rng = np.random.default_rng(seed)
        self.q = rng.normal(size=(nmodes, 2)) * qscale
        self.c = rng.normal(size=(nmodes, NA)) * amp / np.sqrt(nmodes)
        self.ph = rng.uniform(0, 2 * np.pi, size=(nmodes, NA))
        self.envelope = envelope
        self.gen = np.real(1j * TADJ)  # real antisymmetric: (i T^e) = f^e
    def A(self, x):
        a = np.einsum('ma,ma->a', self.c, np.cos((self.q @ x)[:, None] + self.ph))
        if self.envelope is not None: a = a * self.envelope(x)
        return a
    def U(self, x):
        M = np.einsum('a,abc->bc', self.A(np.asarray(x, float)), self.gen)   # real antisymmetric
        lam, V = np.linalg.eigh(1j * M)                                      # iM Hermitian
        return np.real((V * np.exp(-1j * lam)) @ V.conj().T)


def D(Ua, Ub): return np.trace(Ua.T @ Ub)
def F(Ua, Ub, Uc):
    # f^{ijk} f^{lmn} Ua^{il} Ub^{jm} Uc^{kn}
    G = np.tensordot(FABC, Uc, axes=([2], [1]))          # l m k
    G = np.tensordot(G, Ub, axes=([1], [1]))             # l k j
    G = np.tensordot(G, Ua, axes=([0], [1]))             # k j i
    return np.sum(G.transpose(2, 1, 0) * FABC)
def Q(Ua, Ub, Uc, Ud):
    A = Ua.T @ Ub; B = Uc.T @ Ud
    X = np.tensordot(TADJ, A, axes=([2], [0]))           # e i k
    Y = np.tensordot(TADJ, B, axes=([2], [0]))           # e k l -> need sum_e X[e,i,k] Y[e,k,i]
    return np.real(np.einsum('eik,eki->', X, Y))


class Colour:
    """Evaluates the singlet invariants at named points; U cached per point (by coordinates)."""
    def __init__(self, field):
        self.field = field; self.cache = {}
    def u(self, x):
        key = (round(float(x[0]), 15), round(float(x[1]), 15))
        if key not in self.cache:
            if len(self.cache) > 20000: self.cache.clear()
            self.cache[key] = self.field.U(x)
        return self.cache[key]
    def D(self, a, b): return D(self.u(a), self.u(b))
    def F(self, a, b, c): return F(self.u(a), self.u(b), self.u(c))
    def Q(self, a, b, c, d): return Q(self.u(a), self.u(b), self.u(c), self.u(d))

# ---------------------------------------------------------------- one-variable functions
def E_fun(s):
    s = float(s)
    if abs(s) < 1e-14: return 0j
    si, ci = sici(abs(s))
    return ci - GE - np.log(abs(s)) - 1j * np.sign(s) * si

def J_fun(s):
    s = float(s)
    if abs(s) < 1e-8: return 1.0 - 0.5j * s
    return (1 - np.exp(-1j * s)) / (1j * s)

def Ehat(s):
    s = float(s)
    if abs(s) < 1e-300: s = 1e-300
    return np.exp(1j * s) * exp1(1j * s)

def Psi_fun(s):
    s = float(s)
    if abs(s) < 1e-14: return 0j
    return -GE - np.log(abs(s)) - 0.5j * np.pi * np.sign(s) - Ehat(s)

# ---------------------------------------------------------------- transverse kernels
def Kv(v): return v / (v @ v)
def Kf(a, b, z): return Kv(a - z) @ Kv(b - z)

def Mij(x, z, w):
    a = x - w; b = x - z; c = z - w; A = a @ a; B = b @ b; C = c @ c
    L = np.log(B / A); d = np.eye(2)
    t1 = L * (d / (2 * C) - np.outer(a, b) / (A * B)
              + (np.outer(a, c) * B / A + np.outer(c, b) * A / B - np.outer(c, a) - np.outer(b, c)) / (C * (A - B)))
    den = A * B - (a @ b) ** 2
    A0 = np.arccos(np.clip(a @ c / np.sqrt(A * C), -1, 1)) / np.sqrt(max(den, 1e-300))
    t2 = A0 / C * (np.outer(c, a) + np.outer(b, c) - (a @ b) * (np.outer(c, b) / B + np.outer(a, c) / A))
    t3 = -np.log(C / A) / (2 * C) * np.outer(c, b) / B
    return t1 + t2 + t3

def ell(z, w, x): return -0.5 * np.log(((z - w) @ (z - w)) / ((x - w) @ (x - w)))

def Tij(x, z, w):
    return (np.eye(2) / (2 * ((z - w) @ (z - w))) + np.outer(Kv(z - w), Kv(x - z))
            - np.outer(Kv(x - w), Kv(x - z)) - np.outer(Kv(x - w), Kv(z - w)))

# Group I, Row I.III (P = K(x'-w'))
def K_III_z(P):
    x, y, z, w, xp, wp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "w'"))
    Pv = Kv(xp - wp)
    return Pv @ Mij(x, z, w) @ Kv(y - z) + ell(z, w, x) * Kf(y, w, z) * (Pv @ Kv(x - w))

def K_III_z_Mhat_only(P):   # without the ell term (for bookkeeping)
    x, y, z, w, xp, wp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "w'"))
    return Kv(xp - wp) @ Mij(x, z, w) @ Kv(y - z)

def K_III_w(P):
    x, y, z, w, xp, wp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "w'"))
    return -E_fun(KVEC @ (z - w)) * Kf(x, w, z) * (Kv(xp - wp) @ Kv(y - w))

# Group I, Row I.V first part: x, z gluon positions, charges x', y, y'; s = k.(x-z)
def K_V1(P, part=None):
    x, z, y, yp, xp, wp = (P[k] for k in ('x', 'z', 'y', "y'", "x'", "w'"))
    Pv = Kv(xp - wp); s = KVEC @ (x - z)
    t1 = J_fun(s) * (Pv @ Kv(z - x)) * (Kv(y - z) @ Kv(yp - x))
    t2 = -E_fun(s) * (Pv @ Kv(y - z)) * (Kv(z - x) @ Kv(yp - x))
    t3 = -np.exp(-1j * s) * E_fun(-s) * (Pv @ Kv(yp - x)) * (Kv(z - x) @ Kv(y - z))
    return {None: t1 + t2 + t3, 'J': t1, 'E': t2, 'E-': t3}[part]

# Group II and III.II
def K_T(P, swap=False):
    x, y, z, xp = P['x'], P['y'], P['z'], P["x'"]
    w, wp = (P["w'"], P['w']) if swap else (P['w'], P["w'"])
    return (Kv(xp - wp) @ Tij(x, z, w) @ Kv(y - z)) * np.log(((x - z) @ (x - z)) / ((x - w) @ (x - w)))

# Group III (s = k.(z-y'))
def s3(P): return KVEC @ (P['z'] - P["y'"])
def K_Ehat(P):
    x, y, z, w, xp, yp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "y'"))
    return Ehat(s3(P)) * ((Kv(yp - z) @ Kv(xp - yp)) * (Kv(y - w) @ Kv(x - z))
                         - (Kv(yp - z) @ Kv(x - z)) * (Kv(xp - yp) @ Kv(y - w)))
def K_Psi(P):
    x, y, z, w, xp, yp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "y'"))
    return Psi_fun(s3(P)) * (Kv(yp - z) @ Kv(x - z)) * (Kv(xp - yp) @ Kv(y - w))
def g_delta_yp(P):   # coefficient of [pi delta(s) - i P/s] in K_delta^{(y')}
    x, y, z, w, xp, yp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "y'"))
    return (Kv(yp - z) @ Kv(y - w)) * (Kv(x - z) @ Kv(xp - yp))
def g_delta_z(P):    # coefficient of [pi delta(s) - i P/s] in K_delta^{(z)}
    x, y, z, w, xp, yp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "y'"))
    return -(Kv(yp - z) @ Kv(y - w)) * (Kv(x - z) @ Kv(xp - z))
