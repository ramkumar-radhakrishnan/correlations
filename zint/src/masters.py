"""alpha-parameter (Feynman) representations of every master integral needed for
   W = \int d^2 r e^{iq.r} T(r)/D ,   D = (r-s)^2 + M^2 .

Conventions (all 2-vectors):  A = r^2 + m^2 , B = D ,
   alpha A + (1-alpha) B = (r-b)^2 + Dl ,  b = (1-alpha) s ,
   Dl(alpha) = (1-alpha)(alpha s^2 + M^2) + alpha m^2 .
Basic 2D transforms with lam = |q| sqrt(Dl):
   U1 = 2 pi K0,  U2 = pi q K1/sqrt(Dl),  U3 = (pi q/2)[K1/Dl^{3/2} + q K0/(2 Dl)]
   V2^i = i pi q^i K0,  V3^i = i pi q^i q K1/(4 sqrt(Dl))
   T2^{ij} = pi[d^{ij} K0 - q^i q^j sqrt(Dl) K1/q]
   T3^{ij} = (pi/4)[d^{ij} q K1/sqrt(Dl) - q^i q^j K0]
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import k0, k1

def _K(lam):
    return k0(lam), k1(lam)

class Kin:
    def __init__(self, s, q, M2, m2=0.0):
        self.s, self.q, self.M2, self.m2 = np.asarray(s), np.asarray(q), M2, m2
        self.s2 = float(s @ s); self.qn = float(np.hypot(*q)); self.qs = float(q @ s)
        self.D0 = self.s2 + M2
    def Dl(self, a):
        return (1-a)*(a*self.s2 + self.M2) + a*self.m2

def _cint(f, a=0.0, b=1.0):
    re = quad(lambda t: f(t).real, a, b, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
    im = quad(lambda t: f(t).imag, a, b, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
    return re + 1j*im

# ---------------------------------------------------------------- masters (r^2, D)
def G_scalar(K):                     # \int e^{iqr}/((r^2+m^2) D)     [UV divergent]
    def f(a):
        Dl = K.Dl(a); lam = K.qn*np.sqrt(Dl); K0v, K1v = _K(lam)
        return np.exp(1j*(1-a)*K.qs)*np.pi*K.qn*K1v/np.sqrt(Dl)
    return _cint(f)

def Psi(K):                          # \int e^{iqr} r^i/(r^2 D)
    def f(a):
        Dl = K.Dl(a); lam = K.qn*np.sqrt(Dl); K0v, K1v = _K(lam)
        U2 = np.pi*K.qn*K1v/np.sqrt(Dl); V2 = 1j*np.pi*K.q*K0v
        return np.exp(1j*(1-a)*K.qs)*((1-a)*K.s*U2 + V2)
    return np.array([_cint(lambda a, i=i: f(a)[i]) for i in range(2)])

def Nten(K):                         # \int e^{iqr} r^i r^j/(r^2 D)
    def f(a):
        Dl = K.Dl(a); lam = K.qn*np.sqrt(Dl); K0v, K1v = _K(lam)
        U2 = np.pi*K.qn*K1v/np.sqrt(Dl); V2 = 1j*np.pi*K.q*K0v
        T2 = np.pi*(np.eye(2)*K0v - np.outer(K.q, K.q)*np.sqrt(Dl)*K1v/K.qn)
        b = (1-a)*K.s
        return np.exp(1j*(1-a)*K.qs)*(np.outer(b, b)*U2 + np.outer(b, V2)
                                      + np.outer(V2, b) + T2)
    return np.array([[_cint(lambda a, i=i, j=j: f(a)[i, j]) for j in range(2)]
                     for i in range(2)])

def Tten(K):                         # \int e^{iqr}[r^i r^j - d^{ij} r^2/2]/(r^4 D)
    """Traceless projection of M^(2)[r^i r^j]; UV finite, so m = 0 throughout."""
    def X(a):
        Dl = K.Dl(a); lam = K.qn*np.sqrt(Dl); K0v, K1v = _K(lam)
        U3 = np.pi*K.qn/2*(K1v/Dl**1.5 + K.qn*K0v/(2*Dl))
        V3 = 1j*np.pi*K.q*K.qn*K1v/(4*np.sqrt(Dl))
        T3 = np.pi/4*(np.eye(2)*K.qn*K1v/np.sqrt(Dl) - np.outer(K.q, K.q)*K0v)
        b = (1-a)*K.s
        M = 2*a*(np.outer(b, b)*U3 + np.outer(b, V3) + np.outer(V3, b) + T3)
        M = M - np.eye(2)*np.trace(M)/2                      # traceless projection
        return np.exp(1j*(1-a)*K.qs)*M
    return np.array([[_cint(lambda a, i=i, j=j: X(a)[i, j]) for j in range(2)]
                     for i in range(2)])

def Phi(K):                          # \int e^{iqr}/D  (closed form)
    return 2*np.pi*np.exp(1j*K.qs)*k0(np.sqrt(K.M2)*K.qn)

# ------------------------------------------------- masters (r^2, (r-s)^2), massless
def _A_common(K, a):
    Dl = a*(1-a)*K.s2; lam = K.qn*np.sqrt(Dl); K0v, K1v = _K(lam)
    U1 = 2*np.pi*K0v; U2 = np.pi*K.qn*K1v/np.sqrt(Dl); V2 = 1j*np.pi*K.q*K0v
    T2 = np.pi*(np.eye(2)*K0v - np.outer(K.q, K.q)*np.sqrt(Dl)*K1v/K.qn)
    return Dl, U1, U2, V2, T2, np.exp(1j*(1-a)*K.qs)

def A7(K):                           # \int e^{iqr} r.(r-s)/(r^2 (r-s)^2)
    def f(a):
        Dl, U1, U2, V2, T2, E = _A_common(K, a)
        return E*(U1 - 2*Dl*U2 + (1-2*a)*(K.s @ V2))
    return _cint(f)

def A8(K, P):                        # \int e^{iqr} (P.r)(s.(r-s))/(r^2 (r-s)^2)
    def f(a):
        Dl, U1, U2, V2, T2, E = _A_common(K, a)
        return E*(P @ T2 @ K.s - a*K.s2*(P @ V2) + (1-a)*(P @ K.s)*(K.s @ V2)
                  - a*(1-a)*K.s2*(P @ K.s)*U2)
    return _cint(f)

def A9(K, P):                        # \int e^{iqr} (s.r)(P.(r-s))/(r^2 (r-s)^2)
    def f(a):
        Dl, U1, U2, V2, T2, E = _A_common(K, a)
        return E*(K.s @ T2 @ P - a*(P @ K.s)*(K.s @ V2) + (1-a)*K.s2*(P @ V2)
                  - a*(1-a)*K.s2*(P @ K.s)*U2)
    return _cint(f)
