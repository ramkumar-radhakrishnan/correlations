"""Validate the y'=y subtraction scheme in the toy model.

sigma(V)      = int d2C  br(C) f(C) FV(k.C, V)                      (original)
sigma_CT(V)   = int d2C  brF(Chat) f(0) chi(C) FV(k.C, V)           (counterterm)
sigma_reg(V)  = sigma - sigma_CT  (integrand subtracted pointwise)

brF(Chat) = A_i (delta_ij - 2 Chat_i Chat_j) X_j / (A^2 X^2)   [y'->y frozen]
chi(C)    = exp(-C^2/r0^2)

Analytic counterterm (V -> inf):
  int d2C e^{-iq.C} (delta_ij - 2 Chat_i Chat_j) chi(C)
      = -2 pi (delta_ij - 2 qhat_i qhat_j) T(|q|),
  T(q) = (2/q^2)(1 - e^{-q^2 r0^2/4}) - (r0^2/2) e^{-q^2 r0^2/4}
so
  sigma_CT = -2 pi f(0) [(d-2khat khat)_ij A_i X_j/(A^2 X^2)]
             * int_0^1 dxi (1-xi)/xi T(|k|/xi)
(no log V: T(k/xi) ~ 2 xi^2/k^2 at xi->0.)
"""
import numpy as np
from scipy.special import sici
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

A  = np.array([1.0, 0.2]);  X = np.array([-0.6, 0.8]);  sh = np.array([0.13, 0.21])
R  = 0.5; r0 = 0.1; kmag = 1.0
A2 = A @ A; X2 = X @ X
f0 = np.exp(-(sh @ sh))

def br(cx, cy):
    bx, by = X[0]-cx, X[1]-cy
    B2 = bx*bx + by*by; C2 = cx*cx + cy*cy
    if C2 == 0.0: return 0.0
    return (A[0]*bx + A[1]*by)/(A2*B2) - 2.0*(A[0]*cx+A[1]*cy)*(bx*cx+by*cy)/(A2*B2*C2)

def brF(cx, cy):
    C2 = cx*cx + cy*cy
    if C2 == 0.0: return 0.0
    n = np.array([cx, cy])/np.sqrt(C2)
    S = np.eye(2) - 2.0*np.outer(n, n)
    return A @ S @ X / (A2*X2)

def f(cx, cy):  return np.exp(-((cx-sh[0])**2 + (cy-sh[1])**2))
def chi(cx, cy): return np.exp(-(cx*cx+cy*cy)/r0**2)

def FV(a, V):
    if abs(a) < 1e-14: return np.log(V) - 1.0 + 1.0/V
    s1, c1 = sici(abs(a)); sV, cV = sici(abs(a)*V)
    g1 = (cV-c1) - 1j*np.sign(a)*(sV-s1)
    g2 = np.exp(-1j*a) - np.exp(-1j*a*V)/V - 1j*a*g1
    return g1 - g2

def chord(g, cx):
    yl = np.sqrt(max(R*R - cx*cx, 0.0))
    if yl == 0.0: return 0.0
    return quad(lambda cy: g(cx, cy), -yl, yl, limit=400)[0]

grid = np.linspace(-R, R, 4001)
Hfull = CubicSpline(grid, [chord(lambda x,y: br(x,y)*f(x,y), c) for c in grid])
Hct   = CubicSpline(grid, [chord(lambda x,y: brF(x,y)*f0*chi(x,y), c) for c in grid])
Hreg  = CubicSpline(grid, [chord(lambda x,y: br(x,y)*f(x,y) - brF(x,y)*f0*chi(x,y), c)
                           for c in grid])

def I1d(H, V):
    fre = quad(lambda a: FV(kmag*a, V).real*H(a), -R, R, limit=20000,
               points=[0.0], epsabs=1e-9, epsrel=1e-9)[0]
    fim = quad(lambda a: FV(kmag*a, V).imag*H(a), -R, R, limit=20000,
               points=[0.0], epsabs=1e-9, epsrel=1e-9)[0]
    return fre + 1j*fim

# sanity check of the J2 Fourier identity at one (q, component) point
q = 7.3; qv = np.array([q, 0.0])
def integ(comp):
    def gr(r):
        def gth(t):
            n = np.array([np.cos(t), np.sin(t)])
            S = np.eye(2) - 2*np.outer(n, n)
            return (np.cos(q*r*np.cos(t)))*S[comp]*np.exp(-r*r/r0**2)
        return quad(gth, 0, 2*np.pi, limit=200)[0]*r
    return quad(gr, 0, 6*r0, limit=200)[0]
T = lambda qq: (2/qq**2)*(1-np.exp(-qq**2*r0**2/4)) - (r0**2/2)*np.exp(-qq**2*r0**2/4)
num11, num22 = integ((0,0)), integ((1,1))
ana11 = -2*np.pi*T(q)*(1-2*1.0)   # (d-2qhat qhat)_11 with qhat=(1,0) -> -1
ana22 = -2*np.pi*T(q)*(1-0.0)     # _22 -> +1
print(f"FT identity check at q={q}: 11: {num11:+.8f} vs {ana11:+.8f}"
      f" | 22: {num22:+.8f} vs {ana22:+.8f}\n")

# analytic counterterm, V -> infinity
Gcontr = (-A[0]*X[0] + A[1]*X[1])/(A2*X2)        # (d-2khat khat)_ij A_i X_j /(A2 X2), k=(1,0)
Ixi = quad(lambda xi: (1-xi)/xi * T(kmag/xi), 0, 1, limit=2000)[0]
CT_analytic = -2*np.pi*f0*Gcontr*Ixi
print(f"counterterm, analytic (V=inf): {CT_analytic:+.8f}")
print(f"counterterm, numeric 2D (V=1e8): {I1d(Hct, 1e8):+.8f}\n")

print("subtracted integral sigma_reg(V)  [should be flat and finite]:")
for V in [1e2, 1e4, 1e6, 1e8]:
    print(f"   V={V:8.0e}   {I1d(Hreg, V):+.8f}")

tot = I1d(Hfull, 1e8); reg = I1d(Hreg, 1e8); ct = I1d(Hct, 1e8)
print(f"\nconsistency at V=1e8:  reg + CT = {reg+ct:+.8f}")
print(f"                       original = {tot:+.8f}")

print("\npointwise y'-density after xi-integration, V=1e12, C=(eps,0):")
print("   eps        unsubtracted              subtracted")
for eps in [1e-2, 1e-3, 1e-4, 1e-5]:
    u = br(eps,0)*f(eps,0)*FV(kmag*eps, 1e12)
    s = (br(eps,0)*f(eps,0) - brF(eps,0)*f0*chi(eps,0))*FV(kmag*eps, 1e12)
    print(f"   {eps:7.0e}   {u.real:+12.6f}{u.imag:+10.6f}i   {s.real:+12.8f}{s.imag:+11.8f}i")
