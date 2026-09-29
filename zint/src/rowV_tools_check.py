"""Numerical checks of every 'tool' formula used in rowV_2b_step_by_step.pdf (Row V, 2b UV piece)."""
import numpy as np
from scipy import integrate, special
import mpmath as mp
mp.mp.dps = 20

EG = np.euler_gamma
print("T1  polar area element: Int_{rho<1} d^2r = pi :",
      integrate.dblquad(lambda ph, rho: rho, 0, 1, 0, 2*np.pi)[0], np.pi)

print("T2  Gamma: Gamma(z+1)=z Gamma(z) at z=0.37 :", special.gamma(1.37), 0.37*special.gamma(0.37))
x = 1e-4
print("    Gamma(1+x) ~ 1 - gamma_E x  (x=1e-4):", special.gamma(1 + x), 1 - EG*x)
for a in (0.2, 0.5, 0.8):
    lhs = integrate.quad(lambda t: t**(a - 1)/(1 + t), 0, 1)[0] + integrate.quad(lambda t: t**(a - 1)/(1 + t), 1, np.inf)[0]
    print(f"    Beta: Int_0^oo t^(a-1)/(1+t) dt at a={a}: {lhs:.10f}   Gamma(a)Gamma(1-a) = {special.gamma(a)*special.gamma(1-a):.10f}"
          f"   pi/sin(pi a) = {np.pi/np.sin(np.pi*a):.10f}")

for z in (0.3, 2.0):
    ang = integrate.quad(lambda ph: np.cos(z*np.cos(ph)), 0, 2*np.pi)[0]
    print(f"T3  Int_0^2pi dphi e^(i z cos phi) at z={z}: {ang:.10f}   2 pi J0(z) = {2*np.pi*special.j0(z):.10f}")
for q, m in ((1.3, 0.4), (0.5, 2.0)):
    val = float(mp.quadosc(lambda u: u*mp.besselj(0, q*u)/(u**2 + m**2), [0, mp.inf], zeros=lambda n: mp.besseljzero(0, n)/q))
    print(f"    Int_0^oo u J0(qu)/(u^2+m^2) du (q={q}, m={m}): {val:.8f}   K0(qm) = {special.k0(q*m):.8f}")
for z in (1e-3, 1e-5):
    print(f"    K0({z}) = {special.k0(z):.8f}   -log(z/2)-gamma_E = {-np.log(z/2)-EG:.8f} ;"
          f"  z K1(z) = {z*special.k1(z):.8f}")
h = 1e-6; z = 0.7
print("    K1 = -dK0/dz at z=0.7:", special.k1(z), -(special.k0(z + h) - special.k0(z - h))/(2*h))

a, b = 2.3, 0.6
print("T5  Feynman: 1/(ab) =", 1/(a*b), "  Int_0^1 dt/[ta+(1-t)b]^2 =",
      integrate.quad(lambda t: 1/(t*a + (1 - t)*b)**2, 0, 1)[0])

for lam in (1e-3, 1e-6):
    f = lambda xi: np.cos(xi) + 2*xi
    direct = integrate.quad(lambda xi: f(xi)/xi, lam, 1, points=[10*lam], limit=400)[0]
    plus = integrate.quad(lambda xi: (f(xi) - f(0))/xi, 0, 1)[0] + f(0)*np.log(1/lam)
    print(f"T6  plus prescription, f=cos+2xi, lambda={lam:.0e}: direct {direct:.8f}   [.]_+ + f(0)L {plus:.8f}")
    direct = integrate.quad(lambda xi: np.log(1/xi)/xi, lam, 1, points=[10*lam], limit=400)[0]
    print(f"    Int_lambda^1 log(1/xi)/xi = {direct:.8f}   L^2/2 = {0.5*np.log(1/lam)**2:.8f}")

for aa in (0.3, 1.0, 3.0):
    val = float(mp.quad(lambda t: (mp.besselj(0, aa*t) - 1/(1 + t**2))/t, [0, 1])
                + mp.quadosc(lambda t: mp.besselj(0, aa*t)/t, [1, mp.inf], zeros=lambda n: mp.besseljzero(0, n)/aa)
                - mp.quad(lambda t: 1/(t*(1 + t**2)), [1, mp.inf]))
    print(f"R1  Int_0^oo dx/x [J0(ax)-1/(1+x^2)] at a={aa}: {val:.7f}   -1/2 log(a^2 e^(2g)/4) = {-0.5*np.log(aa**2*np.exp(2*EG)/4):.7f}")
