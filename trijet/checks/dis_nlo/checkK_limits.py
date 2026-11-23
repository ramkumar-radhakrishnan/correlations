# Check K: numerical limits of the radial functions used in the paper.
#   G_ab = a^b/(2 zeta zb) int_0^1 dx x^{-1-b} (A_x/Delta_x)^{(a+b)/2} K_{a+b}(sqrt(A_x Delta_x))
#   (1) massless limit: G11 -> Q K1(QD)/(D Y^2); the helicity-flip terms m*G01, xi*om*G10, m*xi*om*G00 -> 0
#   (2) soft limit xi -> 0 at fixed R, Y:  Y^2 G11 -> eps K1(eps R)/(zeta zb R),  Y^2 G01 -> K0(eps R)/(zeta zb)
#   (3) photoproduction limit Q -> 0 is finite for m > 0
#   (4) after-shock-wave emission kernel: om K1(om Y)/Y -> 1/Y^2 for om Y << 1 (dead cone in coordinate space)
import numpy as np
from scipy import integrate, special
kv = special.kv
def pars(z, xi, Q, m):
    zb = 1 - z; eta = 1 - z - xi; a = xi*eta/(z*zb**2)
    return zb, eta, a, np.sqrt(z*zb*Q*Q + m*m), xi*m/zb
def G(al, be, R, Y, z, xi, Q, m):
    zb, eta, a, eps, om = pars(z, xi, Q, m)
    def f(x):
        A = z*zb*Q*Q + m*m*(1 + x*z*xi/eta); Dl = R*R + a*Y*Y/x
        return x**(-1 - be)*(A/Dl)**((al + be)/2)*kv(al + be, np.sqrt(A*Dl))
    t0 = max(0.0, -np.log(a*Y*Y/(R*R)))
    v, _ = integrate.quad(lambda t: np.exp(-t)*f(np.exp(-t)), 0, 200, limit=2000, epsabs=0, epsrel=1e-11, points=[t0])
    return a**be/(2*z*zb)*v
def D(R, Y, z, xi): zb = 1 - z; eta = 1 - z - xi; return np.sqrt(z*zb*R*R + xi*eta*Y*Y/zb)

R, Y, z, xi, Q = 1.0, 0.7, 0.3, 0.2, 1.1
print("(1) massless limit (R=1, Y=0.7, zeta=0.3, xi=0.2, Q=1.1)")
ref = Q*kv(1, Q*D(R, Y, z, xi))/(D(R, Y, z, xi)*Y*Y)
for m in (0.3, 0.1, 0.03, 0.01, 0.0):
    zb, eta, a, eps, om = pars(z, xi, Q, m)
    g11 = G(1, 1, R, Y, z, xi, Q, m)
    print(f"   m={m:5.2f}: G11/[Q K1(QD)/(D Y^2)] = {g11/ref:.8f}   m*Y*G01 = {m*Y*G(0, 1, R, Y, z, xi, Q, m):.3e}"
          f"   xi*om*R*G10 = {xi*om*R*G(1, 0, R, Y, z, xi, Q, m):.3e}   m*xi*om*G00 = {m*xi*om*G(0, 0, R, Y, z, xi, Q, m):.3e}")
print("(2) soft limit (R=1, Y=0.7, zeta=0.3, Q=1.1, m=0.8)")
m = 0.8
for xi_ in (1e-1, 1e-2, 1e-3, 1e-4):
    zb, eta, a, eps, om = pars(z, xi_, Q, m)
    print(f"   xi={xi_:7.0e}: Y^2 G11 / [eps K1(eps R)/(zeta zb R)] = {Y*Y*G(1, 1, R, Y, z, xi_, Q, m)/(eps*kv(1, eps*R)/(z*zb*R)):.6f}"
          f"   Y^2 G01 / [K0(eps R)/(zeta zb)] = {Y*Y*G(0, 1, R, Y, z, xi_, Q, m)/(kv(0, eps*R)/(z*zb)):.6f}")
print("(3) photoproduction limit Q -> 0 (m=0.8, xi=0.2)")
for Q_ in (1.0, 0.3, 0.1, 0.0):
    print(f"   Q={Q_:4.1f}: " + "  ".join(f"G{al}{be}={G(al, be, R, Y, z, xi, Q_, m):.6e}" for al, be in ((1, 1), (1, 0), (0, 1), (0, 0))))
print("(4) dead cone in coordinate space: Y^2 * om K1(om Y)/Y")
for t in (0.01, 0.1, 1.0, 3.0, 10.0):
    print(f"   om*Y={t:5.2f}:  {t*kv(1, t):.6f}")
