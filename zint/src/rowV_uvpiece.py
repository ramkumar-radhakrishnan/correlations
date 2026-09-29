"""Row V part (b): the user's "UV piece"

  d3N^UV_2b/d3k = 1/(2pi)^3 g^4/(8pi^4) 1/k+ Int e^{-ik.w'} e^{ik.x} Int dxi/(2pi) Int d^2r e^{iq.r}
                  1/(xib D) (s.P)/(P^2 r^2) C_UV(xi) [U(x')-U(w')] Nc[U(y)-U(x)] rho rho

Everything reduces to the scalar integral  K(xi) = mu^{-2eps} Int d^dr e^{iq.r} s^2/(xib D r^2).
Checks: prefactor, exact decomposition, UV pole, Bessel closed form vs direct 2D quadrature,
endpoint behaviour in xi, the xi integral with the plus prescription, and the r->s artefact.
"""
import numpy as np, sympy as sp
from scipy import integrate, special

EG = float(sp.EulerGamma)
hr = lambda t: print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# ----------------------------------------------------------------------------------------------
hr("1. PREFACTOR AND DECOMPOSITION")
g, kp = sp.symbols('g k_plus', positive=True)
user = 1/(2*sp.pi)**3 * g**4/(8*sp.pi**4) / kp / (2*sp.pi)
N = 1/(2*sp.pi)**3 * g**4/(16*sp.pi**4) / kp
print("   user prefactor / (N/pi) =", sp.simplify(user/(N/sp.pi)))
xi = sp.symbols('xi', positive=True); xb = 1 - xi
r1, r2, s1, s2_ = sp.symbols('r1 r2 s1 s2', real=True)
rr = r1**2 + r2**2; ss = s1**2 + s2_**2; rs = r1*s1 + r2*s2_
D = (r1 - s1)**2 + (r2 - s2_)**2 + xi/xb*ss
print("   D0 = D(r=0) = s^2/xib :", sp.simplify(D.subs({r1: 0, r2: 0}) - ss/xb) == 0)
lhs = ss/(xb*D*rr); rhs = 1/rr + 2*rs/(rr*D) - 1/D
print("   s^2/(xib D r^2) = 1/r^2 + 2 r.s/(r^2 D) - 1/D :", sp.simplify(lhs - rhs) == 0)

# ----------------------------------------------------------------------------------------------
hr("2. TENSOR vs SCALAR UV TERM (d_perp = d)")
print("   <rhat^i rhat^m>_d = delta/d, so d_perp xi xib (P.r)(s.r)/r^4 -> (d_perp/d) xi xib (P.s)/r^2 = xi xib (P.s)/r^2")
print("   => the user's scalar xi xib term has the same pole; the traceless difference is UV finite")
print("      and belongs to the rest of 2b (it may be evaluated at d_perp = 2 in d = 2).")

# ----------------------------------------------------------------------------------------------
s = np.array([0.7, -0.4]); k = np.array([1.3, 0.5])
s2 = s @ s; kn = np.linalg.norm(k); ks = k @ s
ell = np.log(kn**2*s2*np.exp(2*EG)/4)

def Aterm(x):
    """A/pi = (1/pi) Int d^2r e^{iq.r} 2 r.s/(r^2 D) via Feynman parameter + Bessel"""
    b = 1 - x; q = b*k; qn = b*kn; qs = b*ks
    def f(t, part):
        Del = t*s2*(1 - t + x/b); a = np.sqrt(Del)
        v = np.exp(1j*t*qs)*(2j*qs*special.k0(qn*a) + 2*t*s2*qn/a*special.k1(qn*a))
        return v.real if part == 0 else v.imag
    pts = [1 - 10**(-j) for j in range(1, 9)]
    return sum((1j if p else 1)*integrate.quad(lambda t: f(t, p), 0, 1, points=pts, limit=800,
                                                 epsabs=1e-12, epsrel=1e-11)[0] for p in (0, 1))

def Bterm(x):
    """B/pi = (1/pi) Int d^2r e^{iq.r}/D = 2 e^{iq.s} K0(qM)"""
    b = 1 - x; M = np.sqrt(x/b*s2)
    return 2*np.exp(1j*b*ks)*special.k0(b*kn*M)

def G(x):
    """finite part in coordinate form: K/pi = 1/eps_UV + log(mubar^2 s^2 e^{2g}/4) + G"""
    b = 1 - x
    return -np.log(b**2*kn**2*s2*np.exp(2*EG)/4) + Aterm(x) - Bterm(x)

hr("3. CLOSED FORM vs DIRECT 2D QUADRATURE (absolutely convergent integrand)")
print("   Ktilde = Int d^2r (1/r^2)[e^{iq.r} s^2/(xib D) - s^2/(r^2+s^2)]  should equal  pi*G(xi)")
def Ktilde_direct(x):
    b = 1 - x; q = b*k; M2 = x/b*s2
    def f(rho, ph, part):
        r = rho*np.array([np.cos(ph), np.sin(ph)])
        D = (r - s) @ (r - s) + M2
        v = (np.exp(1j*(q @ r))*s2/(b*D) - s2/(rho**2 + s2))/rho   # d^2r/r^2 = drho dphi / rho
        return v.real if part == 0 else v.imag
    out = 0
    for p in (0, 1):
        tot = 0
        for (a, c) in [(0, 0.5), (0.5, 1.5), (1.5, 6), (6, 60), (60, 2000)]:
            tot += integrate.dblquad(lambda ph, rho: f(rho, ph, p), a, c, 0, 2*np.pi,
                                     epsabs=1e-10, epsrel=1e-10)[0]
        out += (1j if p else 1)*tot
    return out
for x in (0.3, 0.6):
    d = Ktilde_direct(x); c = np.pi*G(x)
    print(f"   xi={x}:  direct {d.real:.8f}{d.imag:+.8f}i   formula {c.real:.8f}{c.imag:+.8f}i")

# ----------------------------------------------------------------------------------------------
hr("4. ENDPOINTS OF G")
def areg():
    """regular part of A/pi at xi = 0 (q = k, Delta = t(1-t)s^2), 2e^{ik.s}/(1-t) subtracted"""
    sn = np.sqrt(s2)
    def f(t, part):
        a = sn*np.sqrt(t*(1 - t))
        v = np.exp(1j*t*ks)*(2j*ks*special.k0(kn*a) + 2*t*s2*kn/a*special.k1(kn*a)) - 2*np.exp(1j*ks)/(1 - t)
        return v.real if part == 0 else v.imag
    pts = [1 - 10**(-j) for j in range(1, 9)]
    return sum((1j if p else 1)*integrate.quad(lambda t: f(t, p), 0, 1, points=pts, limit=800,
                                                 epsabs=1e-12, epsrel=1e-11)[0] for p in (0, 1))
c0 = areg() - ell*(1 - np.exp(1j*ks))
print(f"   predicted: G -> e^(ik.s) log(1/xi) + c0 ,  e^(ik.s) = {np.exp(1j*ks):.8f} ,  c0 = {c0:.8f}")
for x in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    v = G(x); pr = np.exp(1j*ks)*np.log(1/x) + c0
    print(f"   xi={x:.0e}:  G = {v:.8f}   pred = {pr:.8f}   diff = {abs(v - pr):.2e}")
print("   predicted: G -> log(1/xib) + 0   as xi -> 1")
for b in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    v = G(1 - b); pr = np.log(1/b)
    print(f"   xib={b:.0e}: G = {v:.8f}   pred = {pr:.8f}   diff = {abs(v - pr):.2e}")

# ----------------------------------------------------------------------------------------------
hr("5. THE xi INTEGRAL WITH THE PLUS PRESCRIPTION")
C = lambda x: x*(1 - x) + x/(1 - x) + (1 - x)/x
eks = np.exp(1j*ks)
def sub(x):
    return C(x)*G(x) - (eks*np.log(1/x) + c0)/x - np.log(1/(1 - x))/(1 - x)
def cquad(fun, a, b, pts=None):
    re = integrate.quad(lambda t: fun(t).real, a, b, points=pts, limit=800, epsabs=1e-10, epsrel=1e-10)[0]
    im = integrate.quad(lambda t: fun(t).imag, a, b, points=pts, limit=800, epsabs=1e-10, epsrel=1e-10)[0]
    return re + 1j*im
# work in u = log variables near the ends for accuracy
def cquad_log(fun, lam):
    # [lam, 1/2]: xi = e^{-u};  [1/2, 1-lam]: xib = e^{-u}
    f1 = lambda u: fun(np.exp(-u))*np.exp(-u)
    f2 = lambda u: fun(1 - np.exp(-u))*np.exp(-u)
    return cquad(f1, np.log(2), np.log(1/lam)) + cquad(f2, np.log(2), np.log(1/lam))
H = cquad_log(sub, 1e-12)
print(f"   H = Int_0^1 [C_UV G - (e^(ik.s)log(1/xi)+c0)/xi - log(1/xib)/xib] = {H:.8f}")
for lam in (1e-3, 1e-4, 1e-5):
    L = np.log(1/lam)
    direct = cquad_log(lambda x: C(x)*G(x), lam)
    pred = 0.5*(1 + eks)*L**2 + c0*L + H
    print(f"   lambda={lam:.0e}: Int C_UV G = {direct:.6f}   (1+e^iks)L^2/2 + c0 L + H = {pred:.6f}   diff {abs(direct - pred):.1e}")
print("   UV coefficient: Int C_UV = 2L - 11/6 (sympy):",
      sp.limit(sp.integrate(xi*xb + xi/xb + xb/xi, (xi, sp.Symbol('l'), 1 - sp.Symbol('l'))) + 2*sp.log(sp.Symbol('l')),
               sp.Symbol('l'), 0))

# ----------------------------------------------------------------------------------------------
hr("6. THE r -> s ARTEFACT AT xi -> 0 AND ITS CANCELLATION")
P1, P2 = sp.symbols('P1 P2', real=True)
Pr = P1*r1 + P2*r2; Ps = P1*s1 + P2*s2_; PP = P1**2 + P2**2
T4 = xb/xi*(Ps - Pr)/(PP*rr)*(1 - rs/(2*ss))
at_s = {r1: s1, r2: s2_}
print("   full block T4 (the only 1/xi block) at r = s :", sp.simplify(T4.subs(at_s)))
print("   user's (xib/xi)(P.s)/(P^2 r^2) at r = s      :", sp.simplify((xb/xi*Ps/(PP*rr)).subs(at_s)))
print("   => the user's piece has 1/D ~ 1/((r-s)^2+M^2) with a non-zero residue at r = s:")
print("      Int d^2r near s gives pi e^{iq.s} log(1/M^2) -> pi e^{ik.s} log(1/xi); the full T4 vanishes there.")
print("      The rest of 2b therefore carries -W e^{ik.s} log(1/xi)/xi, and the (1/2) e^{ik.s} L^2 cancels in the sum.")
