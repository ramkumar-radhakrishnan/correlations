"""The second piece of the exact split:

     R = int d^2r e^{iq.r} (1/r^2) [ 1/D - 1/D0 ]
       = (1/D0) int d^2r e^{iq.r} (2 r.s - r^2)/(r^2 D)
       = (2 s.Psi - Phi)/D0                                (UV finite: no regulator)

  Psi^i = int d^2r e^{iq.r} r^i/(r^2 D)
        = pi int_0^1 da e^{i(1-a) q.s} [ (1-a) s^i |q| K1(lam)/sqrt(Del) + i q^i K0(lam) ]
  Phi   = int d^2r e^{iq.r}/D = 2 pi e^{iq.s} K0(M|q|)
  Del   = (1-a)(a s^2 + M^2),  lam = |q| sqrt(Del)
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import k0, k1

def sPsi(s, q, M2):
    """s . Psi"""
    s2 = float(s@s); qn = float(np.hypot(*q)); qs = float(q@s)
    def f(a):
        D = (1-a)*(a*s2 + M2); lam = qn*np.sqrt(D)
        return np.exp(1j*(1-a)*qs)*((1-a)*s2*qn*k1(lam)/np.sqrt(D) + 1j*qs*k0(lam))
    re = quad(lambda a: f(a).real, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    im = quad(lambda a: f(a).imag, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return np.pi*(re + 1j*im)

def Phi(s, q, M2):
    return 2*np.pi*np.exp(1j*(q@s))*k0(np.sqrt(M2)*np.hypot(*q))

def R(s, q, M2, D0):
    return (2*sPsi(s, q, M2) - Phi(s, q, M2))/D0

def Itilde(s, q, M2):
    """R = pi * Itilde / D0 ;  Itilde = (2 s.Psi - Phi)/pi"""
    return (2*sPsi(s, q, M2) - Phi(s, q, M2))/np.pi

if __name__ == "__main__":
    import sys; sys.path.insert(0, '.')
    from quad2d import integrate
    from uv_master import Ihat
    g = 0.5772156649015329
    s = np.array([0.7, -0.4]); k = np.array([0.9, 0.6]); s2 = s@s
    print("A. R  vs  direct 2D quadrature of e^{iqr}(1/r^2)[1/D - 1/D0]   (no regulator needed)")
    for xi in (0.3, 0.65, 0.1):
        xb = 1-xi; M2 = xi/xb*s2; D0 = s2/xb; q = xb*k
        an = R(s, q, M2, D0)
        num = integrate(lambda r: (1/(np.sum((r-s)**2,-1)+M2) - 1/D0)/np.sum(r*r,-1),
                        q, split=(np.hypot(*s),), ntheta=2048, nlev=18)
        print(f"   xi={xi:<5} analytic {an:+.9f}   quad {num:+.9f}   rel {abs(an-num)/abs(num):.1e}")
    print()
    print("B. Itilde = (2 s.Psi - Phi)/pi   vs   Ihat + log(e^{2gamma} q^2 D0/4)")
    for xi in (0.65, 0.3, 0.1):
        xb = 1-xi; M2 = xi/xb*s2; D0 = s2/xb; q = xb*k; qn = np.hypot(*q)
        a1 = Itilde(s, q, M2)
        a2 = Ihat(s, k, xi) + np.log(np.exp(2*g)*qn**2*D0/4)
        print(f"   xi={xi:<5} {a1:+.10f}   {a2:+.10f}   rel {abs(a1-a2)/abs(a2):.1e}")
