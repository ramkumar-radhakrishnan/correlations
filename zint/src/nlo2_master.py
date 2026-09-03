"""The single tensor master that the y' (or z) integration reduces to:

      B^{mj}(Q,c) = int d^2r  e^{-i Q.r}  (r^m/r^2) ((c-r)^j/(c-r)^2)

  Feynman:  1/(r^2 (c-r)^2) = int_0^1 da / [ (r - ab c)^2 + a ab c^2 ]^2 ,  ab = 1-a
  shift r = ab c + u :  r^m = u^m + ab c^m ,  (c-r)^j = a c^j - u^j ,  phase e^{-i ab Q.c}

      B^{mj} = int_0^1 da e^{-i ab Q.c} [ -T2^{mj} + a c^j V2^m - ab c^m V2^j
                                          + a ab c^m c^j U2 ]
      U2 = pi|Q| K1(lam)/sqrt(Del),  V2^i = -i pi Q^i K0(lam),
      T2^{ij} = pi[ d^{ij} K0(lam) - Q^i Q^j sqrt(Del) K1(lam)/|Q| ],
      Del = a ab c^2 ,  lam = |Q| sqrt(Del).
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import k0, k1

def Bten(Q, c):
    Qn = float(np.hypot(*Q)); c2 = float(c@c); Qc = float(Q@c)
    def f(a):
        ab = 1-a; Del = a*ab*c2; lam = Qn*np.sqrt(Del)
        K0v, K1v = k0(lam), k1(lam)
        U2 = np.pi*Qn*K1v/np.sqrt(Del)
        V2 = -1j*np.pi*Q*K0v
        T2 = np.pi*(np.eye(2)*K0v - np.outer(Q, Q)*np.sqrt(Del)*K1v/Qn)
        return np.exp(-1j*ab*Qc)*(-T2 + a*np.outer(V2, c) - ab*np.outer(c, V2)
                                  + a*ab*np.outer(c, c)*U2)
    out = np.zeros((2, 2), dtype=complex)
    for m in range(2):
        for j in range(2):
            re = quad(lambda a: f(a)[m, j].real, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
            im = quad(lambda a: f(a)[m, j].imag, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
            out[m, j] = re + 1j*im
    return out

if __name__ == "__main__":
    import sys; sys.path.insert(0, '.')
    from quad2d import integrate
    Q = np.array([0.9, 0.6]); c = np.array([0.7, -0.4])
    B = Bten(Q, c)
    print("B^{mj}(Q,c)   alpha-representation  vs  direct 2D quadrature")
    for (m, j) in [(0,0), (0,1), (1,0), (1,1)]:
        num = integrate(lambda r, m=m, j=j:
                        r[..., m]*(c-r)[..., j]/(np.sum(r*r,-1)*np.sum((c-r)**2,-1)),
                        -Q, split=(np.hypot(*c),), ntheta=4096, nlev=20, R0=12.0)
        print(f"   B^{m}{j}: rep {B[m,j]:+.9f}   quad {num:+.9f}   rel {abs(B[m,j]-num)/abs(num):.1e}")
    print(f"   trace  {np.trace(B):+.9f}")
