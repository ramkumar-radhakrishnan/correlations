"""The single master of the UV-projected row:
      G(xi) = int d^2r e^{i q.r} / [ (r^2+m^2) ( (r-s)^2 + M^2 ) ] ,
   with M^2 = (xi/xibar) s^2, q = xibar k, D0 = s^2 + M^2 = s^2/xibar.

   Exact split:   G = (pi/D0) [ ln(D0/m^2) + Ihat ] ,
      Ihat = D0 * int_0^1 dt [ e^{i t q.s} |q| K1(|q| sqrt(D)) / sqrt(D) - 1/(t D0) ],
      D    = t [ (1-t) s^2 + M^2 ] .
   The integrand is integrable at t -> 0; Ihat is finite and carries no m."""
import numpy as np
from scipy.integrate import quad
from scipy.special import k0, k1

def Ihat(s, k, xi):
    s2 = float(s@s); xb = 1-xi; M2 = xi/xb*s2; D0 = s2/xb
    q = xb*np.asarray(k); qn = float(np.hypot(*q)); qs = float(q@s)
    def f(t):
        D = t*((1-t)*s2 + M2)
        lam = qn*np.sqrt(D)
        core = qn*k1(lam)/np.sqrt(D) if lam < 700 else 0.0
        return np.exp(1j*t*qs)*core - 1/(t*D0)
    re = quad(lambda t: f(t).real, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    im = quad(lambda t: f(t).imag, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return D0*(re + 1j*im)

def G(s, k, xi, m):
    s2 = float(s@s); xb = 1-xi; D0 = s2/xb
    return np.pi/D0*(np.log(D0/m**2) + Ihat(s, k, xi))

if __name__ == "__main__":
    import sys; sys.path.insert(0,'.')
    from quad2d import integrate
    s = np.array([0.7,-0.4]); k = np.array([0.9,0.6])
    print("check:  G  (alpha/t representation)  vs  direct 2D quadrature")
    for xi in (0.3, 0.65, 0.1):
        xb=1-xi; s2=s@s; M2=xi/xb*s2; q=xb*k
        for m in (1e-2, 1e-3):
            an = G(s,k,xi,m)
            num = integrate(lambda r: 1/((np.sum(r*r,-1)+m*m)*(np.sum((r-s)**2,-1)+M2)),
                            q, split=(np.hypot(*s),), ntheta=2048, nlev=18)
            print(f"  xi={xi:<5} m={m:.0e}: rep {an:+.9f}   quad {num:+.9f}   rel {abs(an-num)/abs(num):.1e}")
    print()
    print("Ihat(xi) for these kinematics (s.k = %.3f):" % (s@k))
    for xi in (0.9,0.7,0.5,0.3,0.1,0.03,0.01):
        print(f"   xi={xi:<6} Ihat = {Ihat(s,k,xi):+.6f}")
