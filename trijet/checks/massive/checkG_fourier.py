# Check G: massive coordinate-space functions (antiquark channel; quark channel = same with zeta->eta).
#   zb*FT[{k^j p^n, k^j, p^n, 1}/((k^2+eps^2) den3)] = {R^jY^n G11, i R^j G10, -i Y^n G01, G00},  FT = int d2k d2p/(2pi)^2 e^{ik.R - ip.Y}
#   den3 = xi eta (k^2+eps^2) + zeta zb^2 (p^2+om^2),  eps^2 = zeta zb Q^2 + m^2,  om = xi m/zb,  a = xi eta/(zeta zb^2)
#   Feynman-parameter form:  G_ab = a^b/(2 zeta zb) int_0^1 dx x^{-1-b} (A_x/Delta_x)^{(a+b)/2} K_{a+b}(sqrt(A_x Delta_x))
#       A_x = zeta zb Q^2 + m^2 (1 + x zeta xi/eta),  Delta_x = R^2 + a Y^2/x
#   independent route: p-integral exactly, then 1D Hankel transform in k.
import numpy as np
from scipy import integrate, special
kv=special.kv
def pars(z,xi,Q,m):
    zb=1-z; eta=1-z-xi; a=xi*eta/(z*zb**2); eps=np.sqrt(z*zb*Q*Q+m*m); om=xi*m/zb
    E=np.sqrt(Q*Q+m*m*(1/z+1/eta)); return zb,eta,a,eps,om,E
def G_param(al,be,R,Y,z,xi,Q,m):
    zb,eta,a,eps,om,E=pars(z,xi,Q,m)
    def f(x):
        A=z*zb*Q*Q+m*m*(1+x*z*xi/eta); Dl=R*R+a*Y*Y/x; s=np.sqrt(A*Dl)
        return x**(-1-be)*(A/Dl)**((al+be)/2)*kv(al+be,s)
    # integrate in t = -ln x (integrand x f(x) dt), robust when the peak sits at x ~ a Y^2/R^2 << 1
    v,_=integrate.quad(lambda t: np.exp(-t)*f(np.exp(-t)),0,80,limit=800,epsabs=1e-15,epsrel=1e-12,
                       points=[max(0.0,-np.log(a*Y*Y/(R*R)))]); return a**be/(2*z*zb)*v
def G_hankel(al,be,R,Y,z,xi,Q,m):
    zb,eta,a,eps,om,E=pars(z,xi,Q,m)
    M=lambda k: np.sqrt(om*om+a*(k*k+eps*eps))
    pY = (lambda k: M(k)*kv(1,M(k)*Y)/Y) if be==1 else (lambda k: kv(0,M(k)*Y))
    kR = (lambda k: k*k*special.j1(k*R)/R) if al==1 else (lambda k: k*special.j0(k*R))
    v,_=integrate.quad(lambda k: kR(k)*pY(k)/(k*k+eps*eps),0,np.inf,limit=5000,epsabs=1e-15,epsrel=1e-11)
    return v/(z*zb)
def D_(R,Y,z,xi): zb=1-z; eta=1-z-xi; return np.sqrt(z*zb*R*R+xi*eta*Y*Y/zb)
pts=[(1.0,0.7,0.3,0.2,1.1,0.8),(0.5,1.3,0.25,0.35,0.9,1.4),(1.6,0.35,0.55,0.15,1.7,0.3)]
print("(1) Feynman-parameter form vs independent Hankel evaluation (rel. diff.)")
for (R,Y,z,xi,Q,m) in pts:
    for (al,be) in [(1,1),(1,0),(0,1),(0,0)]:
        p=G_param(al,be,R,Y,z,xi,Q,m); h=G_hankel(al,be,R,Y,z,xi,Q,m)
        print(f"   R={R} Y={Y} zeta={z} xi={xi} Q={Q} m={m}  G{al}{be}: {p:+.10e}  {h:+.10e}  {abs(p-h)/abs(h):.1e}")
print("(2) m=0:  G11 = Q K1(QD)/(D Y^2)  (user's function)")
for (R,Y,z,xi,Q,m) in pts:
    D=D_(R,Y,z,xi); p=G_param(1,1,R,Y,z,xi,Q,0.0); c=Q*kv(1,Q*D)/(D*Y*Y)
    print(f"   {p:.12e} {c:.12e}  {abs(p-c)/c:.1e}")
print("(3) instantaneous: zb*FT[1/den3] = E K1(E D)/D, E^2 = Q^2 + m^2(1/zeta+1/eta)")
for (R,Y,z,xi,Q,m) in pts:
    zb,eta,a,eps,om,E=pars(z,xi,Q,m); D=D_(R,Y,z,xi)
    M=lambda k: np.sqrt(om*om+a*(k*k+eps*eps))
    v,_=integrate.quad(lambda k: k*special.j0(k*R)*kv(0,M(k)*Y),0,np.inf,limit=5000,epsabs=1e-15,epsrel=1e-11)
    h=v/(z*zb); c=E*kv(1,E*D)/D
    print(f"   {h:.12e} {c:.12e}  {abs(h-c)/c:.1e}")
print("(4) UV (Y->0): Y^2 G11 -> eps K1(eps R)/(zeta zb R),  Y^2 G01 -> K0(eps R)/(zeta zb)  [= after-SW limits]")
R,z,xi,Q,m=1.0,0.3,0.2,1.1,0.8; zb,eta,a,eps,om,E=pars(z,xi,Q,m)
for Y in (1e-1,1e-2,1e-3):
    print(f"   Y={Y:g}: {Y*Y*G_param(1,1,R,Y,z,xi,Q,m):.8e} (lim {eps*kv(1,eps*R)/(z*zb*R):.8e})   "
          f"{Y*Y*G_param(0,1,R,Y,z,xi,Q,m):.8e} (lim {kv(0,eps*R)/(z*zb):.8e})")
print("(5) UV (Y->0): G10 - A10 and G00 - A00 stay finite (both ~ ln(1/Y) with the same coefficient)")
def A_after(al,be,R,Y,z,xi,Q,m):
    zb,eta,a,eps,om,E=pars(z,xi,Q,m)
    fR = eps*kv(1,eps*R)/R if al==1 else kv(0,eps*R)
    fY = om*kv(1,om*Y)/Y if be==1 else kv(0,om*Y)
    return fR*fY/(z*zb)
for Y in (1e-2,1e-3,1e-4,1e-5):
    print(f"   Y={Y:g}:  G10-A10 = {G_param(1,0,R,Y,z,xi,Q,m)-A_after(1,0,R,Y,z,xi,Q,m):+.8e}   "
          f"G00-A00 = {G_param(0,0,R,Y,z,xi,Q,m)-A_after(0,0,R,Y,z,xi,Q,m):+.8e}   (A10 = {A_after(1,0,R,Y,z,xi,Q,m):.3e})")
print("(6) after-SW 2D transforms by Hankel: int k^2 J1(kR)/(k^2+e^2) = e K1(eR),  int k J0(kR)/(k^2+e^2) = K0(eR)")
for (R_,e_) in [(0.7,1.3),(1.5,0.4)]:
    v1,_=integrate.quad(lambda k: k*k*special.j1(k*R_)/(k*k+e_*e_),0,np.inf,limit=5000,weight=None) if False else (None,None)
    # oscillatory tails: use the sine/cosine-weighted QAWF via scipy's 'weight' is not available for Bessel; split & sum
    def hankel(f,R_):
        tot=0.0; zeros=special.jn_zeros(1,4000)/R_ if True else None
        a0=0.0
        for b in zeros:
            v,_=integrate.quad(f,a0,b,limit=200); tot+=v; a0=b
        return tot
    h1=hankel(lambda k: (k*k/(k*k+e_*e_)-1.0)*special.j1(k*R_),R_) + 1.0/R_   # int_0^inf J1(kR) dk = 1/R
    h0=hankel(lambda k: k*special.j0(k*R_)/(k*k+e_*e_),R_)
    print(f"   R={R_} e={e_}:  {h1:.8f} vs {e_*kv(1,e_*R_):.8f}   {h0:.6f} vs {kv(0,e_*R_):.6f}")
