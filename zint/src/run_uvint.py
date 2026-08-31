"""The two integrations of the UV-projected row, and the UV/finite separation.
   -> ../RESULTS_uvint.txt"""
import numpy as np, mpmath as mp
from scipy.integrate import quad
from uv_master import Ihat, G
from masters import Kin, Psi, Phi
from quad2d import integrate
out=[]
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)

s=np.array([0.7,-0.4]); k=np.array([0.9,0.6]); P=np.array([1.3,0.5])
s2=s@s; P2=P@P; sP=s@P; ph=np.exp(1j*(k@s)); g=0.5772156649015329
C=lambda x: x*(1-x)+x/(1-x)+(1-x)/x
p(f"kinematics: s = y-x = {s}   k = {k}   P = x'-w' = {P}")
p(f"            s^2 = {s2:.4f}   k.s = {k@s:.4f}   e^(i k.s) = {ph:+.6f}   (s.P)/(P^2 s^2) = {sP/(P2*s2):.6f}")
p("")
p("A. the z integral.  G = int d^2r e^{iq.r}/[(r^2+m^2)((r-s)^2+M^2)] = (pi/D0)[ln(D0/m^2) + Ihat]")
p("   (i) representation vs direct 2D quadrature")
for xi in (0.3,0.65):
    xb=1-xi; M2=xi/xb*s2; q=xb*k; m=1e-3
    an=G(s,k,xi,m)
    num=integrate(lambda r: 1/((np.sum(r*r,-1)+m*m)*(np.sum((r-s)**2,-1)+M2)), q,
                  split=(np.hypot(*s),10*m,100*m), ntheta=2048, nlev=18)
    p(f"      xi={xi:.2f} m={m:.0e}:  rep {an:+.8f}   quad {num:+.8f}   rel {abs(an-num)/abs(num):.1e}")
p("   (ii) Ihat vs the independent route  ln(4 e^-2g /(q^2 D0)) + (2 s.Psi - Phi)/pi")
for xi in (0.65,0.3,0.1):
    xb=1-xi; M2=xi/xb*s2; D0=s2/xb; q=xb*k; qn=np.hypot(*q); K=Kin(s,q,M2)
    alt=np.log(4*np.exp(-2*g)/(qn**2*D0))+(2*(s@Psi(K))-Phi(K))/np.pi; a=Ihat(s,k,xi)
    p(f"      xi={xi:.2f}: {a:+.10f}   {alt:+.10f}   rel {abs(a-alt)/abs(alt):.1e}")
p("")
p("B. Ihat(xi), the finite part left by the z integral")
for xi in (0.9,0.7,0.5,0.3,0.1,0.01):
    p(f"      xi={xi:<5} Ihat = {Ihat(s,k,xi):+.6f}")
p("   xi -> 0 :  Ihat = e^(i k.s) ln(1/xi) + c0 + O(xi ln xi)")
for xi in (1e-3,1e-4,1e-5):
    p(f"      xi={xi:.0e}:  Ihat - e^(iks) ln(1/xi) = {Ihat(s,k,xi)-ph*np.log(1/xi):+.6f}")
p("   xi -> 1 :  Ihat ~ xibar ln(1/xibar) -> 0, so C_UV*Ihat stays integrable (no extra l_k)")
for xb in (1e-2,1e-3,1e-4):
    a=Ihat(s,k,1-xb); p(f"      xibar={xb:.0e}:  Ihat = {a:+.4e}   C_UV*Ihat = {C(1-xb)*a:+.6f}")
p("")
p("C. the xi integrals of C_UV = Pgg/(2Nc)   (mpmath, 30 digits, lambda = 1e-9)")
mp.mp.dps=30; lam=mp.mpf('1e-9'); lk=mp.log(1/lam); Cm=lambda t: t*(1-t)+t/(1-t)+(1-t)/t
p(f"      int C_UV            = {mp.nstr(mp.quad(Cm,[lam,1-lam]),14)}     2 l_k - 11/6             = {mp.nstr(2*lk-mp.mpf(11)/6,14)}")
p(f"      int C_UV ln(1/xib)  = {mp.nstr(mp.quad(lambda t:Cm(t)*mp.log(1/(1-t)),[lam,1-lam]),14)}     l_k^2/2 + pi^2/6 - 67/36 = {mp.nstr(lk**2/2+mp.pi**2/6-mp.mpf(67)/36,14)}")
p("")
p("D. the xi integral of the finite part,  int dxi C_UV(xi) Ihat(xi)")
def IC(l):
    f=lambda x: C(x)*Ihat(s,k,x)
    return (quad(lambda x:f(x).real,l,1-l,limit=200,epsabs=1e-10,epsrel=1e-9)[0]
            +1j*quad(lambda x:f(x).imag,l,1-l,limit=200,epsabs=1e-10,epsrel=1e-9)[0])
lams=[1e-3,1e-4,1e-5,1e-6]; vals=[IC(l) for l in lams]
for l,v in zip(lams,vals): p(f"      lambda={l:.0e}  l_k={np.log(1/l):7.4f} :  {v:+.6f}")
A=np.array([[np.log(1/l)**2/2,np.log(1/l),1.0] for l in lams])
c,*_=np.linalg.lstsq(A,np.array(vals),rcond=None)
c0=Ihat(s,k,1e-5)-ph*np.log(1e5)
p("      fit to  a l_k^2/2 + b l_k + c :")
p(f"         a = {c[0]:+.6f}    e^(i k.s) = {ph:+.6f}    (a must be e^(i k.s))")
p(f"         b = {c[1]:+.6f}    c0        = {c0:+.6f}    (b must be c0)")
p(f"         c = {c[2]:+.6f}")
open('/home/user/correlations/zint/RESULTS_uvint.txt','w').write("\n".join(out)+"\n")
