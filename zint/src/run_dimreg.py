"""Checks for the dim-reg / MSbar note.  -> ../RESULTS_dimreg.txt"""
import mpmath as mp
mp.mp.dps = 30
out=[]
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)
g = mp.euler

p("A. the d-dimensional transform, formula checked where the answer is known (d=3, a=1)")
d,a = 3,1
pref = mp.pi**(mp.mpf(d)/2)*2**(d-2*a)*mp.gamma((mp.mpf(d)-2*a)/2)/mp.gamma(a)
p(f"   pi^(d/2) 2^(d-2a) Gamma((d-2a)/2)/Gamma(a) = {mp.nstr(pref,14)}")
p(f"   known result  int d^3x e^(iqx)/x^2 = 2 pi^2/|q| ,  2 pi^2 = {mp.nstr(2*mp.pi**2,14)}")
p("")
p("B. the epsilon expansion:  pi Gamma(-eps) (q^2/(4 pi mu^2))^eps + pi/eps  ->  pi log(mubar^2/q^2)")
q2, mu2 = mp.mpf('2.7'), mp.mpf('1.3'); mubar2 = 4*mp.pi*mu2*mp.exp(-g)
p(f"   q^2 = {q2},  mu^2 = {mu2},  mubar^2 = 4 pi mu^2 e^-gamma = {mp.nstr(mubar2,12)}")
for e in ['1e-3','1e-4','1e-5','1e-6','1e-7']:
    eps = mp.mpf(e)
    v = mp.pi*mp.gamma(-eps)*(q2/(4*mp.pi*mu2))**eps + mp.pi/eps
    p(f"   eps={e:>5}:  Gamma-form + pi/eps = {mp.nstr(v,14):>16}   pi log(mubar^2/q^2) = {mp.nstr(mp.pi*mp.log(mubar2/q2),14)}")
p("")
p("C. the p+ moments   (lambda = 1e-9, l_k = log(1/lambda))")
lam = mp.mpf('1e-9'); lk = mp.log(1/lam)
C = lambda t: t*(1-t) + t/(1-t) + (1-t)/t
I0 = mp.quad(C,[lam,1-lam])
I1 = mp.quad(lambda t: C(t)*mp.log(1/(1-t)**2),[lam,1-lam])
Iodd = mp.quad(lambda t: C(t)*mp.log(t/(1-t)),[lam,1-lam])
Ic10 = mp.quad(lambda t: C(t)*mp.log(t*(1-t)),[lam,1-lam])
p(f"   int C_UV                = {mp.nstr(I0,14)}    2 l_k - 11/6            = {mp.nstr(2*lk-mp.mpf(11)/6,14)}")
p(f"   int C_UV log(1/xib^2)   = {mp.nstr(I1,14)}    l_k^2 + pi^2/3 - 67/18  = {mp.nstr(lk**2+mp.pi**2/3-mp.mpf(67)/18,14)}")
p(f"   int C_UV log(xi/xib)    = {mp.nstr(Iodd,8)}  (odd under xi <-> xibar, vanishes)")
p(f"   int C_UV log(xi xib)    = {mp.nstr(Ic10,14)}   [paper (C.10)];  sum with the row above: {mp.nstr(Ic10+I1,8)}")
open('/home/user/correlations/zint/RESULTS_dimreg.txt','w').write("\n".join(out)+"\n")

# ---------------------------------------------------------------- momentum space
out.append("")
out.append("D. MS-bar in momentum space: the 1/r^2 term as the one-loop transverse bubble")
cG = lambda e: mp.gamma(1+e)*mp.gamma(1-e)**2/mp.gamma(1-2*e)
out.append("   Gamma(-e)^2/Gamma(-2e) = -(2/e) Gamma(1-e)^2/Gamma(1-2e) :")
for e in ['0.3','0.05','-0.02']:
    e=mp.mpf(e); l=mp.gamma(-e)**2/mp.gamma(-2*e); r=-(2/e)*mp.gamma(1-e)**2/mp.gamma(1-2*e)
    out.append(f"      e={str(e):>6}: {mp.nstr(l,16):>22} {mp.nstr(r,16):>22}   diff {mp.nstr(l-r,4)}")
out.append("   the two routes (q^2=2.7, mu^2=1.3):")
q2,mu2 = mp.mpf('2.7'), mp.mpf('1.3'); mubar2 = 4*mp.pi*mu2*mp.exp(-g)
for e in ['1e-3','1e-5','1e-7']:
    e=mp.mpf(e)
    Tm=(mp.pi/e)*cG(e)*(4*mp.pi*mu2/q2)**e
    Tc=mp.pi*mp.gamma(-e)*(q2/(4*mp.pi*mu2))**e
    out.append(f"      e={str(e):>5}:  momentum finite part {mp.nstr(Tm-mp.pi/e,14):>17}"
               f"   coordinate finite part {mp.nstr(Tc+mp.pi/e,14):>17}")
    out.append(f"              pi log(mubar^2/q^2) = {mp.nstr(mp.pi*mp.log(mubar2/q2),14)}"
               f"    difference of routes {mp.nstr(Tm-Tc,14)}  vs 2 pi/e = {mp.nstr(2*mp.pi/e,14)}")
open('/home/user/correlations/zint/RESULTS_dimreg.txt','w').write("\n".join(out)+"\n")

# ------------------------------------------------ the remainder, and the assembly
import sys; sys.path.insert(0,'/home/user/correlations/zint/src')
import numpy as np
from scipy.integrate import quad as squad
from remainder import R as Rrem, Itilde
from uv_master import Ihat
from quad2d import integrate
gE = 0.5772156649015329
s = np.array([0.7,-0.4]); k = np.array([0.9,0.6]); s2v = float(s@s); ph = np.exp(1j*(k@s))
mu2 = 1.3
out.append("")
out.append("E. the remainder  R = int d^2r e^(iqr)(1/r^2)[1/D - 1/D0] = (2 s.Psi - Phi)/D0")
out.append("   (UV finite, no regulator)   analytic vs direct 2D quadrature:")
for xi in (0.3,0.65,0.1):
    xb=1-xi; M2=xi/xb*s2v; D0=s2v/xb; q=xb*k
    an=Rrem(s,q,M2,D0)
    num=integrate(lambda r:(1/(np.sum((r-s)**2,-1)+M2)-1/D0)/np.sum(r*r,-1),q,
                  split=(np.hypot(*s),),ntheta=2048,nlev=18)
    out.append(f"      xi={xi:<5} {an:+.9f}   {num:+.9f}   rel {abs(an-num)/abs(num):.1e}")
out.append("")
out.append("F. adding the two pieces: the q^2 cancels")
out.append("   T|MSbar/D0 + R  =  (pi/D0)[ log(e^(2 gamma) mu^2 D0 / 4) + Ihat ]")
for xi in (0.65,0.3,0.1):
    xb=1-xi; M2=xi/xb*s2v; D0=s2v/xb; q=xb*k; qn=np.hypot(*q)
    lhs=np.pi*np.log(mu2/qn**2)/D0 + Rrem(s,q,M2,D0)
    rhs=np.pi/D0*(np.log(np.exp(2*gE)*mu2*D0/4)+Ihat(s,k,xi))
    out.append(f"      xi={xi:<5} {lhs:+.10f}   {rhs:+.10f}   rel {abs(lhs-rhs)/abs(rhs):.1e}")
out.append("")
out.append("G. the assembled xi-integral,  int dxi C_UV [ log(e^2g mu^2 s^2/(4 xibar)) + Ihat ]")
Cf = lambda x: x*(1-x)+x/(1-x)+(1-x)/x
full = lambda xi: np.log(np.exp(2*gE)*mu2*s2v/(4*(1-xi))) + Ihat(s,k,xi)
def Iq(l):
    f=lambda x: Cf(x)*full(x)
    return (squad(lambda x:f(x).real,l,1-l,limit=200,epsabs=1e-10,epsrel=1e-9)[0]
            +1j*squad(lambda x:f(x).imag,l,1-l,limit=200,epsabs=1e-10,epsrel=1e-9)[0])
lams=[1e-3,1e-4,1e-5,1e-6]; vals=[Iq(l) for l in lams]
Amat=np.array([[np.log(1/l)**2,np.log(1/l),1.0] for l in lams])
cf,*_=np.linalg.lstsq(Amat,np.array(vals),rcond=None)
LL=np.log(np.exp(2*gE)*mu2*s2v/4)
out.append(f"      l_k^2 : {cf[0]:+.6f}   predicted (1 + e^(i k.(y-x)))/2 = {(1+ph)/2:+.6f}")
out.append(f"      l_k   : {cf[1]:+.6f}   predicted 2 L + c0             = {2*LL+(-0.588515-0.116237j):+.6f}")
out.append(f"      const : {cf[2]:+.6f}   predicted -11/6 L + pi^2/6 - 67/36 + c1 = "
           f"{-11/6*LL+np.pi**2/6-67/36+(0.350743+0.149676j):+.6f}")
out.append(f"      [ L = log(e^(2gamma) mu^2 (y-x)^2/4) = log(mu^2 (y-x)^2) + 2(gamma - log 2) ]")
out.append(f"      2(gamma - log 2) = {2*(gE-np.log(2)):+.8f}   <-- the 2 b (gamma - log 2) of K_JSJ, Eq. (2.64)")
open('/home/user/correlations/zint/RESULTS_dimreg.txt','w').write("\n".join(out)+"\n")
