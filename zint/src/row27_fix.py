"""Row 27b: validity of the closed form, and the subleading p+ primitive."""
import numpy as np, sympy as sp, mpmath as mp
from row25_lib import val, KT, KP, s
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
X=(xp-yp)/s(xp-yp); B=s(x-w); kabs=np.sqrt(KT@KT); kh=KT/kabs; gE=np.euler_gamma
Kc=((x-w)@X)/(2*B); Kk=((x-w)@kh)*(X@kh)/(2*B)
def Iasy(Z,P):
    b=P/KP; zh=Z/np.sqrt((Z**2).sum(1))[:,None]; r2=(Z**2).sum(1)
    core=(1/P)*( (P**2/(KP**2*(P+KP)))*((zh@(x-w))*(zh@X))/(2*B) + Kc/KP )
    return np.exp(-1j*(KT@(yp-w)))*np.exp(-1j*b*(KT@yp))*np.exp(1j*b*(Z@KT))/r2*core
def asy_closed(P,rho0):
    b=P/KP; L=np.log(2*np.exp(-gE)/(b*kabs*rho0))
    return (np.exp(-1j*(KT@(yp-w)))*np.exp(-1j*b*(KT@yp))
            *(np.pi/P)*( 2*Kc/KP*L + (P**2/(KP**2*(P+KP)))*((L+0.5)*Kc - Kk) ))
def asy_num(P,rho0,R,nr=1600,nth=9000):
    lo,hi=np.log(rho0),np.log(R); lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; tot=0j
    for r0 in np.exp(lr): tot+=Iasy(ring(r0,nth),P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot
print("="*78); print("1.  THE CLOSED FORM IS THE SMALL-(q rho0) LIMIT -- and q rho0 = p+|k|rho0/k+"); print("="*78)
print("   %-8s %-9s %22s %22s %10s"%("p+","q rho0","numeric","closed form","rel.diff"))
for P in (1.5,0.37,0.1,0.03,0.01,0.003):
    rho0=1.0; q=P*kabs/KP
    nu=asy_num(P,rho0,300*KP/(P*kabs)); cl=asy_closed(P,rho0)
    print("   %-8.3f %-9.4f %22s %22s %10.2e"%(P,q*rho0,"%.6f%+.6fj"%(nu.real,nu.imag),
          "%.6f%+.6fj"%(cl.real,cl.imag),abs(nu-cl)/abs(cl)),flush=True)
print("   -> the error is O((q rho0)^2) and vanishes exactly where the logs live (p+ -> 0).")
print("      If you want it valid at LARGE p+ too, keep the masters unexpanded:")
print("         M_0 = 2 pi int_{q rho0}^inf dt J_0(t)/t   (exact),  and similarly M_2.")
print()
print("="*78); print("2.  THE REMAINDER FALLS LIKE |z|^-3 -> ABSOLUTELY CONVERGENT"); print("="*78)
print("   %-9s %20s"%("|z|","<|I - I_asy|> x |z|^3"))
for R in (1e3,1e4,1e5,1e6):
    Z=ring(R,20000)
    print("   %-9.0e %20.6f"%(R,abs((val(Z,0.37)-Iasy(Z,0.37)).mean())*R**3))
print("   constant => remainder ~ |z|^-3 , and int d^2z |z|^-3 converges.  No cutoff needed.")
print()
print("="*78); print("3.  THE p+ INTEGRALS"); print("="*78)
p,k_,Xi=sp.symbols('p kplus Xi',positive=True)
F1=-sp.log(Xi/p)**2/2
print("   leading:  int dp (1/p) log(Xi/p) = -(1/2) log^2(Xi/p)")
print("     check:  d/dp[-(1/2)log^2(Xi/p)] - (1/p)log(Xi/p) = %s"%sp.simplify(sp.diff(F1,p)-sp.log(Xi/p)/p))
F2=p*(sp.log(Xi/p)+sp.Rational(3,2)) - k_*((sp.log(Xi)+sp.Rational(1,2))*sp.log(p+k_)
     - (sp.log(p)*sp.log(1+p/k_)+sp.polylog(2,-p/k_)))
print("   subleading:  int dp [p/(p+k+)] (log(Xi/p) + 1/2)")
print("     = p(log(Xi/p) + 3/2) - k+[(log Xi + 1/2) log(p+k+) - log p log(1+p/k+) - Li_2(-p/k+)]")
print("     check:  d/dp - integrand = %s"%sp.simplify(sp.diff(F2,p)-p/(p+k_)*(sp.log(Xi/p)+sp.Rational(1,2))))
mp.mp.dps=20
kv,Xv,lo,hi=mp.mpf(1),mp.mpf(7),mp.mpf('1e-4'),mp.mpf('3')
num=mp.quad(lambda t: t/(t+kv)*(mp.log(Xv/t)+mp.mpf(1)/2),[lo,hi])
Fn=lambda t: t*(mp.log(Xv/t)+mp.mpf(3)/2)-kv*((mp.log(Xv)+mp.mpf(1)/2)*mp.log(t+kv)
              -(mp.log(t)*mp.log(1+t/kv)+mp.polylog(2,-t/kv)))
print("     numeric %s   primitive difference %s   diff %.1e"%(mp.nstr(num,12),mp.nstr(Fn(hi)-Fn(lo),12),abs(num-(Fn(hi)-Fn(lo)))))
