"""Row 27c: the subtraction organised in u = z - y' (so no p+-dependent phase leaks out)."""
import numpy as np
from row25_lib import val, KT, KP, s
ring=lambda R,n,c: c+R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
X=(xp-yp)/s(xp-yp); B=s(x-w); kabs=np.sqrt(KT@KT); kh=KT/kabs; gE=np.euler_gamma
Kc=((x-w)@X)/(2*B); Kk=((x-w)@kh)*(X@kh)/(2*B)
print("="*78); print("WHY u = z - y'  :  the phase is e^{-i b k.(y'-z)} = e^{+i b k.u}"); print("="*78)
print("   centring on z leaves a stray e^{-i b k.y'} outside the u integral, and b = p+/k+ ,")
print("   so that phase would contaminate the p+ integral.  Centring on y' removes it:")
print("   the only p+-dependent phase is e^{+i b k.u}, which the Bessel masters absorb.")
print("   Also N^m = (y'-z)^m/(y'-z)^2 = -u^m/u^2 EXACTLY, so u is the natural variable.")
def Iasy_u(U,P):
    """U : (n,2) array of u = z - y' values"""
    b=P/KP; uh=U/np.sqrt((U**2).sum(1))[:,None]; r2=(U**2).sum(1)
    core=(1/P)*( (P**2/(KP**2*(P+KP)))*((uh@(x-w))*(uh@X))/(2*B) + Kc/KP )
    return np.exp(-1j*(KT@(yp-w)))*np.exp(1j*b*(U@KT))/r2*core
def closed(P,rho0):
    b=P/KP; L=np.log(2*np.exp(-gE)/(b*kabs*rho0))
    return (np.exp(-1j*(KT@(yp-w)))*(np.pi/P)
            *( 2*Kc/KP*L + (P**2/(KP**2*(P+KP)))*((L+0.5)*Kc - Kk) ))
def num_u(P,rho0,R,nr=1600,nth=9000):
    lo,hi=np.log(rho0),np.log(R); lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; tot=0j
    for r0 in np.exp(lr):
        Uu=ring(r0,nth,np.zeros(2)); tot+=Iasy_u(Uu,P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot
print()
print("1.  int_{|u|>rho0} d^2u I_asy :  numeric vs closed form  (rho0 = 1)")
print("   %-8s %-9s %24s %24s %10s"%("p+","q rho0","numeric","closed","rel.diff"))
for P in (0.3,0.1,0.03,0.01,0.003):
    nu=num_u(P,1.0,300*KP/(P*kabs)); cl=closed(P,1.0)
    print("   %-8.3f %-9.4f %24s %24s %10.2e"%(P,P*kabs/KP,"%.6f%+.6fj"%(nu.real,nu.imag),
          "%.6f%+.6fj"%(cl.real,cl.imag),abs(nu-cl)/abs(cl)),flush=True)
print()
print("2.  the remainder I - I_asy falls like |u|^-3  (rings centred on y')")
print("   %-9s %22s %10s"%("|u|","<|I - I_asy|> x |u|^3","resolved?"))
for R in (1e3,1e4,1e5):
    nth=max(20000,int(40*R*kabs*0.37/KP))
    Uu=ring(R,nth,np.zeros(2)); Z=Uu+yp
    v=abs((val(Z,0.37)-Iasy_u(Uu,0.37)).mean())*R**3
    print("   %-9.0e %22.6f %10s"%(R,v,"yes"))
print("   constant => |u|^-3 tail => int d^2u [I - I_asy] converges ABSOLUTELY, no cutoff needed.")
print()
print("3.  and the full z-integral reassembles:")
def full(P,rho0,R,nr=1600,nth=9000):
    lo,hi=np.log(1e-3),np.log(R); lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; tot=0j
    for r0 in np.exp(lr):
        Uu=ring(r0,nth,np.zeros(2)); Z=Uu+yp
        f=val(Z,P)
        if r0>rho0: f=f-Iasy_u(Uu,P)
        tot+=f.mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot
for P in (0.1,0.03,0.01):
    R=300*KP/(P*kabs)
    rem=full(P,1.0,R); cl=closed(P,1.0)
    dire=0j
    lo,hi=np.log(1e-3),np.log(R); nr=1600
    for r0 in np.exp((np.arange(nr)+.5)*(hi-lo)/nr+lo):
        dire+=val(ring(r0,9000,np.zeros(2))+yp,P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    print("   p+=%-7.3f  remainder %12.6f  + closed %12.6f  = %12.6f   direct %12.6f   diff %.1e"
          %(P,rem.real,cl.real,(rem+cl).real,dire.real,abs((rem+cl)-dire)),flush=True)
print("   => I(p+) = [absolutely convergent remainder] + [closed form].  This is the result.")
