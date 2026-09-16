"""Grid refinement: does route A converge onto route B?"""
import numpy as np, mpmath as mp
mp.mp.dps=20
Kp=1.0; Pv=30.0; Lam=1e-9; kt=np.array([0.7,-0.4]); Rcut=300.0
def K(a,b): d=a-b; return d/(d@d)
x=np.array([0.31,-0.77]); xp=np.array([-0.52,0.19])
y=np.array([1.13,0.42]); yp=np.array([-0.31,0.88]); r=yp-y; kap=(kt@r)/Kp
X,Xp=K(x,y),K(xp,yp); XXp=X@Xp; r2=r@r
i0=complex(mp.quad(lambda t: mp.e**(-1j*kap*t),[Lam,Pv]))
i1=np.log(Pv/Lam)+complex(mp.quad(lambda t:(mp.e**(-1j*kap*t)-1)/t,[0,Pv]))
i2=complex(mp.quad(lambda t: t*mp.e**(-1j*kap*t),[Lam,Pv]))
i3=complex(mp.quad(lambda t: t*mp.e**(-1j*kap*t)/(Kp+t)**2,[Lam,Pv]))
T=np.pi/2*(np.log(Rcut**2/r2)+1)*np.eye(2)-np.pi*np.outer(r,r)/r2
B=(np.trace(T)*XXp*(i1+i2/Kp**2)+2*(X@T@Xp)*i3+2/Kp*(Xp@T@X-X@T@Xp)*i0)
print("route B (z in closed form, exact) = %s" % np.round(B,8))
print("%-16s %-30s %s" % ("z grid","route A (p+ first, then z)","rel. diff"))
for nr,nth in ((2000,360),(4000,720),(8000,1440),(16000,2880)):
    rr=np.geomspace(1e-8,Rcut,nr); th=np.linspace(0,2*np.pi,nth,endpoint=False)
    RR,TH=np.meshgrid(rr,th,indexing='ij')
    Z=np.stack([RR*np.cos(TH),RR*np.sin(TH)],-1)
    W=RR*np.gradient(rr)[:,None]*(2*np.pi/nth)
    U=y-Z; V=yp-Z; M=U/(U**2).sum(-1)[...,None]; N=V/(V**2).sum(-1)[...,None]
    dot=lambda a,b:(a*b).sum(-1)
    A=(((dot(M,N)*XXp*(i1+i2/Kp**2)) + 2*(M@X)*(N@Xp)*i3
        + 2/Kp*((M@Xp)*(N@X)-(M@X)*(N@Xp))*i0)*W).sum()
    print("%-16s %-30s %.2e" % ("%dx%d"%(nr,nth), np.round(A,6), abs(A-B)/abs(B)))
