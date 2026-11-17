import numpy as np, mpmath as mp, sys
from scipy import integrate, special
g=np.euler_gamma
def J(kap,th,b=1.0):
    c=kap*b*np.cos(th)
    def f(al,part):
        M=kap*np.sqrt(al*(1-al))
        if b*M<1e-6:
            core=np.pi*(b*b/2)*(np.log(b*M/2)+g-0.5)   # b K1(bM)/M - 1/M^2 for small M
        else:
            core=np.pi*b*special.k1(b*M)/M-np.pi/M**2
        v=np.exp(1j*al*c)*core
        return v.real if part==0 else v.imag
    I=sum(integrate.quad(f,lo,hi,args=(p,),limit=200,epsabs=1e-12,epsrel=1e-11)[0]*(1 if p==0 else 1j) for p in (0,1) for lo,hi in ((0,.5),(.5,1)))
    def Ein(z):
        fr=lambda t:((1-np.exp(-t*z))/t).real; fi=lambda t:((1-np.exp(-t*z))/t).imag
        return integrate.quad(fr,0,1,limit=200)[0]+1j*integrate.quad(fi,0,1,limit=200)[0]
    A=-Ein(-1j*c)-np.exp(1j*c)*Ein(1j*c)
    return -0.5*(kap**2*I+np.pi*A-2*np.pi*(1+np.exp(1j*c))*np.log(2*np.exp(-g)/(kap*b)))
k0=2*np.exp(-g)
def j0(th,Kmax=400.0):
    def integrand(t,part):
        kap=np.exp(t); v=J(kap,th)-(2*np.pi*np.log(k0/kap) if kap<k0 else 0.0)
        return v.real if part==0 else v.imag
    out=0
    edges=[-30,-12,-6,-2,np.log(k0)]+list(np.log(np.linspace(k0,Kmax,60))[1:])
    for p in (0,1):
        s=0
        for lo,hi in zip(edges[:-1],edges[1:]):
            s+=integrate.quad(integrand,lo,hi,args=(p,),limit=200,epsabs=1e-10,epsrel=1e-9)[0]
        out+=s*(1 if p==0 else 1j)
    return out
for th in [0.0,0.4,0.9,np.pi/2,2.2,np.pi]:
    a=j0(th,200.0); b_=j0(th,400.0)
    print(f"theta={th:.4f}  j0(Kmax=200)={a:.6f}  j0(Kmax=400)={b_:.6f}")
