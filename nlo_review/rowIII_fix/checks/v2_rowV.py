# Group III Row V, the -1/p+ term:  z-integral  J(kap,b) = int d^2z e^{-i kap.(y'-z)} K(y'-z).K(x-z),  b=x-y', kap = k p+/k+
# then  calJ(Lam) = int_{Lam}^{inf} dp+/p+ J   and its small-Lam expansion  pi (L+lam')^2 + j0,
# L = ln(k+/Lam), lam' = ln(2 e^{-gamma}/(|k||b|)).
import numpy as np, mpmath as mp
mp.mp.dps=25
g=mp.euler
def J(kap,b,th):
    # kap=|kappa|, b=|b|, th = angle between kappa and b
    c=kap*b*mp.cos(th)
    def f(al):
        M=kap*mp.sqrt(al*(1-al))
        return mp.e**(1j*al*c)*(mp.pi*b*mp.besselk(1,b*M)/M-mp.pi/M**2)
    I=mp.quad(f,[0,0.5,1])
    Ein=lambda z: mp.quad(lambda t:(1-mp.e**(-t*z))/t,[0,1])   # Ein(z)=int_0^1 (1-e^{-tz})/t dt
    A=-Ein(-1j*c)-mp.e**(1j*c)*Ein(1j*c)
    return -0.5*(kap**2*I+mp.pi*A-2*mp.pi*(1+mp.e**(1j*c))*mp.log(2*mp.e**(-g)/(kap*b)))
# direct check of the formula for J with Gaussian damping, one point
def Jdirect(kap,b,th,eta):
    kv=np.array([kap*np.cos(th),kap*np.sin(th)]); bv=np.array([b,0.0])
    from scipy import integrate
    def re_im(part):
        def inner(r):
            def gph(ph):
                rv=np.array([r*np.cos(ph),r*np.sin(ph)]); d=rv-bv
                val=np.exp(1j*kv@rv)*(rv@d)/((rv@rv)*(d@d))*np.exp(-eta*r*r)
                return val.real if part==0 else val.imag
            v,_=integrate.quad(gph,0,2*np.pi,limit=400,points=[0.0])
            return v*r
        tot=0
        for lo,hi in [(0,b),(b,3*b),(3*b,40)]:
            v,_=integrate.quad(inner,lo,hi,limit=400); tot+=v
        return tot
    return re_im(0)+1j*re_im(1)
print("formula J(kap=1.3,b=0.8,th=0.9) =",J(1.3,0.8,0.9))
for eta in [0.05,0.02,0.01]:
    print("  direct with damping eta=",eta,":",Jdirect(1.3,0.8,0.9,eta))
print("small-kappa: J - 2 pi ln(2e^-g/(kap b)) :",[mp.nstr(J(kk,1.0,0.7)-2*mp.pi*mp.log(2*mp.e**(-g)/kk),8) for kk in (1e-2,1e-3,1e-4)])
print("large-kappa: J :",[mp.nstr(J(kk,1.0,0.7),8) for kk in (10,30,100)])
