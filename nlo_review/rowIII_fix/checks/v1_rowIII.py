# Check of the Lambda->0 expansion of the exact Group I Row III pair (from Group_I.tex, N3a 2nd term + N3b last term)
#   P(Lam) = int d^2z (1/v^2) ln{ [k v + Lam (x-z)]^2 / (Lam^2 a^2) } F(z),  v=z-w, a=|x-w|, k+=1
# against
#   2 pi L^2 F(w) + 2L I1 + I2 + (pi^3/6) F(w),   L=ln(1/Lam)
#   I1 = int [F(z) - a^2/(a^2+v^2) F(w)]/v^2 ,  I2 = int ln(v^2/a^2) [F(z) - a^2/(a^2+v^2) F(w)]/v^2
import numpy as np
from scipy import integrate
w=np.array([0.0,0.0]); x=np.array([0.7,-0.4]); y=np.array([-0.5,0.9])
a=np.linalg.norm(x-w)
def F(z):  # generic smooth test function, decaying like 1/|z| (as K(y-z)[U(z)-U(y)] does)
    d=z-y; r2=d@d
    return (d[0]*0.8+d[1]*0.3)/(r2+0.05)*(1.0+0.6*np.exp(-((z[0]-0.2)**2+(z[1]+0.1)**2)/0.5))
Fw=F(w)
def polar(f, rmax=np.inf):
    # integral d^2 v f(v) in polar coords around w, log-radial variable
    def inner(t):
        r=np.exp(t)
        g=lambda ph: f(np.array([r*np.cos(ph), r*np.sin(ph)]))
        val,_=integrate.quad(g,0,2*np.pi,limit=400,epsabs=1e-12,epsrel=1e-11)
        return val*r*r
    val,_=integrate.quad(inner,-40,14,limit=800,epsabs=1e-11,epsrel=1e-10, points=[np.log(a)])
    return val
I1=polar(lambda v: (F(w+v)-a*a/(a*a+v@v)*Fw)/(v@v))
I2=polar(lambda v: np.log((v@v)/a**2)*(F(w+v)-a*a/(a*a+v@v)*Fw)/(v@v))
print("F(w)=",Fw," I1=",I1," I2=",I2)
for Lam in [1e-2,1e-3,1e-4,1e-5,1e-6]:
    def f(v):
        z=w+v; q=v+Lam*(x-z)
        return np.log((q@q)/(Lam**2*a**2))/(v@v)*F(z)
    # split radial range around the scale Lam*a where the log changes behaviour
    def inner(t):
        r=np.exp(t)
        g=lambda ph: f(np.array([r*np.cos(ph), r*np.sin(ph)]))
        val,_=integrate.quad(g,0,2*np.pi,limit=600,epsabs=1e-13,epsrel=1e-11,points=None)
        return val*r*r
    val=0
    edges=[-45,np.log(Lam*a)-3,np.log(Lam*a),np.log(Lam*a)+3,np.log(a),14]
    for lo,hi in zip(edges[:-1],edges[1:]):
        vv,_=integrate.quad(inner,lo,hi,limit=600,epsabs=1e-11,epsrel=1e-10); val+=vv
    L=np.log(1/Lam)
    pred=2*np.pi*L*L*Fw+2*L*I1+I2+np.pi**3/6*Fw
    print(f"Lam={Lam:8.1e}  exact={val:+.8f}  expansion={pred:+.8f}  diff={val-pred:+.2e}")
