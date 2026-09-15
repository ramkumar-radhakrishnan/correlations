"""Just part C's truncated integral: does dF/dlogR grow (=> -l_k logR + (1/2)log^2 R)?"""
import numpy as np
from scipy.special import exp1
from scipy.integrate import quad
def Ei_mi(u):
    u = np.asarray(u, float); return -exp1(1j*u) - 1j*np.pi*np.sign(u)
rng = np.random.default_rng(3)
zz, xx, xpp, yy, ww = rng.normal(size=(5,2))
kk = np.array([0.9,0.6]); d, e = xx-zz, yy-ww
lam, Xi = 1e-6, 37.0
th = np.linspace(0,2*np.pi,1024,endpoint=False)
u = np.stack([np.cos(th),np.sin(th)],-1)
def ring(rho):
    r = u*rho; r2=np.sum(r*r,-1); C=xpp-zz-r; C2=np.sum(C*C,-1)
    T2 = (r@d)/(r2*(d@d)) * (C@e)/(C2*(e@e))
    kap = r@kk
    return rho*np.mean(-(Ei_mi(kap*Xi)-Ei_mi(kap*lam))*T2)*2*np.pi
def cum(a,b,n=50):
    ed=np.geomspace(a,b,n)
    return sum(quad(lambda x: ring(x).real,p,q,limit=50)[0] for p,q in zip(ed[:-1],ed[1:]))
print("   truncated int_{|y'-z|<R}, real part, and its slope")
print(f"   l_k = log(1/lam) = {np.log(1/lam):.4f}")
prev=None
for R in (1e2,1e3,1e4,1e5,1e6):
    v=cum(1e-8,R)
    s = f"   dF/dlogR = {(v-prev)/np.log(10):+11.4f}" if prev is not None else ""
    pred = -np.log(1/lam) + np.log(R/np.sqrt(10))   # -l_k + log R, up to constants
    print(f"   F(R={R:.0e}) = {v:+14.5f}{s}" + (f"   [-l_k + log R ~ {pred:+.3f}] x 2pi x coef" if prev is not None else ""))
    prev=v
