"""Row 14: the six p+ masters.  Phase e^{-i kappa p+} present, so Ei-type."""
import numpy as np, mpmath as mp
from scipy.integrate import quad
mp.mp.dps=25
K=1.0
Ei=lambda u: complex(mp.ei(u))
def prim(p,kap,A,be):
    """primitives of h(p) e^{-i kap p}/(p A + k+ B) = h e^{-i kap p}/(A(p+be))"""
    E =lambda a: np.exp(1j*kap*a)*Ei(-1j*kap*(p+a))
    out={}
    out["h=1"]        = E(be)/A
    out["h=1/p"]      = (Ei(-1j*kap*p) - E(be))/(A*be)
    out["h=p"]        = (np.exp(-1j*kap*p)/(-1j*kap) - be*E(be))/A
    out["h=1/S"]      = (E(K) - E(be))/(A*(be-K))
    out["h=p/S"]      = (-K*E(K) + be*E(be))/(A*(be-K))
    a2=-K/(be-K); a1=be/(be-K)**2; b=-be/(be-K)**2
    I2 = -np.exp(-1j*kap*p)/(p+K) - 1j*kap*E(K)      # int e^{-i kap p}/(p+K)^2
    out["h=p/S^2"]    = (a2*I2 + a1*E(K) + b*E(be))/A
    return out
hf={"h=1":lambda p:1.0+0*p, "h=1/p":lambda p:1/p, "h=p":lambda p:p,
    "h=1/S":lambda p:1/(p+K), "h=p/S":lambda p:p/(p+K), "h=p/S^2":lambda p:p/(p+K)**2}
A,B,Lam,V,kap = 0.62, 1.45, 1e-7, 40.0, 0.83
be=K*B/A; Pv=V-K
print("="*76); print("THE SIX p+ MASTERS,  M[h] = int_Lam^{V-k+} dp h(p) e^{-i kap p}/(pA + k+B)")
print("="*76)
print("   A=%.2f  B=%.2f  beta=k+B/A=%.4f  kappa=%.2f  Lambda=%.0e  V=%.0f"%(A,B,be,kap,Lam,V))
print("   NOTE beta > 0, so there is NO pole on the contour: the integrals are ordinary Ei's.\n")
print("   %-10s %30s %30s %10s"%("h","quadrature","closed form","diff"))
F0,F1=prim(Lam,kap,A,be),prim(Pv,kap,A,be)
for nm in hf:
    f=hf[nm]
    re=quad(lambda u: (f(np.exp(u))*np.exp(-1j*kap*np.exp(u))/(np.exp(u)*A+K*B)*np.exp(u)).real,
            np.log(Lam),np.log(Pv),limit=900)[0]
    im=quad(lambda u: (f(np.exp(u))*np.exp(-1j*kap*np.exp(u))/(np.exp(u)*A+K*B)*np.exp(u)).imag,
            np.log(Lam),np.log(Pv),limit=900)[0]
    q=re+1j*im; c=F1[nm]-F0[nm]
    print("   %-10s %30s %30s %10.1e"%(nm,"%.8f%+.8fj"%(q.real,q.imag),
                                       "%.8f%+.8fj"%(c.real,c.imag),abs(q-c)))
print()
print("   only  h = 1/p  is singular as Lambda -> 0 :  Ei(-i kappa Lambda) -> gamma_E + log(|kappa|Lambda)")
print("   so the RAPIDITY LOG is  -(1/(A beta)) log(1/Lambda) = -(1/(k+ B)) log(1/Lambda),")
print("   i.e. residue  -1/[k+ (x-w)^2]  on the structure  -P^{k'm} .")
print("   all five others are finite at Lambda = 0 and V-stable except h=p (grows like V).")
