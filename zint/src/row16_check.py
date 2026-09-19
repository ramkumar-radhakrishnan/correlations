"""Row 16: the two-piece B_3 row -- check the quoted p+ integration."""
import numpy as np, mpmath as mp
mp.mp.dps=30
K=1.0

def direct(A,B,Lam,V,which):
    """PV quadrature of the three p+ integrals, in units of the common i g^4 f/(16 pi^5)."""
    be=K*B/A
    def pv(f,lo,hi,d=1e-9):
        if lo<be<hi: return mp.quad(f,[lo,be-d])+mp.quad(f,[be+d,hi])
        return mp.quad(f,[lo,hi])
    if which==1:   # piece I on [k+Lam, V]
        f=lambda p:(p+K)/((K-p)**2*(p*A-K*B)); return pv(f,K+Lam,V)
    if which==2:   # piece I on [Lam, k-Lam]
        f=lambda p:(p+K)/((K-p)**2*(p*A-K*B)); return pv(f,Lam,K-Lam)
    if which==3:   # piece III on [Lam, V], and the quoted prefactor carries an extra -1
        f=lambda p:(K*B+p*A)/(p*(p*A-K*B))/K; return -pv(f,Lam,V)

def quoted(A,B,Lam,V,which):
    """the quoted terms, same units."""
    if which==1:
        t1 = 2/(A-B)*(1/Lam + 1/(K-V))
        t2 = -1/K*(A+B)/(A-B)**2*(np.log(V/Lam)+np.log(1-K/V))
        t3 =  1/K*(A+B)/(A-B)**2*np.log(abs((V*A-K*B)/((Lam+K)*A-K*B)))
        return t1+t2+t3
    if which==2:
        t1 = 2/(A-B)*(1/Lam - 1/(K-Lam))
        t2 = -1/K*(A+B)/(A-B)**2*(np.log(V/(K-Lam))-np.log(V/Lam))
        t3 =  1/K*(A+B)/(A-B)**2*np.log(abs(((K-Lam)*A-K*B)/(Lam*A-K*B)))
        return t1+t2+t3
    if which==3:
        return 1/K*np.log(V/Lam) - 2/K*np.log(abs((K*B-V*A)/(K*B-Lam*A)))

print("="*78); print("THE THREE p+ INTEGRALS,  A=(x-z)^2, B=(x-w)^2, beta = k+B/A"); print("="*78)
print("   units: the common  i g^4 f^{abd}/(16 pi^5)  (= the quoted i g^4 f/(8 pi^4) times 1/2pi)")
print()
print("  %-28s %6s %18s %18s %9s"%("(A,B,Lambda,V)","piece","PV quadrature","quoted","rel.diff"))
for A,B,Lam,V in [(0.62,1.45,1e-5,60.),(1.80,0.55,1e-6,120.),(0.90,0.90*2.3,1e-5,40.),
                  (2.40,1.10,1e-6,300.)]:
    be=K*B/A
    for wh in (1,2,3):
        q=float(direct(A,B,Lam,V,wh)); c=quoted(A,B,Lam,V,wh)
        lab="(%.2f,%.2f,%.0e,%.0f)"%(A,B,Lam,V) if wh==1 else ""
        print("  %-28s %6d %18.8f %18.8f %9.1e"%(lab,wh,q,c,abs(q-c)/max(abs(q),1e-30)))
    print("  %-28s beta = k+B/A = %.4f  -> inside [Lam,k-Lam]? %s   inside [k+Lam,V]? %s"
          %("",be,be<K-Lam,be>K+Lam))
    print()
