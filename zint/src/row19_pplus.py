"""Row 19: check the six p+-integrated parts 2a..2f against the exact masters."""
import numpy as np
from scipy.integrate import quad
L=np.log; s=lambda v:v@v; dot=lambda a,b:a@b

def geom(x,xp,z,w,wp):
    return s(x-z),s(x-w),s(z-w),s(xp-z),s(xp-wp),s(z-wp)

def exact_M(A,B,Ap,Bp,k,Lam,V):
    """M_ab = (1/(A A')) int_Lam^{V-k} dp h_ab /[(p+al)(p+al')], by quadrature."""
    al,be=k*B/A,k*Bp/Ap
    h={'23':lambda p:1.0+0*p,'22':lambda p:p/k,'33':lambda p:k/p,
       '13':lambda p:k/(p+k),'12':lambda p:p/(p+k),'11':lambda p:p*k/(p+k)**2}
    out={}
    for nm,fn in h.items():
        out[nm]=quad(lambda u: fn(np.exp(u))*np.exp(u)/((np.exp(u)+al)*(np.exp(u)+be)),
                     L(Lam),L(V-k),limit=900)[0]/(A*Ap)
    return out

def their_coef(A,B,C,Ap,Bp,Cp,k,Lam,V,dperp=2.0):
    """the factor each of your parts puts in front of its transverse bracket."""
    D  = A*Bp-Ap*B                     # (x'-w')^2(x-z)^2 - (x-w)^2(x'-z)^2
    LA = L(((V-k)*A+k*B)/(Lam*A+k*B))  # the unprimed log
    LAp= L(((V-k)*Ap+k*Bp)/(Lam*Ap+k*Bp))
    out={}
    # ---- 2a : the S^2/(p+k+) group,  h = 1
    out['23']=(1/k)*(1/D)*( L(((V-k)*A+k*B)/((V-k)*Ap+k*Bp)) - L((Lam*A+k*B)/(Lam*Ap+k*Bp)) )
    # ---- 2b : the S^2/p+^2 group,  h = k+/p+
    out['33']=( (1/k)*L((V-k)/Lam)/(Bp*A*B-B*B*Ap)*(A-B*Ap/Bp)
               -(1/k)*A/(Bp*A*B-B*B*Ap)*LA
               +(1/k)*Ap/(Bp*Bp*A-B*Ap*Bp)*LAp )
    # ---- 2c : the S^2/k+^2 group,  h = p+/k+
    out['22']=( (1/k)*(Bp/Ap)/D*LAp - (1/k)*(B/A)/D*LA )
    # ---- 2d : the S/p+ group,  h = k+/S
    out['13']=( -(1/k)*L(V/(k+Lam))/((A-B)*(Bp-Ap))
               +(1/k)*A/((A-B)*D)*LA
               +(1/k)*Ap/((Bp-Ap)*D)*LAp )
    # ---- 2e : the S/k+ group,  h = p+/S
    out['12']=(  (1/k)*L(V/(Lam+k))/((A-B)*(Bp-Ap))
               -(1/k)*B/((A-B)*D)*LA
               -(1/k)*Bp/((Bp-Ap)*D)*LAp )
    # ---- 2f : the d_perp group,  h = p+k+/S^2   (you pulled d_perp/4 and 1/(C C') out)
    f2f=( (dperp/4)*(1/k)*L(V/(Lam+k))*(Bp*B-A*Ap)/(C*Cp*(A-B)*(Ap-Bp))
         -(dperp/4)*(1/k)*(A*B)/(C*Cp*(A-B))*(Ap-Bp)/D*LA
         +(dperp/4)*(1/V-1/(k+Lam))/(C*Cp)
         +(dperp/4)*(1/k)*(Bp*Ap)/(C*Cp)*(A-B)/((Ap-Bp)*D)*LAp )
    out['11_full']=f2f
    return out

rng=np.random.default_rng(77); k=1.0
names={'23':"2a  S^2/p+k+  (h=1)",'33':"2b  S^2/p+^2  (h=k+/p+)",'22':"2c  S^2/k+^2  (h=p+/k+)",
       '13':"2d  S/p+      (h=k+/S)",'12':"2e  S/k+      (h=p+/S)",'11':"2f  d_perp    (h=p+k+/S^2)"}
print("="*80); print("YOUR COEFFICIENTS vs THE EXACT p+ MASTERS"); print("="*80)
for trial in range(4):
    x,xp,z,w,wp=[rng.normal(size=2)*1.3 for _ in range(5)]
    A,B,C,Ap,Bp,Cp=geom(x,xp,z,w,wp); Lam,V=10.0**rng.uniform(-8,-4),10.0**rng.uniform(4,8)
    M=exact_M(A,B,Ap,Bp,k,Lam,V); T=their_coef(A,B,C,Ap,Bp,Cp,k,Lam,V)
    t1t1=2.0*(A-B)*(Ap-Bp)/(4*C*Cp)          # tau_1 . tau_1'  with d_perp = 2
    print("\n trial %d   Lam=%.2e V=%.2e   A=%.4f B=%.4f A'=%.4f B'=%.4f"%(trial+1,Lam,V,A,B,Ap,Bp))
    print("   %-26s %18s %18s %10s"%("part","yours","exact","rel.diff"))
    for nm in ('23','33','22','13','12'):
        y,e=T[nm],M[nm]; print("   %-26s %18.10f %18.10f %10.1e"%(names[nm],y,e,abs(y-e)/max(abs(e),1e-30)))
    y,e=T['11_full'],M['11']*t1t1
    print("   %-26s %18.10f %18.10f %10.1e"%(names['11'],y,e,abs(y-e)/max(abs(e),1e-30)))
    print("   %-26s (2f is compared INCLUDING its transverse factor, since you folded it in)"%"")

print(); print("="*80); print("REMOVABLE SINGULARITIES: EACH PART BLOWS UP, THE SUM DOES NOT"); print("="*80)
def total(x,xp,z,w,wp,k,Lam,V):
    A,B,C,Ap,Bp,Cp=geom(x,xp,z,w,wp); T=their_coef(A,B,C,Ap,Bp,Cp,k,Lam,V)
    def taus(X,Z,W):
        AA,BB,CC=s(X-Z),s(X-W),s(Z-W)
        return (np.eye(2)*(AA-BB)/(2*CC), np.outer(X-Z,(Z-W)/CC-(X-W)/(2*BB)),
                np.outer((Z-W)/CC+(X-Z)/(2*AA),X-W))
    t1,t2,t3=taus(x,z,w); u1,u2,u3=taus(xp,z,wp); c=lambda a,b:np.einsum('ji,ji->',a,b)
    parts={'23':c(t2,u3)+c(t3,u2),'33':c(t3,u3),'22':c(t2,u2),'13':c(t1,u3)+c(t3,u1),'12':c(t1,u2)+c(t2,u1)}
    tot=T['11_full']+sum(T[n]*parts[n] for n in parts)
    return tot,{n:T[n]*parts[n] for n in parts},T['11_full']
x,xp,z,w,wp=[rng.normal(size=2)*1.3 for _ in range(5)]; Lam,V=1e-6,1e6
print("  (a) approaching  (x-z)^2 = (x-w)^2   [alpha -> k+], where 2d,2e,2f each have 1/(A-B):")
for eps in (1e-1,1e-3,1e-5,1e-7):
    A0=s(x-z); ww=x+(w-x)*np.sqrt((A0*(1-eps))/s(x-w))
    tot,pp,p11=total(x,xp,z,ww,wp,k,Lam,V)
    print("     A/B-1=%8.1e   2d=%12.3e  2e=%12.3e  2f=%12.3e   SUM=%14.9f"%(eps,pp['13'],pp['12'],p11,tot))
print("  (b) approaching  A B' = A' B   [alpha -> alpha'], where 2a..2f have 1/(AB'-A'B):")
for eps in (1e-1,1e-3,1e-5,1e-7):
    A,B,Ap,Bp=s(x-z),s(x-w),s(xp-z),s(xp-wp)
    wwp=xp+(wp-xp)*np.sqrt((Ap*B/A*(1+eps))/Bp)
    tot,pp,p11=total(x,xp,z,w,wwp,k,Lam,V)
    print("     AB'/A'B-1=%8.1e 2a=%12.3e  2c=%12.3e  2f=%12.3e   SUM=%14.9f"%(eps,pp['23'],pp['22'],p11,tot))
