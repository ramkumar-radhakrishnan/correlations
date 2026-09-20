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

