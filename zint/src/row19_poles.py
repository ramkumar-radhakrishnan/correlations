"""Row 19: are the 1/(A-B) and 1/(AB'-A'B) poles really removable?  High precision."""
import mpmath as mp
mp.mp.dps=60
L=mp.log
def coef(A,B,C,Ap,Bp,Cp,k,Lam,V,dperp=2):
    A,B,C,Ap,Bp,Cp,k,Lam,V=[mp.mpf(v) for v in (A,B,C,Ap,Bp,Cp,k,Lam,V)]
    D=A*Bp-Ap*B; LA=L(((V-k)*A+k*B)/(Lam*A+k*B)); LAp=L(((V-k)*Ap+k*Bp)/(Lam*Ap+k*Bp))
    o={}
    o['2a']=(1/k)/D*( L(((V-k)*A+k*B)/((V-k)*Ap+k*Bp)) - L((Lam*A+k*B)/(Lam*Ap+k*Bp)) )
    o['2b']=( (1/k)*L((V-k)/Lam)/(Bp*A*B-B*B*Ap)*(A-B*Ap/Bp) - (1/k)*A/(Bp*A*B-B*B*Ap)*LA
             +(1/k)*Ap/(Bp*Bp*A-B*Ap*Bp)*LAp )
    o['2c']=( (1/k)*(Bp/Ap)/D*LAp - (1/k)*(B/A)/D*LA )
    o['2d']=( -(1/k)*L(V/(k+Lam))/((A-B)*(Bp-Ap)) + (1/k)*A/((A-B)*D)*LA + (1/k)*Ap/((Bp-Ap)*D)*LAp )
    o['2e']=(  (1/k)*L(V/(Lam+k))/((A-B)*(Bp-Ap)) - (1/k)*B/((A-B)*D)*LA - (1/k)*Bp/((Bp-Ap)*D)*LAp )
    o['2f']=( (dperp/mp.mpf(4))*(1/k)*L(V/(Lam+k))*(Bp*B-A*Ap)/(C*Cp*(A-B)*(Ap-Bp))
             -(dperp/mp.mpf(4))*(1/k)*(A*B)/(C*Cp*(A-B))*(Ap-Bp)/D*LA
             +(dperp/mp.mpf(4))*(1/V-1/(k+Lam))/(C*Cp)
             +(dperp/mp.mpf(4))*(1/k)*(Bp*Ap)/(C*Cp)*(A-B)/((Ap-Bp)*D)*LAp )
    return o
k,Lam,V=1,mp.mpf('1e-6'),mp.mpf('1e6')
A,Ap,C,Cp,Bp=mp.mpf('2.3'),mp.mpf('3.1'),mp.mpf('1.7'),mp.mpf('0.9'),mp.mpf('1.4')
print("="*76); print("(a)  B -> A   (alpha -> k+) :  2d, 2e, 2f each carry 1/(A-B)"); print("="*76)
print("  %-10s %20s %20s %20s"%("A/B-1","2d","2e","2f"))
for e in ('1e-2','1e-6','1e-10','1e-20','1e-40','0'):
    B=A*(1-mp.mpf(e)) if e!='0' else A
    if e=='0': B=A*(1-mp.mpf('1e-50'))
    o=coef(A,B,C,Ap,Bp,Cp,k,Lam,V)
    print("  %-10s %20s %20s %20s"%(e,mp.nstr(o['2d'],12),mp.nstr(o['2e'],12),mp.nstr(o['2f'],12)))
print("  -> all three converge: the poles ARE removable, cancelling inside each part.")
print()
print("="*76); print("(b)  A B' -> A' B   (alpha -> alpha') :  every part carries 1/(AB'-A'B)"); print("="*76)
B=mp.mpf('0.8')
print("  %-10s %16s %16s %16s %16s"%("AB'/A'B-1","2a","2b","2c","2f"))
for e in ('1e-2','1e-6','1e-10','1e-20','1e-40'):
    Bp2=Ap*B/A*(1+mp.mpf(e))
    o=coef(A,B,C,Ap,Bp2,Cp,k,Lam,V)
    print("  %-10s %16s %16s %16s %16s"%(e,mp.nstr(o['2a'],10),mp.nstr(o['2b'],10),
                                          mp.nstr(o['2c'],10),mp.nstr(o['2f'],10)))
print("  -> all finite.  (2d, 2e likewise.)")
print()
print("="*76); print("(c)  WHERE DOUBLE PRECISION BREAKS"); print("="*76)
import numpy as np
from row19_lib import their_coef
print("  %-12s %18s %18s %10s"%("AB'/A'B-1","2b  float64","2b  60-digit","rel.err"))
for e in (1e-2,1e-6,1e-8,1e-10,1e-12):
    Bp2f=float(Ap)*float(B)/float(A)*(1+e)
    f64=their_coef(float(A),float(B),float(C),float(Ap),Bp2f,float(Cp),1.0,1e-6,1e6)['33']
    hi=coef(A,B,C,Ap,Ap*B/A*(1+mp.mpf(repr(e))),Cp,k,Lam,V)['2b']
    print("  %-12.0e %18.9f %18s %10.1e"%(e,f64,mp.nstr(hi,10),abs(f64-float(hi))/abs(float(hi))))
print("  => below ~1e-6 in that ratio, float64 loses the cancellation entirely.")
print("     Expand analytically near A=B, A'=B', AB'=A'B before integrating numerically.")
