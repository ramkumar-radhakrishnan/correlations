"""Row 25: the z-integrated integrand as a function of p+ -- single or double log?"""
import numpy as np
from row25_lib import val, KT, KP, s
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
kabs=np.sqrt(KT@KT)
def zint(P,Rfac=60.0,nr=1400,nth=6000):
    R=max(1e3,Rfac*KP/(P*kabs))
    lo,hi=np.log(1e-3),np.log(R); lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo
    tot=0j
    for r0 in np.exp(lr):
        tot+=val(ring(r0,nth),P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot
print("="*78); print("THE z-INTEGRATED INTEGRAND vs p+   (phase on, R = 60 k+/(p+|k|))"); print("="*78)
print("   %-9s %16s %16s %18s"%("p+","I(p+)","I x p+","I x p+ / log(k+/p+)"))
rows=[]
for P in (1.0,0.3,0.1,0.03,0.01,0.003,0.001):
    I=zint(P).real; rows.append((P,I))
    print("   %-9.0e %16.7f %16.7f %18.7f"%(P,I,I*P,I*P/np.log(KP/P)),flush=True)
print()
import numpy as _np
P=_np.array([r[0] for r in rows]); I=_np.array([r[1] for r in rows])
sm=P<0.2
A=_np.vstack([_np.ones(sm.sum()),_np.log(KP/P[sm])]).T
co,*_=_np.linalg.lstsq(A,(I*P)[sm],rcond=None)
print("   fit over p+ < 0.2 :  I(p+) x p+  =  %.6f  +  %.6f log(k+/p+)"%(co[0],co[1]))
print("   a non-zero log coefficient would mean I ~ (1/p+) log(1/p+)  -> a DOUBLE log after")
print("   the p+ integration;  a coefficient consistent with zero means a single rapidity log.")
print()
print("   ratio of the log term to the constant at p+ = 1e-3 :  %.4f"%(co[1]*np.log(KP/1e-3)/abs(co[0])))
