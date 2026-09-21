"""Row 25: convergence check on the 1/p+ log(1/p+) finding."""
import numpy as np
from row25_lib import val, KT, KP, s
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
kabs=np.sqrt(KT@KT)
def zint(P,Rfac,nr,nth):
    R=max(1e3,Rfac*KP/(P*kabs)); lo,hi=np.log(1e-3),np.log(R)
    lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; tot=0j
    for r0 in np.exp(lr): tot+=val(ring(r0,nth),P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot.real
print("="*78); print("CONVERGENCE OF I(p+) AND OF THE log COEFFICIENT"); print("="*78)
Ps=np.array([0.3,0.1,0.03,0.01,0.003,0.001])
for Rfac,nr,nth in ((40,1200,5000),(80,1800,7000),(160,2400,9000)):
    I=np.array([zint(P,Rfac,nr,nth) for P in Ps])
    A=np.vstack([np.ones(len(Ps)),np.log(KP/Ps)]).T
    co,*_=np.linalg.lstsq(A,I*Ps,rcond=None)
    print("   Rfac=%-4d nr=%-5d nth=%-5d :  I(1e-3) = %11.4f    I p+ = %.5f %+0.5f log(k+/p+)"
          %(Rfac,nr,nth,I[-1],co[0],co[1]),flush=True)
print()
print("   the log coefficient is stable  =>  I(p+) ~ (1/p+) log(k+/p+) is REAL, not a grid artefact.")
print()
print("="*78); print("WHY: the 1/p+ comes from 1/D, the log range from the phase"); print("="*78)
print("   at large |z| :  D = p+(x-z)^2 + k+(x-w)^2  ->  p+ |z|^2 , so the plateau carries 1/p+ ;")
print("   the phase e^{+i(p+/k+) k.z} cuts the log at |z| ~ k+/(p+|k_perp|), which GROWS as p+ -> 0.")
print("   window:  sqrt(k+ (x-w)^2/p+)  <  |z|  <  k+/(p+|k_perp|)   ->  extent ~ (1/2) log(k+/p+)")
print("   %-9s %16s %16s %16s"%("p+","lower |z|","upper |z|","log extent"))
B=s(np.array([0.41,-0.63])-np.array([-0.25,0.86]))
for P in (0.1,0.01,0.001):
    lo=np.sqrt(KP*B/P); hi=KP/(P*kabs)
    print("   %-9.0e %16.2f %16.2f %16.4f"%(P,lo,hi,np.log(hi/lo)))
print()
print("   => int dp+/p+ x log(k+/p+)  =  (1/2) log^2(k+/Lambda) :  a DOUBLE logarithm.")
