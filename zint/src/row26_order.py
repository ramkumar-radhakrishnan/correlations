"""Row 26: what happens if you do the p+ integral BEFORE the z integral."""
import numpy as np, mpmath as mp
from row25_lib import val, KT, KP, s
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
kabs=np.sqrt(KT@KT)
print("="*78); print("1.  THE p+ MASTERS AT FIXED z ARE NOT LOGARITHMS -- THEY ARE Ei's"); print("="*78)
print("   at fixed z the p+ dependence is  (rational in p+) x exp(-i alpha p+) ,")
print("      alpha = k.(y'-z)/k+   comes from the phase e^{-i (p+/k+) k.(y'-z)} .")
print("   so every master is an exponential integral, e.g.")
print("      int dp+ e^{-i a p+}/p+            = -Ei(-i a p+)          (or E_1)")
print("      int dp+ e^{-i a p+}/(p+ A + k+ B) = (1/A) e^{i a k+B/A} Ei(-i a (p+ + k+B/A))")
mp.mp.dps=25
a=mp.mpf('3.7'); A=mp.mpf('0.9'); B=mp.mpf('1.4'); k=mp.mpf(1)
lo,hi=mp.mpf('0.05'),mp.mpf('2.0')
num=mp.quad(lambda p: mp.e**(-1j*a*p)/(p*A+k*B),[lo,hi])
cl =(1/A)*mp.e**(1j*a*k*B/A)*(mp.ei(-1j*a*(hi+k*B/A))-mp.ei(-1j*a*(lo+k*B/A)))
print("      numeric %s   closed %s   diff %.1e"%(mp.nstr(num,10),mp.nstr(cl,10),abs(num-cl)))
print("   -> doing p+ first turns every term into Ei's of a z-dependent argument, and you")
print("      still have to do the five transverse integrals afterwards.")

print(); print("="*78); print("2.  AND THE p+-INTEGRATED RESULT HAS A 1/|z|^2 TAIL"); print("="*78)
def J_ring(R,nth,Lam,P,npn=12000):
    """int_Lam^P dp+ of the integrand, averaged over a ring |z|=R"""
    lo,hi=np.log(Lam),np.log(P); lp=(np.arange(npn)+.5)*(hi-lo)/npn+lo
    acc=np.zeros(nth,dtype=complex); Z=ring(R,nth)
    for pv in np.exp(lp): acc+=val(Z,pv)*pv*((hi-lo)/npn)
    return acc.mean()
print("   <J(z)> x |z|^2 , with Lambda = 1e-6 , p+ up to 5 :")
print("   %-8s %18s %18s"%("|z|","<J> x |z|^2","ratio to previous"))
prev=None
for R in (30.,60.,120.,240.):
    v=(J_ring(R,600,1e-6,5.0)).real*R*R
    print("   %-8.0f %18.7f %18s"%(R,v,"-" if prev is None else "%.3f"%(v/prev)))
    prev=v
print("   -> does NOT fall: a flat (indeed slowly growing) |z|^-2 tail, so the z integral")
print("      done AFTER p+ is logarithmically divergent and needs a transverse cutoff R.")
print()
print("   and that cutoff is not free -- it is tied to Lambda:")
print("   %-10s %18s"%("Lambda","<J> x |z|^2 at |z|=120"))
for Lam in (1e-4,1e-6,1e-8):
    print("   %-10.0e %18.7f"%(Lam,(J_ring(120.,600,Lam,5.0)).real*120**2))
print("   the tail grows with log(1/Lambda): the transverse log and the rapidity log are")
print("   the SAME log seen from two sides.  Cutting z at R and p+ at Lambda independently")
print("   double counts unless you impose R ~ k+/(Lambda |k_perp|).")
