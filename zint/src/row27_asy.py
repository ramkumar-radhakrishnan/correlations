"""Row 27: the analytic large-|z| form, and the Bessel masters for the z integral."""
import numpy as np, mpmath as mp
from row25_lib import val, KT, KP, s
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
X=(xp-yp)/s(xp-yp); B=s(x-w); kabs=np.sqrt(KT@KT)
Kc=((x-w)@X)/(2*B)
print("="*78); print("1.  THE ANALYTIC LARGE-|z| FORM"); print("="*78)
print("   keeping only the O(|z|) parts of the bracket (the others give |z|^-3 or less):")
print("     -p+/k+^2 group : V_w^m (x-z)^k'   ->  -(x-w)^m/(2B) x (-z^k')")
print("     -1/k+   group : (x-z)^m V_w^k'   ->  (-z^m) x (-(x-w)^k'/(2B))")
print("     S/k+ (delta) group : -(x-z).(x-w)/2B  ->  +z.(x-w)/(2B)")
print("   contracting with N^m = -z^m/|z|^2 and X^k', and dividing by D -> p+|z|^2 :")
print("     <I>_zhat |z|^2  =  (K/p+) [ 1/k+ + p+^2/(2 k+^2 S) ] ,   K = (x-w).X / (2B) = %.9f"%Kc)
print()
print("   %-8s %18s %18s %10s"%("p+","measured","analytic","ratio"))
for P in (0.05,0.15,0.37,1.0,3.0):
    S=P+KP
    meas=val(ring(3e4,40000),P,phase=False).mean().real*9e8
    ana=(Kc/P)*(1/KP+P**2/(2*KP**2*S))
    print("   %-8.2f %18.9f %18.9f %10.5f"%(P,meas,ana,meas/ana))
print("   -> the analytic form is confirmed.  At small p+ it is simply  K/(p+ k+).")

print(); print("="*78); print("2.  THE TWO BESSEL MASTERS"); print("="*78)
mp.mp.dps=20
print("   <e^{iq.z}>_zhat = J_0(q|z|) , so")
print("     M_0 = int_{|z|>rho0} d^2z e^{iq.z}/|z|^2 = 2 pi int_{q rho0}^inf dt J_0(t)/t")
print("         = 2 pi [ log(2/(q rho0)) - gamma_E ]          (for q rho0 << 1)")
for e in ('1e-2','1e-4','1e-6'):
    ee=mp.mpf(e); num=mp.quad(lambda t: mp.besselj(0,t)/t,[ee,1,10,50,mp.inf])
    cl=mp.log(2/ee)-mp.euler
    print("       eps=%-6s  numeric %s   log(2/eps)-gamma_E %s   diff %.1e"%(e,mp.nstr(num,12),mp.nstr(cl,12),abs(num-cl)))
print()
print("   tensor version, using  <zhat^m zhat^k e^{i x zhat.qhat}> = (J_1(x)/x) delta^{mk}")
print("                                            + (J_0(x) - 2J_1(x)/x) qhat^m qhat^k :")
a2=mp.quad(lambda t:(mp.besselj(0,t)-2*mp.besselj(1,t)/t)/t,[0,1,10,50,mp.inf])
print("     the qhat qhat piece is FINITE at small t (integrand ~ -t/8) :")
print("        a2 = int_0^inf dt [J_0 - 2J_1/t]/t = %s"%mp.nstr(a2,12))
print("     and the delta piece is half the log, fixed by the trace  2 a1 + a2 = log(2/eps) - gamma_E :")
for e in ('1e-4','1e-6'):
    ee=mp.mpf(e); a1=mp.quad(lambda t: mp.besselj(1,t)/t**2,[ee,1,10,50,mp.inf])
    print("        eps=%-6s  a1 numeric %s   (log(2/eps)-gamma_E-a2)/2 = %s"
          %(e,mp.nstr(a1,12),mp.nstr((mp.log(2/ee)-mp.euler-a2)/2,12)))
print()
print("   => M_2^{mk'} = 2 pi [ a1 delta^{mk'} + a2 khat^m khat^k' ] ,  a1 = (1/2)[log(2/(q rho0)) - gamma_E - a2]")
print("      so the LOG sits in the delta part with half weight -- exactly the angular average 1/2.")

print(); print("="*78); print("3.  THE LOG COEFFICIENT, PREDICTED vs MEASURED"); print("="*78)
pred=2*np.pi*Kc/KP
print("   from the masters:  I(p+) p+  ->  2 pi K/k+ x log(.../p+) ,  2 pi K/k+ = %.6f"%pred)
def zint(P,Rfac=80.0,nr=1600,nth=7000):
    R=max(1e3,Rfac*KP/(P*kabs)); lo,hi=np.log(1e-3),np.log(R)
    lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; tot=0j
    for r0 in np.exp(lr): tot+=val(ring(r0,nth),P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot.real
Ps=[0.1,0.03,0.01,0.003,0.001,3e-4,1e-4]
I=[zint(P) for P in Ps]
print("   %-10s %16s %16s %16s"%("p+","I(p+) p+","local slope","predicted"))
for i,(P,Iv) in enumerate(zip(Ps,I)):
    if i==0: print("   %-10.0e %16.7f %16s %16s"%(P,Iv*P,"-","-"))
    else:
        sl=(I[i]*Ps[i]-I[i-1]*Ps[i-1])/np.log(Ps[i-1]/Ps[i])
        print("   %-10.0e %16.7f %16.7f %16.6f"%(P,Iv*P,sl,pred),flush=True)
print("   -> the LOCAL slope converges to the analytic 2 pi K/k+ .")
print("   (my earlier global least-squares fit, -0.0490, was contaminated by the larger-p+")
print("    points where the asymptotics has not set in; the correct coefficient is the above.)")
