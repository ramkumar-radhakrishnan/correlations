"""Group II Row IV: the p+ integration with the plus prescription."""
import numpy as np, mpmath as mp
mp.mp.dps=25

print("="*78); print("1.  IS THE 1/eps OF  G^{mm'}  UV OR IR ?"); print("="*78)
print("   trace:  delta_mm' G^mm' = int d^2z (y-z).(y'-z)/[(y-z)^2 (y'-z)^2]")
print("   small |y-z| : integrand ~ 1/|y-z| , integrable in 2d   -> NO UV")
print("   large |z|   : integrand ~ 1/z^2                        -> LOG DIVERGENT (infrared)")
def trace_cut(r,R,n=4000,m=1024):
    """int_{|u|<R} d^2u  u.(u+r)/[u^2 (u+r)^2] , y at 0 , y'-y = r along x ; log-graded in |u|"""
    lo,hi=np.log(1e-8*r),np.log(R)
    xg,wg=np.polynomial.legendre.leggauss(n)
    lr=0.5*(hi-lo)*xg+0.5*(hi+lo); wr=0.5*(hi-lo)*wg
    rho=np.exp(lr)
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    ux=np.outer(rho,np.cos(ph)); uy=np.outer(rho,np.sin(ph))
    vx=ux+r; vy=uy
    num=ux*vx+uy*vy; den=(ux**2+uy**2)*(vx**2+vy**2)
    w=(wr*rho*rho)[:,None]*(2*np.pi/m)          # r dr dphi = r^2 dlog r dphi
    return float((num/den*w).sum())
r=1.0
print("   %-10s %14s %16s"%("R","numeric","pi*log(R^2/r^2)"))
for R in (20.,60.,180.,540.):
    print("   %-10.0f %14.6f %16.6f"%(R,trace_cut(r,R),np.pi*np.log(R**2/r**2)))
print("   -> grows like pi log R^2 : the pole is INFRARED (|z| -> infinity), not ultraviolet.")
print("      In d = 2-2eps it is regulated by eps > 0 , so 1/eps = 1/eps_IR .")

print(); print("="*78); print("2.  THE PLUS PRESCRIPTION AT  p+ -> 0"); print("="*78)
print("   the only p+ -> 0 singularity is  1/p+  from the last bracket term;")
print("   its companion is the phase  F(p+) = exp(-i kappa p+/k+) ,  F(0) = 1  FINITE.")
print("   => int_Lam^{V-k+} dp+ F(p+)/p+ = F(0) log[(V-k+)/Lam] + int_0^{V-k+} dp+ [F(p+)-F(0)]/p+")
def lhs(kap,lam,T): return mp.quad(lambda t: mp.e**(-1j*kap*t)/t,[lam,T])
def rhs(kap,lam,T): return mp.log(T/lam)+mp.quad(lambda t:(mp.e**(-1j*kap*t)-1)/t,[0,T])
def closed(kap,T):
    z=1j*kap; return -mp.euler-mp.log(z*T)-mp.e1(z*T)
print("   %-22s %24s %24s"%("(kappa,lam,T)","direct","plus-prescription form"))
for kap,lam,T in ((0.7,1e-5,30.),(2.0,1e-6,15.),(0.2,1e-7,80.)):
    print("   k=%.1f lam=%.0e T=%-6.0f %24s %24s"%(kap,lam,T,mp.nstr(lhs(kap,lam,T),10),mp.nstr(rhs(kap,lam,T),10)))
print()
print("   and the subtracted integral in closed form:")
print("      int_0^T dt [e^{-i kappa t} - 1]/t  =  -gammaE - log(i kappa T) - E1(i kappa T)")
for kap,T in ((0.7,30.),(2.0,15.),(0.2,80.)):
    num=mp.quad(lambda t:(mp.e**(-1j*kap*t)-1)/t,[0,T])
    print("      kappa=%.1f T=%-5.0f  numeric %s   closed %s"%(kap,T,mp.nstr(num,10),mp.nstr(closed(kap,T),10)))

print(); print("="*78); print("3.  THE TWO MASTERS"); print("="*78)
def J1c(kap):
    z=1j*kap; return (1+z)*mp.e1(z)-mp.e**(-z)
def J1n(kap):
    f=lambda s:(s-1)/s**2*mp.e**(-1j*kap*s); T0=1+8*2*mp.pi/kap
    return mp.quad(f,[1,T0])+mp.quadosc(f,[T0,mp.inf],period=2*mp.pi/kap)
print("   J1 = (1+i kappa) E1(i kappa) - e^{-i kappa}")
for kap in (0.4,1.2,3.0):
    print("      kappa=%.1f  numeric %s   closed %s"%(kap,mp.nstr(J1n(kap),10),mp.nstr(J1c(kap),10)))
print()
print("   J2 = e^{-i kappa} { [e^{-i kappa T}(1+i kappa T)-1]/kappa^2 + log(T/lam) - gammaE - log(i kappa T) - E1(i kappa T) }")
def J2n(kap,lam,T): return mp.quad(lambda t:(t+1/t)*mp.e**(-1j*kap*(1+t)),[lam,T])
def J2c(kap,lam,T):
    z=1j*kap
    lin=(mp.e**(-z*T)*(1+z*T)-1)/kap**2
    return mp.e**(-z)*(lin+mp.log(T/lam)-mp.euler-mp.log(z*T)-mp.e1(z*T))
for kap,lam,T in ((0.7,1e-5,30.),(2.0,1e-6,15.)):
    print("      kappa=%.1f lam=%.0e T=%-4.0f numeric %s   closed %s"%(kap,lam,T,mp.nstr(J2n(kap,lam,T),10),mp.nstr(J2c(kap,lam,T),10)))
print()
print("   log(T/lam) - log(i kappa T) = log(1/lam) - log(i kappa) = log(k+/Lambda) - log(i kappa):")
print("   the cutoff T = V/k+ - 1 DROPS OUT of the logarithm, as it must.")
print("   J2 -> e^{-i kappa}[ log(k+/Lambda) - gammaE - log(i kappa) - 1/kappa^2 ] + (oscillating in kappa T)")
