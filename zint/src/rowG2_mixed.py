"""Group II Row IV: the xi-form with []_+ -- is the delta coefficient log[(V-k+)/Lam] or log(k+/Lam)?"""
import mpmath as mp
mp.mp.dps=25
def run(K,Lam,Vee,kap,dperp,XX,RR,epsinv,Lr):
    M=epsinv+1-Lr; x0=K/Vee; x1=K/(K+Lam)
    C =lambda x: x*(1-x)+x/(1-x)+(1-x)/x
    ph=lambda x: mp.e**(-1j*kap/x)
    # exact starting point, in xi, before any prescription
    exact=mp.quad(lambda x: ph(x)/x**2*( C(x)*M*XX - dperp*x*(1-x)*RR
                                         -(x/(1-x)+(1-x)/x)*XX ),[x0,0.5,x1])
    h=lambda x: ph(x)/x                       # the full companion of 1/(1-xi)
    plus=lambda x: (h(x)-h(1))/(1-x)          # == ph/x^2 * xi/[1-xi]_+
    reg =lambda x: ph(x)/x**2*( (x*(1-x)+(1-x)/x)*M*XX - dperp*x*(1-x)*RR - ((1-x)/x)*XX )
    body= mp.quad(reg,[x0,0.5,1]) + (M-1)*XX*mp.quad(plus,[x0,0.5,1])
    Lxi = h(1)*(M-1)*XX*mp.log(K/Lam)         # xi scheme, approximated
    Lpp = h(1)*(M-1)*XX*mp.log((Vee-K)/Lam)   # p+ scheme  (what is written)
    # the same xi decomposition kept EXACT: upper limit xi1, exact log
    bodyE= mp.quad(reg,[x0,0.5,x1]) + (M-1)*XX*mp.quad(plus,[x0,0.5,x1])
    LxiE = h(1)*(M-1)*XX*mp.log((1-x0)/(1-x1))
    return exact, body+Lxi, body+Lpp, bodyE+LxiE
print("="*78); print("THE xi FORM WITH []_+ : WHICH LOG ON THE DELTA TERM ?"); print("="*78)
for a in [(1.0,1e-5,60.0, 0.9,2.0, 1.7,-0.6, 12.0, 2.3),
          (2.5,2e-6,150.0,-1.4,2.0,-0.8, 1.1, -7.0,-1.1),
          (0.7,1e-5,40.0,  2.2,2.0, 2.4, 0.9,  5.0, 0.4)]:
    ex,sxi,spp,sex = run(*a)
    print("   k+=%.1f kappa=%+.1f"%(a[0],a[3]))
    print("      exact                         %s"%mp.nstr(ex,12))
    print("      []_+ in xi  + log(k+/Lam)     %s   diff %s"%(mp.nstr(sxi,12),mp.nstr(abs(ex-sxi),3)))
    print("      []_+ in xi  + log((V-k+)/Lam) %s   diff %s   <-- MIXED SCHEMES"%(mp.nstr(spp,12),mp.nstr(abs(ex-spp),3)))
    print("      same, exact limits and log    %s   diff %s   <-- EXACT"%(mp.nstr(sex,12),mp.nstr(abs(ex-sex),3)))
print()
print("   the []_+ subtraction is in xi, so the delta coefficient must be the xi one,")
print("   log(k+/Lambda).  Keeping log[(V-k+)/Lambda] double counts log[(V-k+)/k+].")
