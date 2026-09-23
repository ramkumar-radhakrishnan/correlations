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

print(); print("="*78); print("WHY: THE TWO []_+ ARE DIFFERENT DISTRIBUTIONS"); print("="*78)
K,Lam,Vee,kap=1.0,1e-5,60.0,0.9; x0=K/Vee; T=Vee-K
F =lambda p: mp.e**(-1j*kap*(K+p)/K)
Ft=lambda x: mp.e**(-1j*kap/x)
a=mp.quad(lambda p:(F(p)-F(0))/p,[0,K,T])                        # [1/p+]_+ in p+
b=mp.quad(lambda x:(Ft(x)-Ft(1))/(x*(1-x)),[x0,0.5,1])           # the same, rewritten in xi
c=mp.quad(lambda x:(Ft(x)/x-Ft(1))/(1-x),[x0,0.5,1])             # xi/[1-xi]_+
print("   [1/p+]_+ in p+                :", mp.nstr(a,12))
print("   [1/p+]_+ rewritten in xi      :", mp.nstr(b,12), "  (weight 1/(xi xibar) on BOTH terms)")
print("   xi/[1-xi]_+                   :", mp.nstr(c,12))
print("   c - a                         :", mp.nstr(c-a,12))
print("   F(0) log(V/k+)                :", mp.nstr(Ft(1)*mp.log(Vee/K),12), "  <-- exactly the difference")
print("   and the two delta coefficients differ by log[(V-k+)/k+], the same up to O(k+/V).")
