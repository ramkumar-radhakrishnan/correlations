"""Group II Row IV: plus prescription in xi, and the split into evolution / DGLAP / finite."""
import mpmath as mp
mp.mp.dps=25

def run(K,Lam,Vee,kap,dperp,XX,RR,epsinv,Lr,EXACT=False):
    """XX=(X.X') , RR=(X.rhat)(X'.rhat) , epsinv=1/epsbar , Lr=log(4 pi^2 mu^2 r^2)"""
    M  = epsinv + 1 - Lr                       # the bracket of (0.20)
    x0 = K/Vee ; x1 = K/(K+Lam) ; lam = 1-x1   # lam = Lambda/(k+ + Lambda)
    C  = lambda x: x*(1-x) + x/(1-x) + (1-x)/x
    ph = lambda x: mp.e**(-1j*kap/x)
    # --- the exact integrand of (0.24), ascending limits ---
    I  = lambda x: ph(x)/x**2*( C(x)*M*XX - dperp*x*(1-x)*RR - (x/(1-x)+(1-x)/x)*XX )
    exact = mp.quad(I,[x0,0.5,x1])
    # --- the three sectors ---
    h  = lambda x: ph(x)/x                     # companion of 1/(1-xi)
    Lev= mp.log((1-x0)/lam)                    # exact log; -> log(k+/Lambda)
    evo = mp.e**(-1j*kap)*(epsinv - Lr)*XX*Lev
    plus= lambda x: (h(x)-h(1))/(1-x)          # == ph/x^2 * xi/[1-xi]_+
    top = x1 if EXACT else 1
    reg = mp.quad(lambda x: ph(x)/x**2*(x*(1-x)+(1-x)/x),[x0,0.5,top]) + mp.quad(plus,[x0,0.5,top])
    dgl = epsinv*XX*reg
    fin = -Lr*XX*reg + mp.quad(lambda x: ph(x)/x**2*x*(1-x)*(XX-dperp*RR),[x0,0.5,top])
    return exact, evo+dgl+fin, evo, dgl, fin

print("="*78); print("EVOLUTION + DGLAP + FINITE  vs  the exact (0.24)"); print("="*78)
for args in [(1.0,1e-4,50.0, 0.9,2.0, 1.7,-0.6, 12.0, 2.3),
             (2.5,5e-5,120.0,-1.4,2.0,-0.8, 1.1, -7.0,-1.1),
             (0.7,1e-3,30.0,  2.2,2.0, 2.4, 0.9,  5.0, 0.4)]:
    ex,su,ev,dg,fi = run(*args,EXACT=True)
    _,su2,_,_,_    = run(*args)
    print("   k+=%.1f kappa=%+.1f"%(args[0],args[3]))
    print("      exact        %s"%mp.nstr(ex,12))
    print("      evo+dglap+fin%s   diff %s"%(mp.nstr(su,12),mp.nstr(abs(ex-su),3)))
    print("      evolution    %s"%mp.nstr(ev,10))
    print("      dglap        %s"%mp.nstr(dg,10))
    print("      finite       %s"%mp.nstr(fi,10))
print()
print("   (the small residual is log[(1-xi0)/lam] vs log(k+/Lambda), i.e. the")
print("    1/V and Lambda/k+ corrections that are dropped when Lambda << k+ << V)")
