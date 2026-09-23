"""Group II Row IV: if the subtraction was done in p+ FIRST, what does the xi bracket mean?"""
import mpmath as mp
mp.mp.dps=25
def run(K,Lam,Vee,kap,dperp,XX,RR,epsinv,Lr):
    M=epsinv+1-Lr; x0=K/Vee; x1=K/(K+Lam); T=Vee-K
    ph=lambda x: mp.e**(-1j*kap/x); ph1=mp.e**(-1j*kap)
    # exact starting point (no prescription), in xi
    C=lambda x: x*(1-x)+x/(1-x)+(1-x)/x
    exact=mp.quad(lambda x: ph(x)/x**2*( C(x)*M*XX - dperp*x*(1-x)*RR
                                         -(x/(1-x)+(1-x)/x)*XX ),[x0,0.5,x1])
    reg=lambda x: ph(x)/x**2*( (x*(1-x)+(1-x)/x)*M*XX - dperp*x*(1-x)*RR - ((1-x)/x)*XX )
    # (a) subtraction done in p+ first : image is 1/(xi xibar) [ph(xi) - ph(1)] , delta = log[(V-k+)/Lam]
    subA=mp.quad(lambda x:(ph(x)-ph1)/(x*(1-x)),[x0,0.5,1])
    A = ph1*(M-1)*XX*mp.log(T/Lam) + mp.quad(reg,[x0,0.5,1]) + (M-1)*XX*subA
    # (b) the literal reading of  xi/[1-xi]_+  with delta = log[(V-k+)/Lam]  (mixed)
    subB=mp.quad(lambda x:(ph(x)/x-ph1)/(1-x),[x0,0.5,1])
    B = ph1*(M-1)*XX*mp.log(T/Lam) + mp.quad(reg,[x0,0.5,1]) + (M-1)*XX*subB
    # (c) the literal reading of  xi/[1-xi]_+  with delta = log(k+/Lam)   (consistent xi scheme)
    Cc= ph1*(M-1)*XX*mp.log(K/Lam) + mp.quad(reg,[x0,0.5,1]) + (M-1)*XX*subB
    return exact,A,B,Cc
print("="*80); print("WHICH READING OF THE xi BRACKET IS RIGHT?"); print("="*80)
for a in [(1.0,1e-5,60.0, 0.9,2.0, 1.7,-0.6, 12.0, 2.3),
          (2.5,2e-6,150.0,-1.4,2.0,-0.8, 1.1, -7.0,-1.1),
          (0.7,1e-5,40.0,  2.2,2.0, 2.4, 0.9,  5.0, 0.4)]:
    ex,A,B,Cc=run(*a)
    print("   k+=%.1f kappa=%+.1f"%(a[0],a[3]))
    print("      exact                                         %s"%mp.nstr(ex,12))
    print("   (a) p+ subtraction kept + log[(V-k+)/Lam]        %s   diff %s"%(mp.nstr(A,12),mp.nstr(abs(ex-A),3)))
    print("   (b) literal xi/[1-xi]_+   + log[(V-k+)/Lam]      %s   diff %s"%(mp.nstr(B,12),mp.nstr(abs(ex-B),3)))
    print("   (c) literal xi/[1-xi]_+   + log[k+/Lam]          %s   diff %s"%(mp.nstr(Cc,12),mp.nstr(abs(ex-Cc),3)))
print()
print("   (a) and (c) are both right (residual = O(Lambda/k+) from extending the")
print("   regular pieces to xi=1).  (b), the literal symbol with the p+ log, is not.")
print()
print("   image of [1/p+]_+ in xi :  int dxi [F(xi) - F(1)] / (xi xibar)")
print("   the symbol xi/[1-xi]_+  :  int dxi [F(xi)/xi - F(1)] / (1-xi)")
print("   they differ by F(1) log(V/k+).")
