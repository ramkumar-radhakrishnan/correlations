"""Row 11: p+ integrals for terms I and III -- exact primitives, verified symbolically."""
import sympy as sp, numpy as np, mpmath as mp
mp.mp.dps=30
p,K,A,B,Lam,V = sp.symbols('p kplus A B Lambda V', positive=True)
be = K*B/A

print("="*78); print("TERM I :  integrand  (p+k+)/[(k+-p)^2 (pA - k+B)]"); print("="*78)
fI = (p+K)/((K-p)**2*(p*A-K*B))
FI = -2/((A-B)*(p-K)) + (A+B)/(K*(A-B)**2)*sp.log((p-be)/(p-K))
print("   primitive  F_I =  -2/[(A-B)(p-k+)] + (A+B)/(k+(A-B)^2) log[(p-beta)/(p-k+)]")
print("   dF/dp - integrand =", sp.simplify(sp.diff(FI,p)-fI), "   <- EXACT")
print()
print("   on [Lambda, k+-Lambda] U [k++Lambda, V], principal value at beta = k+B/A:")
print("      1/(p-k+) part  ->  4/(Lambda(A-B)) - 2/(A-B)[1/k+ + 1/(V-k+)]  + O(Lambda)")
print("      log part       ->  (A+B)/(k+(A-B)^2) [ log(A/B) + log((V-beta)/(V-k+)) ]")
print("      the log(Lambda) from the two strip edges CANCELS.")
print()
print("   => T_I = 4/(Lambda (A-B))  -  2/(A-B)[1/k+ + 1/(V-k+)]")
print("            + (A+B)/(k+(A-B)^2)[ log(A/B) + log((V-beta)/(V-k+)) ]   +  O(Lambda)")
print("      plus the on-shell piece   -/+ i pi (A+B)/(k+ (A-B)^2)   from beta.")
FIn = sp.lambdify((p,K,A,B), -2/((A-B)*(p-K)) + (A+B)/(K*(A-B)**2)*sp.log(sp.Abs((p-K*B/A)/(p-K))),'numpy')
def TI_exact(Kn,An,Bn,Ln,Vn):
    return (FIn(Kn-Ln,Kn,An,Bn)-FIn(Ln,Kn,An,Bn)) + (FIn(Vn,Kn,An,Bn)-FIn(Kn+Ln,Kn,An,Bn))
def TI_asym(Kn,An,Bn,Ln,Vn):
    ben=Kn*Bn/An
    return (4/(Ln*(An-Bn)) - 2/(An-Bn)*(1/Kn+1/(Vn-Kn))
            + (An+Bn)/(Kn*(An-Bn)**2)*(np.log(An/Bn)+np.log((Vn-ben)/(Vn-Kn))))
def TI_mp(Kn,An,Bn,Ln,Vn):
    f=lambda t:(t+Kn)/((Kn-t)**2*(t*An-Kn*Bn)); ben=Kn*Bn/An
    segs=[(Ln,Kn-Ln),(Kn+Ln,Vn)]; tot=mp.mpf(0)
    for lo,hi in segs:
        if lo<ben<hi:
            tot+=mp.quad(f,[lo,ben-1e-9])+mp.quad(f,[ben+1e-9,hi])
        else: tot+=mp.quad(f,[lo,hi])
    return tot
print()
print("   %-26s %17s %17s %17s"%("(A,B,Lam,V)","PV quadrature","exact primitive","asymptotic form"))
for An,Bn,Ln,Vn in [(0.36,2.09,1e-4,50.),(1.7,0.40,1e-5,200.),(2.5,1.10,1e-4,80.)]:
    q=float(TI_mp(1.0,An,Bn,Ln,Vn)); e=TI_exact(1.0,An,Bn,Ln,Vn); a=TI_asym(1.0,An,Bn,Ln,Vn)
    print("   %-26s %17.6f %17.6f %17.6f"%("(%.2f,%.2f,%.0e,%.0f)"%(An,Bn,Ln,Vn),q,e,a))

print(); print("="*78); print("TERM III :  (1/k+) integrand (k+B + pA)/[p (pA - k+B)]"); print("="*78)
fIII=(K*B+p*A)/(p*(p*A-K*B))/K
FIII=(-sp.log(p)+2*sp.log(p-be))/K
print("   primitive  F_III = (1/k+)[ -log p + 2 log|p - beta| ]")
print("   dF/dp - integrand =", sp.simplify(sp.diff(FIII,p)-fIII), "   <- EXACT")
print()
print("   T_III = (1/k+)[ -log(V/Lambda) + 2 log| (V-beta)/(Lambda-beta) | ]     [PV, exact]")
print("        -> (1/k+) log( V Lambda A^2/(k+^2 B^2) )   as Lambda->0, V->infinity")
print("      plus the on-shell piece  -/+ 2 i pi / k+ .")
def TIII_ex(Kn,An,Bn,Ln,Vn):
    ben=Kn*Bn/An
    return (-np.log(Vn/Ln)+2*np.log(abs((Vn-ben)/(Ln-ben))))/Kn
def TIII_mp(Kn,An,Bn,Ln,Vn):
    f=lambda t:(Kn*Bn+t*An)/(t*(t*An-Kn*Bn))/Kn; ben=Kn*Bn/An
    return mp.quad(f,[Ln,ben-1e-9])+mp.quad(f,[ben+1e-9,Vn])
print()
print("   %-26s %17s %17s %17s"%("(A,B,Lam,V)","PV quadrature","exact primitive","asymptotic form"))
for An,Bn,Ln,Vn in [(0.36,2.09,1e-6,50.),(1.7,0.40,1e-7,200.),(2.5,1.10,1e-6,80.)]:
    q=float(TIII_mp(1.0,An,Bn,Ln,Vn)); e=TIII_ex(1.0,An,Bn,Ln,Vn)
    a=np.log(Vn*Ln*An**2/(1.0*Bn**2))
    print("   %-26s %17.6f %17.6f %17.6f"%("(%.2f,%.2f,%.0e,%.0f)"%(An,Bn,Ln,Vn),q,e,a))
