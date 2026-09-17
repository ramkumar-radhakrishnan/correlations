"""Row 9: all p+ master integrals in explicit closed form, with numerical checks."""
import numpy as np
from scipy.integrate import quad

def masters(p, k, a, ap):
    L = lambda q: np.log(np.abs(q)); da = ap - a; out = {}
    out["T2T3'  h=1"]       = L((p+a)/(p+ap))/da
    out["T2T2'  h=p/k"]     = (-a*L(p+a) + ap*L(p+ap))/(k*da)
    out["T3T3'  h=k/p"]     = k*( L(p)/(a*ap) + L(p+a)/(a*(a-ap)) + L(p+ap)/(ap*(ap-a)) )
    out["T1T3'  h=k/(p+k)"] = k*( L(p+k)/((a-k)*(ap-k)) + L(p+a)/((k-a)*(ap-a))
                                  + L(p+ap)/((k-ap)*(a-ap)) )
    out["T1T2'  h=p/(p+k)"] = ( -k*L(p+k)/((a-k)*(ap-k)) - a*L(p+a)/((k-a)*(ap-a))
                                - ap*L(p+ap)/((k-ap)*(a-ap)) )
    A2 = -k/((a-k)*(ap-k)); A1 = (a*ap - k*k)/((a-k)**2*(ap-k)**2)
    Bc = -a/((k-a)**2*(ap-a)); Bp = -ap/((k-ap)**2*(a-ap))
    out["T1T1'  h=pk/(p+k)^2"] = k*( -A2/(p+k) + A1*L(p+k) + Bc*L(p+a) + Bp*L(p+ap) )
    return out

hfun = {"T2T3'  h=1":          lambda p,k: 1.0+0*p,
        "T2T2'  h=p/k":        lambda p,k: p/k,
        "T3T3'  h=k/p":        lambda p,k: k/p,
        "T1T3'  h=k/(p+k)":    lambda p,k: k/(p+k),
        "T1T2'  h=p/(p+k)":    lambda p,k: p/(p+k),
        "T1T1'  h=pk/(p+k)^2": lambda p,k: p*k/(p+k)**2}

k,A,B,Ap,Bp_,Lam,V = 1.0, 0.8, 1.7, 1.3, 0.6, 1e-8, 1e7
a, ap = k*B/A, k*Bp_/Ap
lo, hi = Lam, V-k
F0, F1 = masters(lo,k,a,ap), masters(hi,k,a,ap)
Q = lambda f,l,h: quad(lambda u: f(np.exp(u),k)*np.exp(u)/((np.exp(u)+a)*(np.exp(u)+ap)),
                       np.log(l), np.log(h), limit=900)[0]/(A*Ap)
print("="*76); print("N2: six p+ masters,  M = (1/(A A'))[F(V-k+) - F(Lambda)]"); print("="*76)
print("   k+=1, A=0.8, B=1.7, A'=1.3, B'=0.6, Lambda=1e-8, V=1e7")
print("   alpha = k+B/A = %.6f,  alpha' = k+B'/A' = %.6f\n" % (a,ap))
print("  %-22s %17s %17s %10s" % ("structure","quadrature","closed form","abs.diff"))
for nm in hfun:
    num = Q(hfun[nm], lo, hi); cf = (F1[nm]-F0[nm])/(A*Ap)
    print("  %-22s %17.9f %17.9f %10.2e" % (nm, num, cf, abs(num-cf)))
print()
print("="*76); print("WHICH MASTER CARRIES WHAT"); print("="*76)
print("  %-22s %22s %18s" % ("structure","d/dlog(1/Lambda)","predicted"))
for nm in hfun:
    d = (Q(hfun[nm],1e-10,hi) - Q(hfun[nm],1e-6,hi))/np.log(1e4)
    pr = 1/(k*B*Bp_) if "k/p" in nm else 0.0
    print("  %-22s %22.9f %18.9f" % (nm, d, pr))
print()
print("  %-22s %22s %18s" % ("structure","d/dlog V","predicted"))
for nm in hfun:
    d = (Q(hfun[nm],lo,1e9) - Q(hfun[nm],lo,1e5))/np.log(1e4)
    pr = 1/(k*A*Ap) if nm.startswith("T2T2") else 0.0
    print("  %-22s %22.9f %18.9f" % (nm, d, pr))
print()
print("  => RAPIDITY log 1/Lambda lives only in T3T3' (the (S/p+)x(S/p+) term),")
print("     residue 1/[k+ (x-w)^2 (x'-w')^2].")
print("  => log V lives only in T2T2' (the (S/k+)x(S/k+) term),")
print("     residue 1/[k+ (x-z)^2 (x'-z)^2].   All four others are V-stable and Lambda-safe.")
