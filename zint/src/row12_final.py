"""Row 12: the large-|z| coefficient, and the three p+ masters."""
import numpy as np
from scipy.integrate import quad
K=1.0; Pp=0.37; S=Pp+K
s=lambda v:v@v
x=np.array([0.31,-0.77]); y=np.array([0.62,0.41]); z0=np.array([0.05,-0.23])
w=np.array([1.13,0.42]); xp=np.array([-0.52,0.19]); wp=np.array([-0.31,0.88])
def integ(z,P=Pp):
    A=s(x-z); B=s(x-w); C=s(z-w); Ss=P+K
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)
    br=np.eye(2)*(A-B)/(2*Ss*C) + Pm/P + Km/K
    Ki=(xp-wp)/s(xp-wp); Kj=(y-z)/s(y-z)
    return np.einsum('ij,i,j->',br,Ki,Kj)/(P*A+K*B)
print("="*78); print("1.  WHERE THE DIVERGENCE IS:  |z| -> infinity"); print("="*78)
print("   mechanism: the one-rho vertex degenerates,")
print("      bracket/[p+(x-z)^2+k+(x-w)^2]  ->  -(1/(2 p+ k+)) (x-w)^i (x-z)^j/[(x-w)^2(x-z)^2]")
print("   and (y-z)^j/(y-z)^2 -> -zhat^j/|z| , (x-z)^j/(x-z)^2 -> -zhat^j/|z| , so the")
print("   j-contraction gives 1/|z|^2 with NO angular cancellation.")
print()
pred = -1/(2*Pp*K)*((xp-wp)@(x-w))/(s(xp-wp)*s(x-w))
print("   predicted  -(1/(2 p+ k+)) (x'-w').(x-w)/[(x'-w')^2 (x-w)^2] = %.8f" % pred)
print("   %-10s %18s" % ("|z|","angular avg x |z|^2"))
for R in (1e2,1e3,1e4,1e5,1e6):
    acc=0.0
    for th in np.linspace(0,2*np.pi,2048,endpoint=False):
        acc+=integ(R*np.array([np.cos(th),np.sin(th)]))
    print("   %-10.0e %18.8f" % (R, acc/2048*R**2))
print("   -> flat, and equal to the prediction.  The divergence sits ENTIRELY at large |z|.")
print()
print("   other variables:  x  log divergent too, but x is a rho position -- cut by the")
print("                        support of the source;")
print("                     w, w' log divergent only for the BARE integrand; the phase")
print("                        e^{-ik(w'-w)} supplies Bessel R^{-1/2} damping;")
print("                     y, x' convergent (first angular moments vanish).")
print("   => the one unprotected transverse divergence is |z| -> infinity.  It is")
print("      INFRARED, not UV, and its coefficient is a product of two WW kernels --")
print("      the same structure as the 4-rho term, so combine them before regulating.")

print(); print("="*78); print("2.  THE THREE p+ MASTERS,  D = p+ A + k+ B  (NO on-contour pole:"); print("    both terms are positive, unlike the B_3 row)"); print("="*78)
A=s(x-z0); B=s(x-w); C=s(z0-w); Lam=1e-8; V=1e7; Pv=V-K
Q=lambda f,l,h: quad(lambda t:f(np.exp(t))*np.exp(t),np.log(l),np.log(h),limit=900)[0]
rows=[("J1 = int dp/[p D]      (the 1/p+ bracket term)", lambda t:1/(t*(t*A+K*B)),
       np.log(Pv*K*B/(Lam*(Pv*A+K*B)))/(K*B),
       "(1/(k+B)) log[ (V-k+) k+B / (Lambda((V-k+)A+k+B)) ]  ->  (1/(k+B)) log(k+B/(Lambda A))"),
      ("J2 = int dp/D          (the 1/k+ bracket term)", lambda t:1/(t*A+K*B),
       np.log((Pv*A+K*B)/(Lam*A+K*B))/A,
       "(1/A) log[ ((V-k+)A+k+B)/(Lambda A+k+B) ]            ->  (1/A) log( V A/(k+B) )"),
      ("J3 = int dp/[(p+k+) D] (the delta^ij term)",     lambda t:1/((t+K)*(t*A+K*B)),
       np.log(V*B/(Pv*A+K*B))/(K*(B-A)),
       "(1/(k+(B-A))) log[ V B/((V-k+)A+k+B) ]               ->  (1/(k+(B-A))) log(B/A)")]
print("   A=(x-z)^2=%.4f  B=(x-w)^2=%.4f  C=(z-w)^2=%.4f  Lambda=1e-8  V=1e7\n"%(A,B,C))
for nm,f,cf,form in rows:
    print("   %s"%nm); print("      %s"%form)
    print("      quad %16.9f   closed %16.9f   diff %.1e\n"%(Q(f,Lam,Pv),cf,abs(Q(f,Lam,Pv)-cf)))
print("   the delta^ij term collapses exactly, the (A-B) cancelling the (B-A):")
lhs=(A-B)/(2*C)*rows[2][2]; rhs=-np.log(V*B/(Pv*A+K*B))/(2*K*C)
print("      (A-B)/(2C) x J3 = %.12f   =  -(1/(2 k+ C)) log[V B/((V-k+)A+k+B)] = %.12f"%(lhs,rhs))
print("      -> as V -> infinity:   - delta^ij log[(x-w)^2/(x-z)^2] / [2 k+ (z-w)^2]")
print()
print("   ATTRIBUTION:   rapidity log(1/Lambda)  only in J1, residue 1/[k+ (x-w)^2]")
print("                  log V                   only in J2, residue 1/(x-z)^2")
print("                  the delta^ij term is free of BOTH.")
