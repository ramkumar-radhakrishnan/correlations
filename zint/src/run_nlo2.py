"""Every number for the note on the four-WW-kernel row.  -> ../RESULTS_nlo2.txt"""
import numpy as np, sympy as sp, sys
sys.path.insert(0,'/home/user/correlations/zint/src')
from nlo2_master import Bten
from quad2d import integrate
from masters import Kin, A7
out=[]
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)
g=0.5772156649015329

p("A. the bracket collapses to three products of two scalar dipole kernels  (sympy, exact)")
yp,z,xp,x,y,w = [sp.Matrix(sp.symbols(f'{n}1 {n}2', real=True))
                 for n in ('yp','z','xp','xx','yy','ww')]
xi=sp.symbols('xi', positive=True)
Kf=lambda a,b: (a-b)/(((a-b).T*(a-b))[0,0])
K1,K2,K3,K4 = Kf(yp,z),Kf(xp,yp),Kf(y,w),Kf(x,z)
dot=lambda u,v:(u.T*v)[0,0]
br=0
for m in range(2):
  for kp in range(2):
    for i in range(2):
      for j in range(2):
        d=lambda a,b: 1 if a==b else 0
        co=(1/(1+xi))*d(kp,m)*d(i,j) - (1/xi)*d(j,m)*d(i,kp) - d(i,m)*d(j,kp)
        if co!=0: br += K1[m]*K2[kp]*K3[i]*K4[j]*co
T=dot(K1,K2)*dot(K3,K4)/(1+xi) - dot(K1,K4)*dot(K2,K3)/xi - dot(K1,K3)*dot(K2,K4)
p("   bracket == T1 + T2 + T3 :", sp.simplify(br-T)==0)
p("")
p("B. power counting.  Each kernel is 1/|sep|; no pair of points appears in two")
p("   kernels, so nowhere is there a 1/sep^2.  In 2D  int d^2(sep) ~ int rho drho/rho")
p("   converges, so the z and y' integrations are ULTRAVIOLET FINITE.")
p("   (Contrast the earlier row, where the light-cone vertex supplied a second")
p("    1/(x-z)^2 on top of the kernel (z-x)^m/(z-x)^2, giving 1/r^2 and a log.)")
p("")
p("C. the tensor master  B^{mj}(Q,c) = int d^2r e^{-iQ.r} (r^m/r^2)((c-r)^j/(c-r)^2)")
Q=np.array([0.9,0.6]); c=np.array([0.7,-0.4]); B=Bten(Q,c)
p(f"   Q = {Q}   c = {c}")
p("   alpha-representation vs direct 2D quadrature:")
for (m,j) in [(0,0),(0,1),(1,0),(1,1)]:
    num=integrate(lambda r,m=m,j=j: r[...,m]*(c-r)[...,j]/(np.sum(r*r,-1)*np.sum((c-r)**2,-1)),
                  -Q, split=(np.hypot(*c),), ntheta=4096, nlev=20, R0=12.0)
    p(f"      B^{m}{j}: rep {B[m,j]:+.9f}   quad {num:+.9f}   rel {abs(B[m,j]-num)/abs(num):.1e}")
p("   (the quadrature is tail-limited here: the integrand falls only as 1/r^2 with a")
p("    nonzero angular average.  The trace gives an exact check instead:)")
a7=A7(Kin(c,Q,1.0))
p(f"      trace B  = {np.trace(B):+.12f}")
p(f"      -conj(A7)= {-np.conj(a7):+.12f}    rel {abs(np.trace(B)+np.conj(a7))/abs(a7):.1e}")
p("      (A7 is the master already verified against quadrature in the first note)")
p("")
p("D. the exact separation")
p("   B^{mj} = -(pi/2) d^{mj} log(4 e^{-2gamma}/(Q^2 c^2))          <- the logarithm")
p("            - pi d^{mj} + pi Q^m Q^j/Q^2 + pi c^m c^j/c^2        <- constants")
p("            + Bcal^{mj}(Q,c)                                     <- remainder, O(|Q||c|)")
def pred(Q,c):
    Qn=np.hypot(*Q); cn=np.hypot(*c)
    return (-np.pi/2*np.eye(2)*np.log(4*np.exp(-2*g)/(Qn**2*cn**2))
            - np.pi*np.eye(2) + np.pi*np.outer(Q,Q)/(Q@Q) + np.pi*np.outer(c,c)/(c@c))
p("   |Q||c|      max |Bcal^{mj}|     max |B^{mj}|")
for sc in (1.0,0.3,0.1,0.03,0.01,0.003,0.001):
    Qs=sc*Q; Bs=Bten(Qs,c); Qn=np.hypot(*Qs); cn=np.hypot(*c)
    p(f"   {Qn*cn:9.5f}   {np.abs(Bs-pred(Qs,c)).max():13.3e}   {np.abs(Bs).max():13.5f}")
p("   trace of the log+constant part = -pi log(4 e^{-2gamma}/(Q^2 c^2)) :")
for sc in (0.1,0.01,0.001):
    Qs=sc*Q; t=np.trace(Bten(Qs,c)); Qn=np.hypot(*Qs); cn=np.hypot(*c)
    p(f"      |Q||c|={Qn*cn:8.5f}: trace {t:+.6f}   predicted {-np.pi*np.log(4*np.exp(-2*g)/(Qn**2*cn**2)):+.6f}")
open('/home/user/correlations/zint/RESULTS_nlo2.txt','w').write("\n".join(out)+"\n")
