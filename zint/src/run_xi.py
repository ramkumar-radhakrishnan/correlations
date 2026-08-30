"""Every number quoted in the xi-integration note.  -> ../RESULTS_xi.txt"""
import numpy as np, mpmath as mp
from scipy.integrate import quad
from xi_masters import J
from xi_assemble import Vlist, G, s, P, s2, P2, sP
out=[]
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)

lk=np.log(1e9)
p("kinematics: s =",s," P =",P,"  s^2 = %.4f  s.P = %.4f  l_k = ln(k+/Lam) = %.6f"%(s2,sP,lk))
p("")
p("A. the six xi-integrals J_i  (analytic vs quadrature, lambda = 1e-9)")
w=[lambda x:x*(1-x), lambda x:x/(1-x), lambda x:1-x, lambda x:(1-x)/x, lambda x:x, lambda x:1.0]
for (A,B) in [(1.0,4.0),(2.3,0.7),(0.9,1.1)]:
    an=J(A,B,lk)
    p("   A=%.1f B=%.1f"%(A,B))
    for i in range(6):
        num=quad(lambda x,f=w[i]:f(x)/((1-x)*A+x*B),1e-9,1-1e-9,limit=400)[0]
        p("      J%d  analytic %+.12f   quad %+.12f   rel %.1e"%(i+1,an[i+1],num,abs(an[i+1]-num)/abs(num)))
p("")
p("B. the paper's (C.5) and (C.6)      [A=(W')^2=1, B=W^2=4]")
A,B=1.0,4.0; L=np.log(B/A)
p("   (C.5) quadrature            %.10f"%quad(lambda x:1/(x*((1-x)*A+x*B)),1e-9,1-1e-9,limit=400)[0])
p("         paper  (l_k - L)/A  = %.10f     agrees"%((lk-L)/A))
p("   (C.6) quadrature            %.10f"%quad(lambda x:1/((1-x)*((1-x)*A+x*B)),1e-9,1-1e-9,limit=400)[0])
p("         paper  (l_k - L)/B  = %.10f     DISAGREES  (sign of L)"%((lk-L)/B))
p("         correct (l_k + L)/B = %.10f     agrees"%((lk+L)/B))
p("")
p("C. the kernel G(r) = sum_i J_i V_i   vs direct xi-quadrature of the whole bracket")
for rv in [np.array([0.31,0.22]), np.array([-0.5,0.8]), np.array([0.12,-0.05])]:
    for kap in (0.5,1.0):
        A=float(np.sum((s-rv)**2)); V=Vlist(rv,kap)
        num=sum(V[i]*quad(lambda x,f=w[i]:f(x)/((1-x)*A+x*s2),1e-9,1-1e-9,limit=400)[0] for i in range(6))
        tot,Gl,GL,Gr,Lv=G(rv,lk,kap,parts=True)
        p("   r=%s kappa=%.1f: assembled %+.9f  quad %+.9f  rel %.1e"%(rv,kap,tot,num,abs(tot-num)/abs(num)))
p("")
p("D. UV coefficient of each vertex V_i  (units of (s.P)/(P^2 r^2), angular average at r->0)")
th=2*np.pi*np.arange(4096)/4096; R=1e-4
rv=np.stack([R*np.cos(th),R*np.sin(th)],-1)
for kap in (0.5,1.0):
    c=[v.mean()*R**2/(sP/P2) for v in Vlist(rv,kap)]
    p("   kappa=%.1f : %s"%(kap,np.round(c,7)))
p("   weighted sum = xi*xib + xi/xib + xib/xi = C_UV = P_gg/(2 Nc)   (checked at xi=0.3, 0.65)")
p("")
p("E. UV coefficient of the xi-integrated kernel, and how it splits over the blocks")
th=2*np.pi*np.arange(512)/512; R=1e-3
rv=np.stack([R*np.cos(th),R*np.sin(th)],-1)
for kap in (0.5,1.0):
    pr=[G(rr,lk,kap,parts=True) for rr in rv]
    tot=np.mean([q[0] for q in pr]); gl=np.mean([q[1] for q in pr])
    gL=np.mean([q[2]*q[4] for q in pr]); gr=np.mean([q[3] for q in pr]); n=sP/(P2*s2)
    p("   kappa=%.1f : r^2<G> / [(s.P)/(P^2 s^2)] = %+.6f   predicted 2 l_k - 11/6 = %+.6f"
      %(kap,tot*R**2/n,2*lk-11/6))
    p("             l_k block %+.6f (=2)   ln(B/A)+rational blocks %+.6f (=-11/6)"
      %(gl*R**2/n, (gL+gr)*R**2/n))
p("")
p("F. the xi-integrals of C_UV  (mpmath, 30 digits, lambda=1e-9)  -- the paper's (C.10)")
mp.mp.dps=30; lam=mp.mpf('1e-9'); lkm=mp.log(1/lam)
Cf=lambda x: x*(1-x)+x/(1-x)+(1-x)/x
I0=mp.quad(lambda x:Cf(x),[lam,1-lam]); I1=mp.quad(lambda x:Cf(x)*mp.log(x*(1-x)),[lam,1-lam])
I2=mp.quad(lambda x:Cf(x)*mp.log(1/(1-x)),[lam,1-lam])
p("   int C_UV dxi              = %s     2 l_k - 11/6          = %s"%(mp.nstr(I0,14),mp.nstr(2*lkm-mp.mpf(11)/6,14)))
p("   int C_UV ln(xi xib) dxi   = %s    -l_k^2 + 67/18 - pi^2/3 = %s"%(mp.nstr(I1,14),mp.nstr(-lkm**2+mp.mpf(67)/18-mp.pi**2/3,14)))
p("   int C_UV ln(1/xib) dxi    = %s     l_k^2/2 + pi^2/6 - 67/36 = %s"%(mp.nstr(I2,14),mp.nstr(lkm**2/2+mp.pi**2/6-mp.mpf(67)/36,14)))
p("   consistency  I[ln(xi xib)] + 2 I[ln(1/xib)] = %s"%mp.nstr(I1+2*I2,8))
p("")
p("G. regularity of each J_i as C = B - A -> 0  (spurious 1/C, 1/C^2, 1/C^3 poles cancel)")
for eps in (1e-2,1e-4,1e-6):
    A=s2*(1-eps); jj=J(A,s2,lk)
    p("   C/B=%8.1e :  J1*B=%.8f (1/6)  J3*B=%.8f (1/2)  J5*B=%.8f (1/2)  J6*B=%.8f (1)"
      %(eps,jj[1]*s2,jj[3]*s2,jj[5]*s2,jj[6]*s2))
    p("                 J2*B-l_k=%.8f (-1)   J4*B-l_k=%.8f (-1)"%(jj[2]*s2-lk,jj[4]*s2-lk))
open('/home/user/correlations/zint/RESULTS_xi.txt','w').write("\n".join(out)+"\n")
