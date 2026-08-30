"""G(r) = int dxi [bracket]/(xibar D) = l_k*G_l + ln(B/A)*G_L + G_rat,
   and its UV (r->0) separation."""
import numpy as np
from scipy.integrate import quad
from xi_masters import J

s=np.array([0.7,-0.4]); P=np.array([1.3,0.5]); s2=s@s; P2=P@P; sP=s@P; dperp=2.0

def Vlist(r,kap=0.5):
    r2=np.sum(r**2,-1); u=r-s; u2=np.sum(u**2,-1); Pr=r@P; sr=r@s
    return np.array([dperp/2*(2*Pr*sr/(P2*r2**2)-Pr/(P2*r2)),
        sP/(P2*r2)+sP*(sr-r2)/(2*P2*r2*u2),
        Pr/(P2*r2)*(2-2*sr/r2-sr/(2*s2)),
        (sP-Pr)/(P2*r2)*(1-kap*sr/s2),
        -Pr/P2*(2*sr/r2**2+(s2-sr)/(2*r2*u2)-1/(2*r2)),
        2*Pr*sr/(P2*r2**2)-Pr/(P2*r2)+sP/(2*P2*s2)-sP*sr/(2*P2*r2*s2)
          +sr*(sP-Pr)/(2*P2*r2*u2)])

def G(r,lk,kap=0.5,parts=False):
    A=float(np.sum((s-r)**2)); B=s2; C=B-A; L=np.log(B/A)
    V=Vlist(r,kap); Jd=J(A,B,lk)
    Gl  = V[1]/B + V[3]/A
    GL  = -A*B/C**3*V[0] - A/(B*C)*V[1] + B/C**2*V[2] - B/(A*C)*V[3] - A/C**2*V[4] + V[5]/C
    Gr  = (1/(2*C)+A/C**2)*V[0] - V[2]/C + V[4]/C
    tot = sum(Jd[i+1]*V[i] for i in range(6))
    if parts: return tot, Gl, GL, Gr, L
    return tot

lk=np.log(1e9)
print("A. G(r) assembled from J_i   vs   direct xi-quadrature of the bracket/(xibar D)")
for rv in [np.array([0.31,0.22]), np.array([-0.5,0.8]), np.array([0.12,-0.05])]:
    for kap in (0.5,1.0):
        A=float(np.sum((s-rv)**2)); B=s2
        w=[lambda x:x*(1-x), lambda x:x/(1-x), lambda x:1-x, lambda x:(1-x)/x,
           lambda x:x, lambda x:1.0]
        V=Vlist(rv,kap)
        num=sum(V[i]*quad(lambda x,f=w[i]:f(x)/((1-x)*A+x*B),1e-9,1-1e-9,limit=400)[0]
                for i in range(6))
        tot,Gl,GL,Gr,L = G(rv,lk,kap,parts=True)
        print(f"  r={rv} kappa={kap}: assembled {tot:+.10f}  quad {num:+.10f}"
              f"   rel {abs(tot-num)/abs(num):.1e}")
        print(f"      = l_k*({Gl:+.7f}) + ln(B/A)*({GL:+.7f}) + ({Gr:+.7f})"
              f"   [l_k={lk:.4f}, ln(B/A)={L:+.5f}]")
print()
print("B. UV separation:  r^2 * <G>_angle  ->  (2 l_k - 11/6) (s.P)/(P^2 s^2)")
th=2*np.pi*np.arange(8192)/8192
for kap in (0.5,1.0):
    for R in (1e-2,1e-3,1e-4):
        rv=np.stack([R*np.cos(th),R*np.sin(th)],-1)
        vals=np.array([G(rr,lk,kap) for rr in rv[::32]])
        print(f"  kappa={kap} R={R:.0e}: <G> r^2 = {vals.mean()*R**2:+.8f}"
              f"   predicted {(2*lk-11/6)*sP/(P2*s2):+.8f}")

print()
print("C. how the UV coefficient 2 l_k - 11/6 is distributed over the three blocks")
for kap in (0.5,1.0):
    R=1e-3; rv=np.stack([R*np.cos(th),R*np.sin(th)],-1)[::16]
    gl=np.array([G(rr,lk,kap,parts=True)[1] for rr in rv])
    gL=np.array([G(rr,lk,kap,parts=True)[2]*G(rr,lk,kap,parts=True)[4] for rr in rv])
    gr=np.array([G(rr,lk,kap,parts=True)[3] for rr in rv])
    n=sP/(P2*s2)
    print(f"  kappa={kap}:  l_k-block  {gl.mean()*R**2/n:+.6f}  (expect 2)")
    print(f"            ln(B/A)-block {gL.mean()*R**2/n:+.6f}")
    print(f"            rational      {gr.mean()*R**2/n:+.6f}")
    print(f"            log+rational  {(gL.mean()+gr.mean())*R**2/n:+.6f}  (expect -11/6 = {-11/6:.6f})")
