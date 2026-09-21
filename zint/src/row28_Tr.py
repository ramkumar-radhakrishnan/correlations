"""Row 28: check the r,s,P,xi rewriting T(r), and extract the UV coefficient at r -> 0."""
import numpy as np
s_=lambda v:v@v; dd=np.eye(2)
def Wjk(y,xv,zv,Pv,Kv):
    A=s_(y-xv); Bv=s_(y-zv); C=s_(xv-zv)
    return ( dd*(A-Bv)/(2*C)
            +(Kv/(Kv-Pv))*( np.outer(y-xv,xv-zv)/C - np.outer(y-xv,y-zv)/(2*Bv) )
            +(Kv/Pv)     *( np.outer(xv-zv,y-zv)/C + np.outer(y-xv,y-zv)/(2*A) ) )
def Rim(y,xv,zv,Pv,Kv):
    W=Wjk(y,xv,zv,Pv,Kv)
    return dd*np.trace(W) - (Kv/Pv)*W.T - (Kv/(Kv-Pv))*W
def exact_T(r,s,Pv_,xi,dperp=2.0):
    """xi xibar x (P^i r^m /P^2 r^2) R^{im} , in the r,s,P variables (k+=1)"""
    K=1.0; P=xi*K
    xv=np.zeros(2); zv=r; y=s                      # x at origin => r=z-x=z , s=y-x=y
    R=Rim(y,xv,zv,P,K)
    return xi*(1-xi)*np.einsum('i,m,im->',Pv_/s_(Pv_),r/s_(r),R)
def yours_T(r,s,Pv_,xi,dperp=2.0):
    xb=1-xi; P2=s_(Pv_); r2=s_(r); s2=s_(s); sr=s@r; Pr=Pv_@r; Ps=Pv_@s; smr=s-r; smr2=s_(smr)
    t14=dperp*xi*xb/2*( 2*Pr*sr/(P2*r2**2) - Pr/(P2*r2) )
    t15=(xi/xb)*( Ps/(P2*r2) + Ps*(sr-r2)/(2*P2*r2*smr2) )
    t16=xb*(Pr/(P2*r2))*( 2 - 2*sr/r2 - sr/(2*s2) )
    t18=(xb/xi)*((Ps-Pr)/(P2*r2))*( 1 - sr/s2 )
    t19=-xi*(Pr/P2)*( 2*sr/r2**2 + (s2-sr)/(2*r2*smr2) - 1/(2*r2) )
    t110=( 2*Pr*sr/(P2*r2**2) - Pr/(P2*r2) + Ps/(2*P2*s2) - Ps*sr/(2*P2*r2*s2)
          + sr/(2*P2*r2*smr2)*(Ps-Pr) )
    return t14+t15+t16+t18+t19+t110
rng=np.random.default_rng(11)
print("="*80); print("1.  YOUR T(r) vs THE EXACT CONTRACTION"); print("="*80)
print("   %-30s %16s %16s %10s"%("(xi, random r,s,P)","yours","exact","diff"))
for t in range(5):
    r,s,Pv_=[rng.normal(size=2)*1.1 for _ in range(3)]; xi=rng.uniform(0.2,0.8)
    a=yours_T(r,s,Pv_,xi); b=exact_T(r,s,Pv_,xi)
    print("   xi=%.3f                         %16.8f %16.8f %10.2e"%(xi,a,b,abs(a-b)))
print()
print("   NOTE on eq (1.9) of your note: the middle term reads (s^2 - s.r)/[2(s-r)^2] ,")
print("   which is DIMENSIONLESS while its neighbours are 1/length^2.  It must be")
print("   (s^2 - s.r)/[2 r^2 (s-r)^2] ; that is what is used above.")

print(); print("="*80); print("2.  THE UV AT r -> 0 :  <T(r)>_rhat x r^2"); print("="*80)
print("   %-8s %14s %14s %14s   %s"%("xi","rho=1e-2","1e-4","1e-6","-> c_UV(xi)"))
s0=np.array([0.7,-0.4]); P0=np.array([-0.3,0.9])
for xi in (0.1,0.25,0.5,0.75,0.9):
    v=[]
    for rho in (1e-2,1e-4,1e-6):
        acc=0.0
        for th in np.linspace(0,2*np.pi,2000,endpoint=False):
            rr=rho*np.array([np.cos(th),np.sin(th)]); acc+=exact_T(rr,s0,P0,xi)
        v.append(acc/2000*rho*rho)
    print("   %-8.2f %14.8f %14.8f %14.8f"%(xi,*v))
print()
print("   compare with candidate xi-structures, normalised at xi = 0.5 :")
def cuv(xi,rho=1e-5,n=4000):
    acc=0.0
    for th in np.linspace(0,2*np.pi,n,endpoint=False):
        rr=rho*np.array([np.cos(th),np.sin(th)]); acc+=exact_T(rr,s0,P0,xi)
    return acc/n*rho*rho
xs=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
c=np.array([cuv(x) for x in xs])
Pgg=lambda z: z/(1-z)+(1-z)/z+z*(1-z)
print("   %-7s %14s %14s %14s %14s"%("xi","c_UV","c_UV/C_UV","c_UV/[xi xibar]","c_UV/1"))
for xx,cc in zip(xs,c):
    print("   %-7.2f %14.8f %14.8f %14.8f %14.8f"%(xx,cc,cc/Pgg(xx),cc/(xx*(1-xx)),cc))
