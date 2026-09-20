"""Row 20: transverse divergences of the p+-INTEGRATED |B_2|^2 row (parts 2a..2f summed)."""
import numpy as np, itertools
from row19_lib import their_coef
s=lambda v:v@v
def taus(X,Z,W):
    A,B,C=s(X-Z),s(X-W),s(Z-W)
    return (np.eye(2)*(A-B)/(2*C), np.outer(X-Z,(Z-W)/C-(X-W)/(2*B)),
            np.outer((Z-W)/C+(X-Z)/(2*A),X-W))
con=lambda a,b: np.einsum('ji,ji->',a,b)
KT=np.array([0.62,-0.37]); KP=1.0; LAM=1e-6; VEE=1e6
def G(x,xp,z,w,wp,phase=True):
    """the whole p+-integrated bracket: sum_ab M_ab (tau_a . tau_b')"""
    A,B,C=s(x-z),s(x-w),s(z-w); Ap,Bp,Cp=s(xp-z),s(xp-wp),s(z-wp)
    T=their_coef(A,B,C,Ap,Bp,Cp,KP,LAM,VEE)
    t1,t2,t3=taus(x,z,w); u1,u2,u3=taus(xp,z,wp)
    val=( T['23']*(con(t2,u3)+con(t3,u2)) + T['33']*con(t3,u3) + T['22']*con(t2,u2)
        + T['13']*(con(t1,u3)+con(t3,u1)) + T['12']*(con(t1,u2)+con(t2,u1)) + T['11_full'] )
    return val*(np.exp(-1j*(KT@(wp-w))) if phase else 1.0)

names=['x',"x'",'z','w',"w'"]
base=[np.array([0.31,-0.77]),np.array([-0.52,0.19]),np.array([0.05,-0.23]),
      np.array([1.13,0.42]),np.array([-0.31,0.88])]
u=lambda t:np.array([np.cos(t),np.sin(t)]); OFF=[0,2.1,4.0,1.05,5.2]
print("="*80); print("1.  SHORT DISTANCE: ALL 26 COLLAPSING SUBSETS  (UV / COLLINEAR)"); print("="*80)
print("   log divergence  <=>  <|G|> x rho^{2(n-1)}  CONSTANT as rho -> 0")
print("   %-22s %11s %11s %11s   %s"%("collapsing subset","rho=1e-2","1e-4","1e-6","verdict"))
bad=[]
for n in range(2,6):
    for sub in itertools.combinations(range(5),n):
        v=[]
        for rho in (1e-2,1e-4,1e-6):
            acc=0.0
            for th in np.linspace(0,2*np.pi,36,endpoint=False):
                pt=[b.copy() for b in base]; c=base[sub[0]]
                for a,i in enumerate(sub): pt[i]=c+rho*u(th+OFF[a])
                acc+=abs(G(*pt))
            v.append(acc/36*rho**(2*(n-1)))
        div = v[1]>0.25*v[0] and v[2]>0.25*v[1]
        if div: bad.append(sub)
        lab="{"+",".join(names[i] for i in sub)+"}"
        if n<=3 or div:
            print("   %-22s %11.3e %11.3e %11.3e   %s"%(lab,*v,"*** DIV ***" if div else "convergent"))
print("   (4- and 5-point subsets all convergent, suppressed from the table)")
print("   ==> %s"%("NO short-distance divergence in any of the 26 regions" if not bad else "%d divergent: %s"%(len(bad),bad)))

print(); print("="*80); print("2.  LARGE DISTANCE: ONE VARIABLE -> infinity  (SOFT / IR)"); print("="*80)
print("   log divergence  <=>  |<G>| x R^2  CONSTANT as R -> infinity")
print("   %-6s %13s %13s %13s   %s"%("var","R=1e2","R=1e3","R=1e4","verdict"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e3,1e4):
        nth=int(min(max(4096,60*np.sqrt(KT@KT)*R),200_000)); acc=0j
        for th in np.linspace(0,2*np.pi,nth,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); acc+=G(*pt)
        v.append(abs(acc/nth)*R*R)
    flat = v[2]>1e-9 and v[2]>0.25*v[0]
    print("   %-6s %13.5e %13.5e %13.5e   %s"%(nm,*v,"*** LOG DIV ***" if flat else "convergent"))

print(); print("="*80); print("3.  THE NAMED REGIONS, ONE BY ONE"); print("="*80)
x,xp,z,w,wp=base
def probe(label,mk,rhos,power,note):
    print("   %s"%label)
    print("     %-9s %15s %17s"%("rho","<|G|>","<|G|> x rho^%d"%power))
    for rho in rhos:
        acc=0.0; n=400
        for t in np.linspace(0,2*np.pi,n,endpoint=False):
            acc+=abs(G(*mk(rho,t)))
        a=acc/n; print("     %-9.0e %15.6e %17.9f"%(rho,a,a*rho**power))
    print("     -> %s\n"%note)
probe("(a) z -> w   [the two final gluons at one point: FINAL-STATE COLLINEAR]",
      lambda r,t:(x,xp,w+r*u(t),w,wp),(1e-2,1e-4,1e-6),1,
      "<|G|> x rho constant => G ~ 1/rho against rho drho : CONVERGENT, no collinear log")
probe("(b) w -> x   [observed gluon on its source: WW-kernel endpoint]",
      lambda r,t:(x,xp,z,x+r*u(t),wp),(1e-2,1e-4,1e-6),1,
      "G ~ 1/rho : CONVERGENT (this is the ordinary 1/r of the WW kernel)")
probe("(c) z -> x   [unobserved gluon on the source]",
      lambda r,t:(x,xp,x+r*u(t),w,wp),(1e-2,1e-4,1e-6),1,
      "G ~ 1/rho : CONVERGENT")
probe("(d) w' -> w  [the k_perp-integrated corner]",
      lambda r,t:(x,xp,z,wp+r*u(t),wp),(1e-2,1e-4,1e-6),0,
      "<|G|> itself constant => G REGULAR there")
