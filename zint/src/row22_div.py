"""Row 22: transverse divergences of the reduced two-C row (before the p+ integration)."""
import numpy as np, itertools
s=lambda v:v@v; KT=np.array([0.62,-0.37]); KP=1.0; PP=0.37
def integ(x,xp,y,yp,z,P=PP,k=KP,phase=True):
    S=P+k; b1=2.0*P/S**2-2/k; b2=P/k**2+1/P; b3=2/k
    N=(y-z)/s(y-z); Np=(yp-z)/s(yp-z); X=(x-y)/s(x-y); Xp=(xp-yp)/s(xp-yp)
    val=b1*(N@X)*(Np@Xp)+b2*(N@Np)*(X@Xp)+b3*(N@Xp)*(Np@X)
    return val*(np.exp(-1j*(KT@(yp-y))*S/k) if phase else 1.0)
names=['x',"x'",'y',"y'",'z']
base=[np.array([0.41,-0.63]),np.array([-0.58,0.22]),np.array([0.09,0.51]),
      np.array([1.07,-0.34]),np.array([-0.25,0.86])]
u=lambda t:np.array([np.cos(t),np.sin(t)]); OFF=[0,2.1,4.0,1.05,5.2]
print("="*80); print("1.  SHORT DISTANCE: ALL 26 COLLAPSING SUBSETS"); print("="*80)
print("   %-20s %11s %11s %11s   %s"%("subset","rho=1e-2","1e-4","1e-6","verdict"))
bad=[]
for n in range(2,6):
    for sub in itertools.combinations(range(5),n):
        v=[]
        for rho in (1e-2,1e-4,1e-6):
            acc=0.0
            for th in np.linspace(0,2*np.pi,36,endpoint=False):
                pt=[b.copy() for b in base]; c=base[sub[0]]
                for a,i in enumerate(sub): pt[i]=c+rho*u(th+OFF[a])
                acc+=abs(integ(*pt))
            v.append(acc/36*rho**(2*(n-1)))
        dv=v[1]>0.25*v[0] and v[2]>0.25*v[1]
        if dv: bad.append(sub)
        if n<=3 or dv:
            print("   %-20s %11.3e %11.3e %11.3e   %s"%("{"+",".join(names[i] for i in sub)+"}",*v,
                  "*** DIV ***" if dv else "convergent"))
print("   (4- and 5-point subsets all convergent)")
print("   ==> %s"%("NO short-distance divergence" if not bad else "%d divergent"%len(bad)))
print()
print("   the two kernel endpoints, in detail:")
for lab,mk in (("z -> y   (C kernel)",lambda r,t:(base[0],base[1],base[2],base[3],base[2]+r*u(t))),
               ("y -> x   (A kernel, the collinear candidate)",
                lambda r,t:(base[0],base[1],base[0]+r*u(t),base[3],base[4]))):
    print("     %s"%lab)
    for rho in (1e-2,1e-4,1e-6):
        acc=0.0
        for t in np.linspace(0,2*np.pi,300,endpoint=False): acc+=abs(integ(*mk(rho,t)))
        print("       rho=%-8.0e  <|I|> x rho = %14.9f"%(rho,acc/300*rho))
    print("       -> ~ 1/rho against rho drho : convergent.  A log needs 1/rho^2, i.e. the")
    print("          amplitude AND conjugate kernel singular on the SAME separation.")

print(); print("="*80); print("2.  LARGE DISTANCE"); print("="*80)
print("   %-6s %13s %13s %13s   %s"%("var","R=1e2","R=1e3","R=1e4","verdict"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e3,1e4):
        nth=int(min(max(4096,60*np.sqrt(KT@KT)*R),200_000)); acc=0j
        for th in np.linspace(0,2*np.pi,nth,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); acc+=integ(*pt)
        v.append(abs(acc/nth)*R*R)
    print("   %-6s %13.5e %13.5e %13.5e   %s"%(nm,*v,
          "*** LOG DIV ***" if (v[2]>1e-9 and v[2]>0.25*v[0]) else "convergent"))
print()
x,xp,y,yp,_=base; S=PP+KP
b1=2.0*PP/S**2-2/KP; b2=PP/KP**2+1/PP; b3=2/KP
X=(x-y)/s(x-y); Xp=(xp-yp)/s(xp-yp)
pred=(b1+b3)*0.5*(X@Xp)+b2*(X@Xp)     # <N N'> -> d^mm'/(2|z|^2) etc.
print("   large-|z| coefficient  (angle-averaged, phase stripped):")
print("     <N^m N'^m'> -> delta^{mm'}/(2|z|^2)  =>  coefficient = [(b1+b3)/2 + b2] (X.X')")
acc=0.0; n=6000
for t in np.linspace(0,2*np.pi,n,endpoint=False): acc+=integ(x,xp,y,yp,3e4*u(t))
print("     measured %.9f    predicted %.9f"%((acc/n/np.exp(-1j*(KT@(yp-y))*S/KP)).real*9e8,pred))
print("   -> THIS is the one divergence: int d^2z/|z|^2, i.e. the z-gluon collinear to the")
print("      parent.  In d = 2-2eps it is the 1/eps of the master z-integral.")
