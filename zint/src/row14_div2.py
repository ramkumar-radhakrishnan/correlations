"""Row 14: UV collapse scan + large-distance + the + prescription test (fast)."""
import numpy as np, itertools
K=1.0; kt=np.array([0.7,-0.4]); s=lambda v:v@v; kn=np.sqrt(kt@kt)
def Tt(x,z,w,P):
    Ss=K+P; A=s(x-z); B=s(x-w); C=s(z-w)
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)
    M=np.eye(2)*(A-B)/(2*Ss*C) + Km/K + Pm/P
    return (P/Ss)*np.eye(2)*np.trace(M) - M - (P/K)*M.T, A, B
def F(x,xp,yp,z,w,P=0.37):
    T,A,B=Tt(x,z,w,P); kap=(kt@(yp-z))/K
    ph=np.exp(-1j*(kt@(yp-w)))*np.exp(-1j*kap*P)
    return ph*np.einsum('km,k,m->',T,(xp-yp)/s(xp-yp),(yp-z)/s(yp-z))/(P*A+K*B)
names=['x',"x'","y'",'z','w']
base=[np.array([0.31,-0.77]),np.array([-0.52,0.19]),np.array([0.62,0.41]),
      np.array([0.05,-0.23]),np.array([1.13,0.42])]
u=lambda t:np.array([np.cos(t),np.sin(t)]); OFF=[0,2.1,4.0,1.05,5.2]
print("="*78); print("1.  UV: ALL 26 COLLAPSING SUBSETS"); print("="*78)
bad=[]
for n in range(2,6):
    for sub in itertools.combinations(range(5),n):
        v=[]
        for rho in (1e-2,1e-4,1e-6):
            acc=0.0
            for th in np.linspace(0,2*np.pi,28,endpoint=False):
                pt=[b.copy() for b in base]; c=base[sub[0]]
                for a,i in enumerate(sub): pt[i]=c+rho*u(th+OFF[a])
                acc+=abs(F(*pt))
            v.append(acc/28*rho**(2*(n-1)))
        if v[1]>0.3*v[0] and v[2]>0.3*v[1]:
            bad.append(sub); print("   {%s}: %.3e %.3e %.3e *** DIV ***"%(",".join(names[i] for i in sub),*v))
print("   ->","NO UV in any of the 26 regions" if not bad else "%d divergent"%len(bad))

print(); print("="*78); print("2.  LARGE-DISTANCE END (R up to 1e4)"); print("="*78)
print("  %-5s %14s %14s %14s   %s"%("var","R=1e2","R=1e3","R=1e4","verdict"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e3,1e4):
        nth=int(min(max(4096,60*kn*R),400_000)); acc=0j
        for th in np.linspace(0,2*np.pi,nth,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); acc+=F(*pt)
        v.append(abs(acc/nth)*R**2)
    flat=v[2]>1e-8 and v[2]>0.2*v[0]
    print("  %-5s %14.5e %14.5e %14.5e   %s"%(nm,*v,"*** LOG DIV ***" if flat else "convergent"))

print(); print("="*78); print("3.  THE + PRESCRIPTION TEST"); print("="*78)
x,xp,yp,w=base[0],base[1],base[2],base[4]
def core(z,P,phase=True):
    A=s(x-z); B=s(x-w); C=s(z-w)
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)
    val=-np.einsum('km,k,m->',Pm,(xp-yp)/s(xp-yp),(yp-z)/s(yp-z))
    if not phase: return val/(K*B)                       # g(0): no phase, D at p+=0
    return val*np.exp(-1j*((kt@(yp-z))/K)*P)/(P*A+K*B)
print("   the 1/p+ structure carries  g(p+) = e^{-i kappa p+}/(p+A + k+B),  kappa = k.(y'-z)/k+")
print("   the + prescription replaces it by g(0) = 1/(k+B), which DELETES the z-phase.")
print()
print("   %-9s %20s %20s %20s"%("|z|","subtraction g(0)","full, p+=0.37","full, p+=0.01"))
for R in (1e2,1e3,1e4,1e5):
    nth=int(min(max(8192,60*kn*R*0.37),600_000)); th=np.linspace(0,2*np.pi,nth,endpoint=False)
    a=b=c=0j
    for t in th:
        zz=R*u(t); a+=core(zz,0,False); b+=core(zz,0.37); c+=core(zz,0.01)
    print("   %-9.0e %20.9f %20.4e %20.4e"%(R,(a/nth*R**2).real,abs(b/nth)*R**2,abs(c/nth)*R**2))
pred=-0.5*((xp-yp)@(x-w))/(s(xp-yp)*K*s(x-w))
print()
print("   predicted plateau of the subtraction term = -(1/2)(x'-y').(x-w)/[(x'-y')^2 k+ (x-w)^2]")
print("                                             = %.9f"%pred)
