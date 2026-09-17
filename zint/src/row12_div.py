"""Row 12: locate the transverse divergence of the 3-rho piece, before and after p+."""
import numpy as np, itertools
K=1.0; Pp=0.37; S=Pp+K; Lam=1e-8; Pv=1e7; kt=np.array([0.7,-0.4])
s=lambda v:v@v
def pieces(x,y,z,w,xp,wp):
    A=s(x-z); B=s(x-w); C=s(z-w)
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)          # P^{ij}
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)          # K^{ij}
    return A,B,C,Pm,Km
def F_un(x,y,z,w,xp,wp,P=Pp):
    A,B,C,Pm,Km=pieces(x,y,z,w,xp,wp); Ss=P+K
    br = np.eye(2)*(A-B)/(2*Ss*C) + Pm/P + Km/K
    Ki=(xp-wp)/s(xp-wp); Kj=(y-z)/s(y-z)
    return np.exp(-1j*(kt@(wp-w)))*np.einsum('ij,i,j->',br,Ki,Kj)/(P*A+K*B)
def J123(A,B):
    J1=np.log(Pv*K*B/(Lam*(Pv*A+K*B)))/(K*B)
    J2=np.log((Pv*A+K*B)/(Lam*A+K*B))/A
    J3=np.log((Pv+K)*B/(Pv*A+K*B))/(K*(B-A))
    return J1,J2,J3
def F_int(x,y,z,w,xp,wp):
    A,B,C,Pm,Km=pieces(x,y,z,w,xp,wp); J1,J2,J3=J123(A,B)
    br = np.eye(2)*(A-B)/(2*C)*J3 + Pm*J1 + Km/K*J2
    Ki=(xp-wp)/s(xp-wp); Kj=(y-z)/s(y-z)
    return np.exp(-1j*(kt@(wp-w)))*np.einsum('ij,i,j->',br,Ki,Kj)
names=['x','y','z','w',"x'","w'"]
base=[np.array([0.31,-0.77]),np.array([0.62,0.41]),np.array([0.05,-0.23]),
      np.array([1.13,0.42]),np.array([-0.52,0.19]),np.array([-0.31,0.88])]
u=lambda t:np.array([np.cos(t),np.sin(t)]); OFF=[0,2.1,4.0,1.05,5.2,3.3]

for tag,F in (("BEFORE the p+ integral",F_un),("AFTER the p+ integral",F_int)):
    print("="*78); print("COLLAPSE SCAN, %s"%tag); print("="*78)
    found=[]
    for n in range(2,7):
        for sub in itertools.combinations(range(6),n):
            v=[]
            for rho in (1e-2,1e-4,1e-6):
                acc=0.0
                for th in np.linspace(0,2*np.pi,24,endpoint=False):
                    pt=[b.copy() for b in base]; c=base[sub[0]]
                    for a,i in enumerate(sub): pt[i]=c+rho*u(th+OFF[a])
                    acc+=abs(F(*pt))
                v.append(acc/24*rho**(2*(n-1)))
            if v[1]>0.3*v[0] and v[2]>0.3*v[1]:
                found.append(sub)
                print("   {%s}: %.4e %.4e %.4e  *** DIVERGENT ***"
                      %(",".join(names[i] for i in sub),*v))
    print("   ->", "NO UV in any of the 57 regions" if not found else "%d divergent"%len(found))
    print()

print("="*78); print("LARGE-DISTANCE, one variable at a time (bare, phase set to 1)"); print("="*78)
print("  %-5s %14s %14s %14s   %s"%("var","R=1e2","R=1e4","R=1e6","verdict"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e4,1e6):
        acc=0j
        for th in np.linspace(0,2*np.pi,2048,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th)
            A,B,C,Pm,Km=pieces(*pt); br=np.eye(2)*(A-B)/(2*S*C)+Pm/Pp+Km/K
            Ki=(pt[4]-pt[5])/s(pt[4]-pt[5]); Kj=(pt[1]-pt[2])/s(pt[1]-pt[2])
            acc+=np.einsum('ij,i,j->',br,Ki,Kj)/(Pp*A+K*B)
        v.append(abs(acc/2048)*R**2)
    flat = v[2]>1e-9 and v[2]>0.2*v[0]
    print("  %-5s %14.6e %14.6e %14.6e   %s"%(nm,*v,"*** LOG DIVERGENT ***" if flat else "convergent"))
print()
print("="*78); print("JOINT LARGE-DISTANCE (two variables to infinity together)"); print("="*78)
print("  %-10s %14s %14s %14s   %s"%("pair","R=1e2","R=1e4","R=1e6","verdict"))
for i,j in itertools.combinations(range(6),2):
    v=[]
    for R in (1e2,1e4,1e6):
        acc=0j
        for th in np.linspace(0,2*np.pi,512,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); pt[j]=R*u(th+1.7)
            A,B,C,Pm,Km=pieces(*pt); br=np.eye(2)*(A-B)/(2*S*C)+Pm/Pp+Km/K
            Ki=(pt[4]-pt[5])/s(pt[4]-pt[5]); Kj=(pt[1]-pt[2])/s(pt[1]-pt[2])
            acc+=np.einsum('ij,i,j->',br,Ki,Kj)/(Pp*A+K*B)
        v.append(abs(acc/512)*R**4)      # two points -> measure R^4 dR/R
    flat = v[2]>1e-9 and v[2]>0.2*v[0]
    if flat:
        print("  %-10s %14.6e %14.6e %14.6e   *** DIVERGENT ***"%("{%s,%s}"%(names[i],names[j]),*v))
print("  (only pairs flagged are printed)")
