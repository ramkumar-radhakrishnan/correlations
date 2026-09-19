"""Row 17: divergences of N_2, and whether the log carries P_gg."""
import numpy as np, itertools
K=1.0; kt=np.array([0.7,-0.4]); s=lambda v:v@v
def Wji(x,z,w,P):
    S=P+K; xz,xw,zw=x-z,x-w,z-w
    return (np.eye(2)*(s(xz)-s(xw))/(2*s(zw))
            + (S/K)*(np.outer(xz,zw)/s(zw) - np.outer(xz,xw)/(2*s(xw)))
            + (S/P)*(np.outer(zw,xw)/s(zw) + np.outer(xz,xw)/(2*s(xz))))
Dd=lambda x,z,w,P: P*s(x-z)+K*s(x-w)
def F(x,xp,z,w,wp,P=0.37):
    S=P+K
    return (np.exp(-1j*(kt@(wp-w)))*np.einsum('ji,ji->',Wji(x,z,w,P),Wji(xp,z,wp,P))
            *(P*K/S**2)/(Dd(x,z,w,P)*Dd(xp,z,wp,P)))
names=['x',"x'",'z','w',"w'"]
base=[np.array([0.31,-0.77]),np.array([-0.52,0.19]),np.array([0.05,-0.23]),
      np.array([1.13,0.42]),np.array([-0.31,0.88])]
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
print("   (this includes x,z,w -> point, where the energy denominator p+(x-z)^2+k+(x-w)^2")
print("    vanishes like rho^2: integrand rho^-2 against a three-point measure rho^3 drho.)")

print(); print("="*78); print("2.  LARGE-DISTANCE END"); print("="*78)
print("  %-5s %13s %13s %13s   %s"%("var","R=1e2","R=1e3","R=1e4","verdict"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e3,1e4):
        nth=int(min(max(4096,60*np.sqrt(kt@kt)*R),300_000)); acc=0j
        for th in np.linspace(0,2*np.pi,nth,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); acc+=F(*pt)
        v.append(abs(acc/nth)*R**2)
    flat=v[2]>1e-9 and v[2]>0.2*v[0]
    print("  %-5s %13.5e %13.5e %13.5e   %s"%(nm,*v,"*** LOG DIV ***" if flat else "convergent"))

print(); print("="*78); print("3.  THE LARGE-|z| COEFFICIENT, AND ITS zeta DEPENDENCE"); print("="*78)
x,xp,w,wp=base[0],base[1],base[3],base[4]
Kv=lambda a,b:(a-b)/s(a-b)
kk=Kv(x,w)@Kv(xp,wp); ph=np.exp(-1j*(kt@(wp-w)))
def plateau(P,R=3e4,nth=6000):
    acc=0j
    for t in np.linspace(0,2*np.pi,nth,endpoint=False):
        acc+=F(x,xp,R*u(t),w,wp,P)
    return acc/nth*R**2
Cuv=lambda zt: zt/(1-zt)+(1-zt)/zt+zt*(1-zt)
print("   mechanism: W_{ji}/[S D] -> -(1/(2 p+ k+)) (x-w)^i (x-z)^j/[..], so the product of")
print("   the two vertices gives a constant x 1/|z|^2 with no angular cancellation.")
print()
print("   %-7s %-9s %16s %16s %13s %11s"%("p+","zeta","plateau/phase","1/(4 p+ k+) x K","C_UV(zeta)","ratio"))
vals=[]
for Pp in (0.15,0.37,0.70,1.50,3.00,6.00):
    zt=K/(K+Pp); v=(plateau(Pp)/ph).real; pred=kk/(4*Pp*K); c=Cuv(zt)
    vals.append((zt,v,c)); print("   %-7.2f %-9.5f %16.9f %16.9f %13.6f %11.6f"%(Pp,zt,v,pred,c,v/c))
r=[v/c for _,v,c in vals]
print()
print("   ratio to C_UV spans %.5f .. %.5f  -> NOT constant: the log does not carry P_gg."%(min(r),max(r)))
print("   what IS flat is the ratio to the SOFT piece zeta/zetabar = k+/p+ alone:")
print("   %-9s %16s"%("zeta","plateau x p+"))
for (zt,v,c),Pp in zip(vals,(0.15,0.37,0.70,1.50,3.00,6.00)):
    print("   %-9.5f %16.9f"%(zt,v*Pp))
print("   constant  =>  plateau = (1/(4 k+ p+)) K  ~  zeta/zetabar , the zeta -> 1 (soft)")
print("   end of P_gg only, with no zetabar/zeta and no zeta zetabar partner.")
