"""Row 13: check the quoted p+-integrated 3-rho result against direct integration."""
import numpy as np
from scipy.integrate import quad
K=1.0; Lam=1e-7; V=1e6; Pv=V-K
s=lambda v:v@v
rng=np.random.default_rng(4)

def setup(x,y,z,w,xp,wp):
    A=s(x-z); B=s(x-w); C=s(z-w)
    Ki=(xp-wp)/s(xp-wp); Kj=(y-z)/s(y-z)
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)     # P^{ij}  (the 1/p+ bracket)
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)     # K^{ij}  (the 1/k+ bracket)
    return A,B,C,Ki,Kj,Pm,Km

def direct(x,y,z,w,xp,wp):
    """(1/k+) x int_Lambda^{V-k+} dp+ of the contracted pre-integration integrand,
       i.e. exactly the brace the quoted result multiplies by -i g^4/(16 pi^5 k+)."""
    A,B,C,Ki,Kj,Pm,Km=setup(x,y,z,w,xp,wp)
    dd=np.einsum('ij,i,j->',np.eye(2),Ki,Kj); pp=np.einsum('ij,i,j->',Pm,Ki,Kj)
    kk=np.einsum('ij,i,j->',Km,Ki,Kj)
    f=lambda p: (dd*(A-B)/(2*(p+K)*C) + pp/p + kk/K)/(p*A+K*B)
    return quad(lambda t:f(np.exp(t))*np.exp(t), np.log(Lam), np.log(Pv), limit=900)[0]

def quoted(x,y,z,w,xp,wp):
    """the quoted three groups, divided by (-i g^4/(16 pi^5)) -- so it should equal 'direct'."""
    A,B,C,Ki,Kj,Pm,Km=setup(x,y,z,w,xp,wp)
    L1=np.log((Pv*A+K*B)/(Lam*A+K*B))
    L2=np.log(V/(Lam+K))
    L3=np.log((V-K)/Lam)
    E1 = (np.einsum('ij,i,j->',np.eye(2),Ki,Kj)/(2*C)
          + np.einsum('ij,i,j->',Km,Ki,Kj)/A
          - np.einsum('ij,i,j->',Pm,Ki,Kj)/B)
    E2 = np.einsum('ij,i,j->',np.eye(2),Ki,Kj)/C
    E3 = np.einsum('ij,i,j->',Pm,Ki,Kj)/B
    return (E1*L1 - 0.5*E2*L2 + E3*L3)/K

print("="*76); print("DIRECT p+ INTEGRAL  vs  THE QUOTED THREE GROUPS"); print("="*76)
print("   (both divided by the common -i g^4/(16 pi^5); Lambda=1e-7, V=1e6)")
print("   %-6s %18s %18s %12s"%("config","direct","quoted","rel.diff"))
for t in range(6):
    pts=[rng.normal(size=2)*1.2 for _ in range(6)]
    a,b=direct(*pts),quoted(*pts)
    print("   %-6d %18.9f %18.9f %12.2e"%(t,a,b,abs(a-b)/max(abs(a),1e-30)))

print(); print("="*76); print("WHY IT WORKS -- the exact decomposition"); print("="*76)
print("   J3 = int dp/[(p+k+)D] = (1/(k+(B-A))) [ log((p+k+)/(pA+k+B)) ]_Lambda^{V-k+}")
print("                         = (1/(k+(B-A))) log[ V(Lambda A+k+B) / ((Lambda+k+)((V-k+)A+k+B)) ]")
print("   so   (A-B)/(2C) J3 = -(1/(2 k+ C)) [ log(V/(Lambda+k+)) - log(((V-k+)A+k+B)/(Lambda A+k+B)) ]")
print("        -> your delta^ij/(2C) in group 1  AND  your group 2.   Exact endpoints kept. ")
print()
print("   J1 = int dp/[p D] = (1/(k+B)) log[ (V-k+) k+B / (Lambda((V-k+)A+k+B)) ]")
print("                     = (1/(k+B)) [ log((V-k+)/Lambda) - log(((V-k+)A+k+B)/(Lambda A+k+B)) ]")
print("                       + O(Lambda)")
print("        -> your group 3 (the blue one) AND the -P^{ij}/B piece of group 1.")
print()
print("   J2 = int dp/D = (1/A) log[((V-k+)A+k+B)/(Lambda A+k+B)]")
print("        -> the +K^{ij}/A piece of group 1.   All five bracket terms accounted for.")
E=lambda q:q
print()
print("="*76); print("BOOKKEEPING OF THE LOGS"); print("="*76)
print("   rapidity log(1/Lambda) : ONLY in the blue group-3 term, residue P^{ij}/(k+ B).")
print("                            -> apply [1/p+]_+ to that term alone.")
print("   log V : appears in group 1 (via L1 ~ log(V A/(k+B))) and group 2 (via L2 ~ log(V/k+)).")
print("           For the delta^ij structure the two CANCEL:")
print("              log(V/(Lambda+k+)) - log(((V-k+)A+k+B)/(Lambda A+k+B)) -> log(k+B/(k+ A)) = log(B/A)")
print("           leaving  -delta^ij log[(x-w)^2/(x-z)^2] / (2 k+ (z-w)^2) , free of V and Lambda.")
print("           For K^{ij}/A the log V SURVIVES: residue 1/(x-z)^2 -- a genuine large-p+ log.")
for t in range(2):
    pts=[rng.normal(size=2)*1.2 for _ in range(6)]
    A,B,C,Ki,Kj,Pm,Km=setup(*pts)
    combo=np.log(V/(Lam+K))-np.log((Pv*A+K*B)/(Lam*A+K*B))
    print("      check: log(V/(Lam+k+)) - L1 = %12.8f   vs  log(B/A) = %12.8f"%(combo,np.log(B/A)))

print()
print("="*76); print("DOES THE LARGE-|z| DIVERGENCE SURVIVE THE p+ INTEGRATION?"); print("="*76)
x=np.array([0.31,-0.77]); y=np.array([0.62,0.41]); w=np.array([1.13,0.42])
xp=np.array([-0.52,0.19]); wp=np.array([-0.31,0.88])
print("   %-10s %20s"%("|z|","angular avg x |z|^2"))
for R in (1e2,1e3,1e4,1e5):
    acc=0.0
    for th in np.linspace(0,2*np.pi,1024,endpoint=False):
        acc+=quoted(x,y,R*np.array([np.cos(th),np.sin(th)]),w,xp,wp)
    print("   %-10.0e %20.8f"%(R,acc/1024*R**2))
print("   -> yes: flat, so the |z| -> infinity logarithm is untouched by the p+ integral.")
print("      It still has to be regulated (d = 2-2eps) and combined with the 4-rho term.")
