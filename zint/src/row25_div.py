"""Row 25: colour, and the transverse divergences of the A^dag Cbar^dag x B_2 row."""
import numpy as np, itertools
from scipy.linalg import expm
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3): f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(19)
def Uadj(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
print("="*80); print("1.  COLOUR: DOES ANYTHING COLLAPSE TO N_c HERE?"); print("="*80)
Uy,Ux,Uxp,Uz,Uw=[Uadj(rng.normal(size=8)) for _ in range(5)]
print("   bra (A^dag Cbar^dag) : free indices a on f, b on U(z)")
print("   ket (B_2)            : free indices a on U(w)/f , b on f/U(z)")
print("   The four products, with index orders made uniform:")
lab=["bra1 x ket1","bra1 x ket2","bra2 x ket1","bra2 x ket2"]
E=[np.einsum('xya,yb,ex,bcd,ac->ed',f,Uz,Uy,f,Uw),
   np.einsum('xya,yb,ex,acd,bc,dg->eg',f,Uz,Uy,f,Uz,Ux),
   np.einsum('xya,yb,ex,bcd,ac->ed',f,Uz,Uxp,f,Uw),
   np.einsum('xya,yb,ex,acd,bc,dg->eg',f,Uz,Uxp,f,Uz,Ux)]
for n,M in zip(lab,E):
    print("     %-12s : ||M||_max = %8.4f ,  ||M - N_c x (anything proportional to 1)|| tested:"%(n,np.abs(M).max()))
    best=min([("N_c I",np.abs(M-Nc*np.eye(8)).max()),("-N_c I",np.abs(M+Nc*np.eye(8)).max()),
              ("N_c U(x)",np.abs(M-Nc*Ux).max()),("N_c U(z)",np.abs(M-Nc*Uz).max())],key=lambda t:t[1])
    print("        closest simple N_c form: %-10s  residual %.3f"%best)
print()
print("   sum_b U^{d'b}(z) U^{bc}(z) = [U(z)U(z)]^{d'c}  (NOT a delta) :  ||.-1|| = %.3f"%np.abs(Uz@Uz-np.eye(8)).max())
print("   sum_b U^{bd'}(z) U^{bc}(z) = delta^{d'c}                     :  ||.-1|| = %.1e"%np.abs(Uz.T@Uz-np.eye(8)).max())
print("   => as in the previous rows, your bra writes U^{d'b}(z) in term 1 and U^{bd'}(z) in term 2.")
print("      Make them uniform.  But even then NO overall N_c appears here, because on the bra side")
print("      a sits on f and b on U(z), while on the ket side a sits on U(w) (term 1) -- the two free")
print("      indices are never contracted through a COMMON summed index of two f's.")
print("      Only bra2 x ket2 has both U(z)'s: there sum_b gives delta and f^{c'd'a}f^{acd} -> N_c delta^{c'd},")
print("      residual to N_c U(x')^T U(x) : %.2e"
      %np.abs(np.einsum('xya,yb,xe,acd,bc,dg->eg',f,Uz,Uxp,f,Uz,Ux)-Nc*(Uxp.T@Ux)).max())

print(); print("="*80); print("2.  THE INTEGRAND AFTER THE w' INTEGRATION"); print("="*80)
s=lambda v:v@v; dd=np.eye(2); KT=np.array([0.6,-0.35]); KP=1.0
def bracket(x,z,w,P,K,dperp=2.0):
    S=P+K; A,B,C=s(x-z),s(x-w),s(z-w); xz,xw,zw=x-z,x-w,z-w
    R =dd*dperp*P/(2*C*S**2)*(A-B)
    R+=dd*(P/(K*S))*((xz@zw)/C-(xz@xw)/(2*B)-(A-B)/(2*C))
    R+=dd*(1/S)     *((xw@zw)/C+(xz@xw)/(2*A)-(A-B)/(2*C))
    R+=-(1/P)   *( np.outer(zw,xw)/C + np.outer(xz,xw)/(2*A) )          # CORRECTED term 4
    R+=-(P/K**2)*( np.outer(zw,xz)/C - np.outer(xw,xz)/(2*B) )
    R+=-(1/K)   *( np.outer(xz,zw)/C - np.outer(xz,xw)/(2*B)
                  +np.outer(xw,zw)/C + np.outer(xw,xz)/(2*A) )
    return R
def integ(x,xp,yp,z,w,P=0.37,K=KP,phase=True):
    N=(yp-z)/s(yp-z); X=(xp-yp)/s(xp-yp); D=P*s(x-z)+K*s(x-w)
    val=np.einsum('m,k,mk->',N,X,bracket(x,z,w,P,K))/D
    if phase: val=val*np.exp(-1j*(KT@(yp-w)))*np.exp(-1j*(P/K)*(KT@(yp-z)))
    return val
names=['x',"x'","y'",'z','w']
base=[np.array([0.41,-0.63]),np.array([-0.58,0.22]),np.array([0.09,0.51]),
      np.array([1.07,-0.34]),np.array([-0.25,0.86])]
u=lambda t:np.array([np.cos(t),np.sin(t)]); OFF=[0,2.1,4.0,1.05,5.2]
print("   variables: x, x', y', z, w   (w' integrated out, y never appears)")
print()
print("="*80); print("3.  SHORT DISTANCE: ALL 26 COLLAPSING SUBSETS"); print("="*80)
print("   %-18s %11s %11s %11s   %s"%("subset","rho=1e-2","1e-4","1e-6","verdict"))
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
            print("   %-18s %11.3e %11.3e %11.3e   %s"%("{"+",".join(names[i] for i in sub)+"}",*v,
                  "*** DIV ***" if dv else "convergent"))
print("   (4- and 5-point subsets all convergent)")
print("   ==> %s"%("NO short-distance divergence" if not bad else "%d divergent: %s"%(len(bad),bad)))
print()
print("="*80); print("4.  LARGE DISTANCE"); print("="*80)
print("   %-6s %13s %13s %13s   %s"%("var","R=1e2","R=1e3","R=1e4","verdict"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e3,1e4):
        nth=int(min(max(4096,80*np.sqrt(KT@KT)*R),200_000)); acc=0j
        for th in np.linspace(0,2*np.pi,nth,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); acc+=integ(*pt)
        v.append(abs(acc/nth)*R*R)
    print("   %-6s %13.5e %13.5e %13.5e   %s"%(nm,*v,
          "*** LOG DIV ***" if (v[2]>1e-9 and v[2]>0.25*v[0]) else "convergent"))
print()
print("   and with the phase switched off (to expose what oscillation is hiding):")
print("   %-6s %13s %13s %13s"%("var","R=1e2","R=1e3","R=1e4"))
for i,nm in enumerate(names):
    v=[]
    for R in (1e2,1e3,1e4):
        acc=0.0
        for th in np.linspace(0,2*np.pi,4096,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(th); acc+=integ(*pt,phase=False)
        v.append(abs(acc/4096)*R*R)
    print("   %-6s %13.5e %13.5e %13.5e"%(nm,*v))
