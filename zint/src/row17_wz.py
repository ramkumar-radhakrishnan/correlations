"""Row 17: does the identity f f U(w) U(z) = N_c U act at w -> z, and can it make a divergence?"""
import numpy as np, itertools
from scipy.linalg import expm
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3):
    f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(2024)
def Uadj(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])

print("="*78); print("1.  THE IDENTITY IS TRUE"); print("="*78)
V=Uadj(rng.normal(size=8)*1.3)
print("   U orthogonal:      max|U^T U - 1|                      = %.2e"%np.abs(V.T@V-np.eye(8)).max())
print("   f^{abc}U^{cd}(z) = U^{aa'}(z)U^{bb'}(z) f^{a'b'd}      = %.2e"
      %np.abs(np.einsum('abc,cd->abd',f,V)-np.einsum('ax,by,xyd->abd',V,V,f)).max())
print("   f^{acd}f^{bcd}   = N_c delta^{ab}                      = %.2e"
      %np.abs(np.einsum('acd,bcd->ab',f,f)-Nc*np.eye(8)).max())
print("   Both are needed.  The COLLAPSE  f f U U -> N_c U  requires two Wilson lines")
print("   AT THE SAME POINT contracted through a COMMON SUMMED INDEX of the two f's,")
print("   so that sum_a U^{a c'} U^{a c} = delta^{c' c} can fire first.")

print(); print("="*78); print("2.  WHERE IT ACTUALLY FIRES IN YOUR ROW"); print("="*78)
Ux,Uxp,Uz=Uadj(rng.normal(size=8)),Uadj(rng.normal(size=8)),V
Uw,Uwp=Uadj(rng.normal(size=8)),Uadj(rng.normal(size=8))
C22=np.einsum('axy,azw,bx,fy,bz,we->fe',f,f,Uz,Uxp,Uz,Ux)
print("   (a)  (2)x(2):   f^{ac'd'}U^{bc'}(z) x f^{acd}U^{bc}(z)")
print("        both lines at z, common summed index b  ->  delta^{c'c}  ->  N_c delta^{d'd}")
print("        C22 = N_c [U(x')U(x)]^{e'e} :  max diff = %.2e   *** EXACT, for ANY w ***"%np.abs(C22-Nc*(Uxp@Ux)).max())
C11w=np.einsum('xbp,ax,ybd,ay->pd',f,Uw,f,Uw)
C11g=np.einsum('xbp,ax,ybd,ay->pd',f,Uwp,f,Uw)
print("   (b)  (1)x(1) at w'=w:  f^{c'bd'}U^{ac'}(w') x f^{cbd}U^{ac}(w)")
print("        common summed index a  ->  delta^{c'c}  ->  f^{c'bd'}f^{c'bd}=N_c delta^{d'd}")
print("        C11(w'=w) = N_c delta^{d'd} :  max diff = %.2e   (||C11|| generic = %.4f)"
      %(np.abs(C11w-Nc*np.eye(8)).max(),np.abs(C11g).max()))

print(); print("="*78); print("3.  IT DOES **NOT** FIRE AT w -> z"); print("="*78)
print("   the cross terms are   f^{ac'd'}U^{bc'}(z)  x  f^{cbd}U^{ac}(w) .")
print("   U(w) carries the OBSERVED gluon's colour a ; U(z) carries the UNOBSERVED")
print("   gluon's colour b .  Those are two DIFFERENT free indices, summed independently,")
print("   one against each f.  Setting w=z makes the two matrices equal but leaves the")
print("   index structure untouched: no sum_a U^{ac'}U^{ac}, hence no delta, hence no N_c.")
Sig=lambda Uw: np.einsum('axy,bx,cbd,ac->yd',f,Uz,f,Uw)
gen=np.array([np.abs(Sig(Uadj(rng.normal(size=8)))).max() for _ in range(500)])
print()
print("   ||Sigma(w=z)||_max          = %.6f"%np.abs(Sig(Uz)).max())
print("   ||Sigma(w)||  over 500 random w :  min %.4f  mean %.4f  max %.4f"%(gen.min(),gen.mean(),gen.max()))
for n,c in (('N_c d^{d\'d}',Nc*np.eye(8)),('-N_c d',-Nc*np.eye(8)),('N_c U(z)',Nc*Uz),('-N_c U(z)',-Nc*Uz),
            ('(N_c/2)U(z)',Nc/2*Uz),('-(N_c/2)U(z)',-Nc/2*Uz)):
    print("   Sigma(w=z) - [%-14s] : %.3e"%(n,np.abs(Sig(Uz)-c).max()))
print("   -> no candidate N_c form fits, and Sigma(w=z) sits INSIDE the generic range.")
print("      The w -> z point is colour-wise completely unremarkable.")

print(); print("="*78); print("4.  AND EVEN IF IT DID: THE TRANSVERSE INTEGRAND FORBIDS A LOG"); print("="*78)
s=lambda v:v@v
def vert(x,z,w,P,K):
    S=P+K; A,B,C=s(x-z),s(x-w),s(z-w)
    Vw=(z-w)/C-(x-w)/(2*B); Vz=(z-w)/C+(x-z)/(2*A)
    return (np.eye(2)*(A-B)/(2*C)+(S/K)*np.outer(x-z,Vw)+(S/P)*np.outer(Vz,x-w)), P*A+K*B
xs,xps,wp0,z0=[rng.normal(size=2) for _ in range(4)]
kt=np.array([0.6,-0.35]); P,K=0.37,1.0
base=rng.normal(size=8)*0.8
Uof=lambda pt: Uadj(base*np.exp(-s(pt)/6.0))
rhox,rhoxp=rng.normal(size=8),rng.normal(size=8)
def colour(z,w,wp):
    Uw,Uwp,Uzz=Uof(w),Uof(wp),Uof(z)
    ket=np.einsum('cbd,ac,d->ab',f,Uw,rhox)-np.einsum('acd,bc,de,e->ab',f,Uzz,Uof(xs),rhox)
    bra=np.einsum('cbd,ac,d->ab',f,Uwp,rhoxp)-np.einsum('acd,bc,ed,e->ab',f,Uzz,Uof(xps),rhoxp)
    return np.einsum('ab,ab->',bra,ket)
def full(z,w,wp):
    Wa,D=vert(xs,z,w,P,K); Wb,Dp=vert(xps,z,wp,P,K)
    return np.exp(-1j*(kt@(wp-w)))*np.einsum('ji,ji->',Wa,Wb)*(P*K/(P+K)**2)/(D*Dp)*colour(z,w,wp)
print("   full integrand INCLUDING Wilson lines and colour, as rho = |w-z| -> 0 :")
print("   %-9s %15s %16s %18s"%("rho","<|full|>","<|full|> x rho","|<full>| signed"))
for rho in (1e-1,1e-2,1e-3,1e-4,1e-5,1e-6):
    A_=0.0; Sg=0j; n=720
    for t in np.linspace(0,2*np.pi,n,endpoint=False):
        v=full(z0,z0+rho*np.array([np.cos(t),np.sin(t)]),wp0); A_+=abs(v); Sg+=v
    print("   %-9.0e %15.6e %16.9f %18.9f"%(rho,A_/n,A_/n*rho,abs(Sg/n)))
print("   <|full|> x rho -> const : integrand ~ 1/rho against d^2rho ~ rho drho => int drho.")
print("   CONVERGENT.  The signed average tends to a CONSTANT, not 1/rho, because every")
print("   1/rho piece of W^{ji} is odd in (z-w) and averages away: softer still.")
print("   A log needs rho^-2.  Colour is a bounded factor (contractions of ORTHOGONAL")
print("   matrices): it can only make a region more convergent, never less.")

print(); print("="*78); print("5.  THE REGION WHERE THE IDENTITY *DOES* FIRE:  w' -> w"); print("="*78)
print("   %-9s %15s %16s %18s"%("rho","<|full|>","<|full|> x rho","|<full>| signed"))
for rho in (1e-1,1e-2,1e-3,1e-4,1e-5):
    A_=0.0; Sg=0j; n=720
    for t in np.linspace(0,2*np.pi,n,endpoint=False):
        v=full(z0,wp0+rho*np.array([np.cos(t),np.sin(t)]),wp0); A_+=abs(v); Sg+=v
    print("   %-9.0e %15.6e %16.9f %18.9f"%(rho,A_/n,A_/n*rho,abs(Sg/n)))
print("   <|full|> itself tends to a CONSTANT: the integrand is REGULAR at w'=w.")
print("   So the N_c collapse there multiplies a perfectly convergent region.  (It is")
print("   also the k_perp-integrated corner, suppressed by the phase in any case.)")
