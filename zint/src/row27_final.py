"""Row 27: the z integration done first -- asymptotic subtraction, closed form, p+ integral."""
import numpy as np, itertools
from scipy.linalg import expm
from row25_lib import val, KT, KP, s
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
X=(xp-yp)/s(xp-yp); B=s(x-w); kabs=np.sqrt(KT@KT); kh=KT/kabs; gE=np.euler_gamma
Kc =((x-w)@X)/(2*B)                    # K   = (x-w).X / 2B
Kk =((x-w)@kh)*(X@kh)/(2*B)            # K_k = [(x-w).khat][X.khat] / 2B
print("="*78); print("1.  THE SUBTRACTION FUNCTION"); print("="*78)
print("   I_asy(z) = e^{-ik.(y'-w)} e^{-i b k.y'} (e^{i b k.z}/|z|^2) (1/p+) x")
print("              { (p+^2/(k+^2 S)) (x-w)^m X^k' zhat^m zhat^k' /2B  +  (1/k+) K } ,  b = p+/k+")
print("   K = (x-w).X/2B = %.9f ,  K_k = [(x-w).khat][X.khat]/2B = %.9f"%(Kc,Kk))
def Iasy(Z,P):
    b=P/KP; zh=Z/np.sqrt((Z**2).sum(1))[:,None]; r2=(Z**2).sum(1)
    core=(1/P)*( (P**2/(KP**2*(P+KP)))*((zh@(x-w))*(zh@X))/(2*B) + Kc/KP )
    return np.exp(-1j*(KT@(yp-w)))*np.exp(-1j*b*(KT@yp))*np.exp(1j*b*(Z@KT))/r2*core
print()
print("   check it really is the tail of the full integrand:  <I - I_asy> x |z|^3")
print("   %-9s %16s %16s"%("|z|","<|I|> x |z|^2","<|I - I_asy|> x |z|^3"))
for R in (1e3,1e4,1e5):
    Z=ring(R,20000)
    print("   %-9.0e %16.8f %16.8f"%(R,abs(val(Z,0.37,phase=False).mean())*R*R,
          abs((val(Z,0.37)-Iasy(Z,0.37)).mean())*R**3))
print("   -> the remainder falls one power faster, so int d^2z [I - I_asy] converges ABSOLUTELY.")

print(); print("="*78); print("2.  THE CLOSED FORM OF int_{|z|>rho0} d^2z I_asy"); print("="*78)
print("   M_0 = 2 pi L ,   M_2^{mk'} = pi [ (L + 1/2) delta^{mk'} - khat^m khat^k' ] ,")
print("   with  L = log[ 2 e^{-gamma_E} k+ / (p+ |k| rho0) ]   (a2 = -1/2 exactly).")
def asy_closed(P,rho0):
    b=P/KP; L=np.log(2*np.exp(-gE)/(b*kabs*rho0))
    return (np.exp(-1j*(KT@(yp-w)))*np.exp(-1j*b*(KT@yp))
            *(np.pi/P)*( 2*Kc/KP*L + (P**2/(KP**2*(P+KP)))*((L+0.5)*Kc - Kk) ))
def asy_numeric(P,rho0,R,nr=1500,nth=9000):
    lo,hi=np.log(rho0),np.log(R); lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; tot=0j
    for r0 in np.exp(lr): tot+=Iasy(ring(r0,nth),P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot
print("   %-7s %-8s %22s %22s %9s"%("p+","rho0","numeric","closed form","rel.diff"))
for P,rho0 in ((0.37,1.0),(0.37,0.3),(0.05,1.0),(1.5,1.0)):
    R=200*KP/(P*kabs); nu=asy_numeric(P,rho0,R); cl=asy_closed(P,rho0)
    print("   %-7.2f %-8.2f %22s %22s %9.2e"%(P,rho0,"%.6f%+.6fj"%(nu.real,nu.imag),
          "%.6f%+.6fj"%(cl.real,cl.imag),abs(nu-cl)/abs(cl)),flush=True)

print(); print("="*78); print("3.  THE COLOUR AT LARGE |z| :  U(z) -> 1"); print("="*78)
Nc=3
lam=np.zeros((8,3,3),dtype=complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1; lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=lam/2; f=np.zeros((8,8,8))
for a,b,c in itertools.product(range(8),repeat=3): f[a,b,c]=(-2j*np.trace((T[a]@T[b]-T[b]@T[a])@T[c])).real
rng=np.random.default_rng(31)
def U(al):
    V=expm(1j*sum(al[a]*T[a] for a in range(8)))
    return np.array([[(2*np.trace(T[a]@V@T[b]@V.conj().T)).real for b in range(8)] for a in range(8)])
Uy,Ux,Uxp,Uw=[U(rng.normal(size=8)) for _ in range(4)]; Dif=Uy-Uxp
print("   setting U(z) = 1 in the (uniform-ordering) colour factor:")
print("     bra -> f^{c'ba} [U(y')-U(x')]^{c'e'} rho^{e'}(x')")
print("     ket -> f^{cbd} U^{ac}(w) rho^d(x)  -  f^{abd} U^{de}(x) rho^e(x)")
piece2=-np.einsum('xba,abd,de,xf->fe',f,f,Ux,Dif)
print("   the SECOND product collapses:  -sum_{a,b} f^{c'ba} f^{abd} = + N_c delta^{c'd}")
print("     residual to  N_c [ (U(y')-U(x'))^T U(x) ] :  %.2e"%np.abs(piece2-Nc*(Dif.T@Ux)).max())
piece1=np.einsum('xba,cbd,ac,xf->fd',f,f,Uw,Dif)
print("   the FIRST does not collapse (sum_b f^{c'ba}f^{cbd} is an irreducible 4-index object);")
print("     ||first piece|| = %.4f  -- keep it as it stands."%np.abs(piece1).max())

print(); print("="*78); print("4.  THE p+ INTEGRAL OF THE LOG PIECE, IN CLOSED FORM"); print("="*78)
print("   leading term:  (2 pi K/(p+ k+)) log(Xi/p+) ,   Xi = 2 e^{-gamma_E} k+ /(|k| rho0)")
print("   int_Lambda^{V-k+} dp+ (1/p+) log(Xi/p+) = (1/2)[ log^2(Xi/Lambda) - log^2(Xi/(V-k+)) ]")
import sympy as sp
p,Lm,Vv,Xi=sp.symbols('p Lambda vee Xi',positive=True)
prim=sp.integrate(sp.log(Xi/p)/p,p)
print("   primitive of log(Xi/p)/p :  %s   (= -(1/2) log^2(Xi/p))"%sp.simplify(prim))
print("   check d/dp of -(1/2)log^2(Xi/p) - log(Xi/p)/p = %s"%sp.simplify(sp.diff(-sp.log(Xi/p)**2/2,p)-sp.log(Xi/p)/p))
print()
print("   subleading term:  (pi p+/(k+^2 S)) [ (log(Xi/p+) + 1/2) K - K_k ] -- finite at both ends;")
k_=sp.Symbol('kplus',positive=True); Kk_,Kc_=sp.symbols('K_k K')
sub=sp.integrate(p/(p+k_)*(sp.log(Xi/p)+sp.Rational(1,2)),p)
print("   int dp p/(p+k+) (log(Xi/p)+1/2) = %s"%sp.simplify(sub))
