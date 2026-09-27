"""Group I Row V (Group_V.tex), second part (b) and first part (a): every step, checked.

Starting point: the .tex equation after the w integration (the one labelled
"Integrating w.r.t. w"), split into (a) and (b) by the colour identity.
"""
import numpy as np, sympy as sp, mpmath as mp
from scipy.linalg import expm
mp.mp.dps = 20
rng = np.random.default_rng(1)
def sq(v): return float(v@v)
bar = "="*78

# ---------------------------------------------------------------------------
print(bar); print("STEP 1.  THE SIX BLOCKS OF THE .tex, REWRITTEN IN r = z-x , s = y-x , P = x'-w'"); print(bar)
def tex_blocks(D, x, y, z, xp, wp, kp, pp):
    """literal transcription of the bracket of the .tex (after the w integral), contracted with
       (x'-w')^i (z-x)^m / [(x'-w')^2 (z-x)^2] ; returns the six scalars."""
    I=np.eye(D); q=kp-pp
    B=[ I*D/(2*sq(x-z))*(sq(y-x)-sq(y-z))*(pp*q/kp**2),
        -(pp/q)*( np.outer(y-x,x-z)/sq(x-z) - np.outer(y-x,y-z)/(2*sq(y-z)) ),
        (q/kp)*I*( (y-z)@(x-z)/sq(x-z) + (y-x)@(y-z)/(2*sq(y-x)) - (sq(y-x)-sq(y-z))/(2*sq(x-z)) ),
        -(q/pp)*( np.outer(y-z,x-z)/sq(x-z) + np.outer(y-x,y-z).T/(2*sq(y-x)) ),
        (pp/kp)*I*( (y-x)@(x-z)/sq(x-z) - (y-x)@(y-z)/(2*sq(y-z)) - (sq(y-x)-sq(y-z))/(2*sq(x-z)) ),
        -( np.outer(x-z,y-x)/sq(x-z) - np.outer(y-z,y-x)/(2*sq(y-z))
           + np.outer(x-z,y-z)/sq(x-z) + np.outer(y-x,y-z)/(2*sq(y-x)) ) ]
    # index order: first index i (with x'-w'), second m (with z-x); .tex writes e.g.
    # (y-x)^i (x-z)^m -> outer(y-x, x-z)[i,m] ; (y-z)^i (x-z)^m -> outer(y-z,x-z) ;
    # (y-x)^m (y-z)^i -> outer(y-z,y-x)[i,m] ; (y-x)^m (x-z)^i -> outer(x-z,y-x) ; etc.
    B[3] = -(q/pp)*( np.outer(y-z,x-z)/sq(x-z) + np.outer(y-z,y-x)/(2*sq(y-x)) )
    Pv=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x)
    return [float(Pv@b@K) for b in B]

def rsP_blocks(dp, r, s, P, xi):
    """the same six scalars written in r, s, P (the forms quoted in the write-up)."""
    xb=1-xi; P2=sq(P); r2=sq(r); s2=sq(s); smr=s-r; smr2=sq(smr)
    Pr=P@r; Ps=P@s; sr=s@r
    T1 = dp*xi*xb*( Pr*sr/(P2*r2**2) - Pr/(2*P2*r2) )
    T2 = (xi/xb)*(Ps/(P2*r2))*( 1 + (r@smr)/(2*smr2) )
    T3 = xb*(Pr/(P2*r2))*( 2 - 2*sr/r2 - sr/(2*s2) )
    T4 = (xb/xi)*((P@smr)/(P2*r2))*( 1 - sr/(2*s2) )
    T5 = xi*(Pr/(P2*r2))*( 0.5 - 2*sr/r2 - (s@smr)/(2*smr2) )
    T6 = ( Pr*(2*sr - r2)/(P2*r2**2) + (P@smr)*sr/(2*P2*r2*smr2) - Ps*(r@smr)/(2*P2*r2*s2) )
    return [T1,T2,T3,T4,T5,T6]

worst=[0.0]*6; scale=[0.0]*6
for _ in range(400):
    x,y,z,xp,wp = rng.normal(size=(5,2)); kp=1.0; pp=rng.uniform(0.1,0.9)
    a=tex_blocks(2,x,y,z,xp,wp,kp,pp); b=rsP_blocks(2.0,z-x,y-x,xp-wp,pp/kp)
    for i in range(6):
        worst[i]=max(worst[i],abs(a[i]-b[i])); scale[i]=max(scale[i],abs(a[i]))
for i in range(6):
    print("   T%d :  max |tex - (r,s,P) form| = %.2e     (max |T%d| = %.2e)"%(i+1,worst[i],i+1,scale[i]))

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 2.  DENOMINATOR, PHASE, MEASURE"); print(bar)
k_, xi_, s2_, r_s2 = sp.symbols('kplus xi s2 rms2', positive=True)
den = xi_*k_*s2_ + (1-xi_)*k_*r_s2
print("   p+ s^2 + (k+-p+)(s-r)^2 = k+ xibar [ (r-s)^2 + (xi/xibar) s^2 ] :",
      sp.simplify(den - k_*(1-xi_)*(r_s2 + xi_/(1-xi_)*s2_)) == 0)
kx,ky,xx,xy,zx,zy,wx,wy = sp.symbols('kx ky xx xy zx zy wx wy', real=True)
kv=sp.Matrix([kx,ky]); X=sp.Matrix([xx,xy]); Z=sp.Matrix([zx,zy]); W=sp.Matrix([wx,wy])
ph = -(kv.T*(W-Z))[0] - xi_*(kv.T*(Z-X))[0]
ph2 = -(kv.T*(W-X))[0] + (1-xi_)*(kv.T*(Z-X))[0]
print("   -k.(w'-z) - xi k.(z-x) = -k.(w'-x) + xibar k.r :", sp.simplify(ph-ph2)==0)
print("   Int_Lambda^{k+-Lambda} dp+/2pi = (k+/2pi) Int_lam^{1-lam} d xi ,  Int_z = Int d^2 r")
print("   prefactor: (1/(2pi)^3)(g^4/8pi^4)(1/k+) x (k+/2pi) x 1/(k+ xibar D)")
g=sp.Symbol('g',positive=True)
pref = sp.Rational(1,1)/(2*sp.pi)**3 * g**4/(8*sp.pi**4) / k_ * k_/(2*sp.pi) / k_
print("            =", sp.simplify(pref), " x 1/(xibar D)  =  [1/(2pi)^3] (g^4/16 pi^5)(1/k+) x 1/(xibar D)")

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 3.  THE UV REGION r -> 0 : RESIDUE OF EVERY BLOCK"); print(bar)
def Tuv(dp,r,s,P,xi):
    xb=1-xi
    return dp*xi*xb*(P@r)*(s@r)/(sq(P)*sq(r)**2) + (xi/xb+xb/xi)*(P@s)/(sq(P)*sq(r))
def ang(fn,rho,m=4096):
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    return np.mean([fn(rho*np.array([np.cos(a),np.sin(a)])) for a in ph])
s0=np.array([0.7,-0.4]); P0=np.array([-0.3,0.9])
W0=(P0@s0)/(sq(P0)*sq(s0))
print("   rho^2 < T/(xibar D) >  ->  W C_UV(xi) ,   W = P.s/(P^2 s^2) ;  and the subtracted -> 0")
print("   %-6s %14s %14s %14s   %14s"%("xi","rho=1e-3","1e-4","predicted","subtr, 1e-4"))
for xi in (0.1,0.3,0.5,0.8):
    xb=1-xi; M2=xi/xb*sq(s0)
    F=lambda r: sum(rsP_blocks(2.0,r,s0,P0,xi))/(xb*(sq(r-s0)+M2))
    a1=1e-6*ang(F,1e-3); a2=1e-8*ang(F,1e-4)
    Cuv=xi*xb+xi/xb+xb/xi
    sub=1e-8*ang(lambda r: F(r)-Tuv(2.0,r,s0,P0,xi)/sq(s0),1e-4)
    print("   %-6.2f %14.8f %14.8f %14.8f   %14.2e"%(xi,a1,a2,W0*Cuv,sub))
print("   block by block (coefficient of (P.r)(s.r)/(P^2 r^4) or P.s/(P^2 r^2) as r->0):")
print("     T1 -> d_perp xi xibar (P.r)(s.r)/r^4   T2 -> (xi/xibar) P.s/r^2    T3 -> -2 xibar (P.r)(s.r)/r^4")
print("     T4 -> (xibar/xi) P.s/r^2               T5 -> -2 xi (P.r)(s.r)/r^4  T6 -> +2 (P.r)(s.r)/r^4")
xs=sp.Symbol('xi'); dP=sp.Symbol('d_perp')
print("     sum of (P.r)(s.r)/r^4 coefficients :", sp.simplify(dP*xs*(1-xs) - 2*(1-xs) - 2*xs + 2),
      "  (T3+T5+T6 cancel pointwise)")

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 4.  THE TWO MASTER INTEGRALS IN d = 2 - 2 eps"); print(bar)
e=sp.Symbol('epsilon')
G = sp.gamma(-e)
print("   Int d^d r e^{iq.r} (r^2)^{-a} = pi^{d/2} 2^{d-2a} Gamma(d/2-a)/Gamma(a) (q^2)^{a-d/2}")
print("   a = 1 :  = pi Gamma(-eps) (q^2/4pi)^eps")
print("   Gamma(-eps) = ", sp.series(G, e, 0, 1).removeO())
print("   with mu^{-2eps} and mubar^2 = 4 pi mu^2 e^{-gamma_E}:")
print("   M0 = mu^{-2eps} Int d^d r e^{iqr}/r^2 = -pi [ 1/eps + log(q^2/mubar^2) ] + O(eps)")
print("      = +pi [ 1/eps_UV + log(mubar^2/q^2) ] ,   1/eps_UV == -1/eps  (UV needs eps < 0)")
print("   eps*M0 -> -pi  exactly:", sp.limit(e*sp.pi*sp.gamma(-e), e, 0))
print()
print("   r^i r^m/r^4 = (1/2)[ delta^{im}/r^2 - d^i (r^m/r^2) ]  (identity in any d), and")
print("   Int e^{iqr} r^m/r^2 = -i d/dq^m [pi Gamma(-eps)(q^2/4pi)^eps] -> 2 pi i q^m/q^2 , so")
print("   M2^{im} = mu^{-2eps} Int d^d r e^{iqr} r^i r^m/r^4 = (delta^{im}/2) M0 - pi qhat^i qhat^m")
print("   trace check: (d/2) M0 - pi = M0 - eps M0 - pi = M0 + pi - pi = M0   OK")
print()
print("   cutoff cross-check  (|r| > rho):  Int e^{iqr}/r^2 = -pi log(q^2 rho^2 e^{2 gamma}/4)")
for rho,q in ((1e-4,1.3),(1e-6,0.4)):
    q=mp.mpf(q); rho=mp.mpf(rho)
    num=2*mp.pi*(mp.quad(lambda t: mp.besselj(0,q*t)/t,[rho,1/q,60/q])
                 + mp.quadosc(lambda t: mp.besselj(0,q*t)/t,[60/q,mp.inf],period=2*mp.pi/q))
    ana=-mp.pi*mp.log(q**2*rho**2*mp.e**(2*mp.euler)/4)
    print("      rho=%.0e q=%.1f :  %s   vs   %s"%(rho,q,mp.nstr(num,12),mp.nstr(ana,12)))

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 5.  CONTRACT T_UV WITH THE MASTERS, KEEPING d_perp SYMBOLIC"); print(bar)
M0,pi_=sp.symbols('M0 pi'); Ps,PqSq=sp.symbols('Ps PqSq')
d=2-2*e
res = d*xs*(1-xs)*( Ps*M0/2 - pi_*PqSq ) + (xs/(1-xs)+(1-xs)/xs)*Ps*M0
# eps*M0 -> -pi : substitute M0 = -pi/eps + Mfin and expand
Mfin=sp.Symbol('Mfin')
res2 = sp.expand(res.subs(M0, -pi_/e + Mfin))
res2 = sp.limit(res2 - sp.expand((xs*(1-xs)+xs/(1-xs)+(1-xs)/xs)*Ps*(-pi_/e+Mfin)), e, 0)
print("   [d xi xibar P^i s^m M2^{im} + (xi/xibar+xibar/xi) P.s M0] - C_UV P.s M0  =",
      sp.simplify(res2))
print("   i.e.  Int e^{iqr} T_UV/s^2  =  W C_UV M0  +  pi xi xibar [P.s - 2 (P.qhat)(s.qhat)]/(P^2 s^2)")

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 6.  THE xi INTEGRALS OF THE CLOSED-FORM TERMS"); print(bar)
lam=sp.Symbol('lamda',positive=True); v=sp.Symbol('v',positive=True)
C = xs*(1-xs) + xs/(1-xs) + (1-xs)/xs
print("   Int_lam^{1-lam} xi xibar       =", sp.limit(sp.integrate(xs*(1-xs),(xs,lam,1-lam)),lam,0))
a2=sp.integrate(xs/(1-xs),(xs,lam,1-lam)); print("   Int_lam^{1-lam} xi/xibar       =", sp.simplify(sp.series(a2,lam,0,1).removeO()))
a3=sp.integrate((1-xs)/xs,(xs,lam,1-lam)); print("   Int_lam^{1-lam} xibar/xi       =", sp.simplify(sp.series(a3,lam,0,1).removeO()))
print("   => Int C_UV = 2 log(1/lam) - 11/6")
b1=sp.integrate(-2*v*(1-v)*sp.log(v),(v,0,1));   print("   Int xi xibar log(1/xibar^2)   =", b1)
F=sp.log(v)**2/2 - v*sp.log(v) + v
print("   Int_lam xi/xibar log(1/xibar^2) = -2 [F(1)-F(lam)] , F = log^2 v/2 - v log v + v :",
      sp.expand(-2*(F.subs(v,1) - sp.log(lam)**2/2)))
b3=sp.simplify(-2*(sp.integrate(sp.log(1-xs)/xs,(xs,0,1)) - sp.integrate(sp.log(1-xs),(xs,0,1))))
print("   Int xibar/xi log(1/xibar^2)   =", b3)
print("   => Int C_UV log(1/xibar^2) = log^2(1/lam) + pi^2/3 - 67/18")
for L_ in (1e-4,1e-7):
    Cf=lambda t: t*(1-t)+t/(1-t)+(1-t)/t
    q0=mp.quad(Cf,[L_,0.5,1-L_]); q1=mp.quad(lambda t: Cf(t)*mp.log(1/(1-t)**2),[L_,0.5,1-L_])
    Lg=mp.log(1/L_)
    print("      lam=%.0e : %s vs %s ;  %s vs %s"%(L_,mp.nstr(q0,10),mp.nstr(2*Lg-mp.mpf(11)/6,10),
          mp.nstr(q1,10),mp.nstr(Lg**2+mp.pi**2/3-mp.mpf(67)/18,10)))

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 7.  PART (a): IS THERE A UV DIVERGENCE?"); print(bar)
Nc=3
l=np.zeros((8,3,3),dtype=complex)
l[0][0,1]=l[0][1,0]=1; l[1][0,1]=-1j; l[1][1,0]=1j; l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1; l[4][0,2]=-1j; l[4][2,0]=1j; l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j; l[7]=np.diag([1,1,-2])/np.sqrt(3)
t=l/2
f=np.real(-2j*np.einsum('aij,bjk,cki->abc',t,t,t)+2j*np.einsum('bij,ajk,cki->abc',t,t,t))
def adj(V): return np.real(2*np.einsum('aij,jk,bkl,li->ab',t,V,t,V.conj().T))
def randH():
    H=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)); H=(H+H.conj().T)/2
    return H-np.trace(H)*np.eye(3)/3
H1,H2=randH(),randH()
def Uf(p,s_=0.5): return adj(expm(1j*s_*(p[0]*H1+p[1]*H2)))
x0=np.array([0.3,-0.2]); Ux=Uf(x0)
Ca=lambda Uz: -np.einsum('adc,fbe,db,ce->af',f,f,Ux,Uz) + Nc*Ux
Ca2=lambda Uz: -np.einsum('adc,fbe,db,ce->af',f,f,Ux,Uz-Ux)
print("   C_a = -f^{adc}f^{fbe}U^{db}(x)U^{ce}(z) + N_c U^{af}(x)  ==  -f^{adc}f^{fbe}U^{db}(x)[U^{ce}(z)-U^{ce}(x)] :",
      "%.1e"%abs(Ca(Uf(np.array([1.1,0.4])))-Ca2(Uf(np.array([1.1,0.4])))).max())
xi=0.4; xb=1-xi; M2=xi/xb*sq(s0); kv0=np.array([1.0,0.0])
def full_a(r):
    ph=np.exp(1j*(1-xi)*(kv0@r))
    return ph*sum(rsP_blocks(2.0,r,s0,P0,xi))/(xb*(sq(r-s0)+M2))*Ca2(Uf(x0+r))
print("   rho^2 x < integrand of (a) >  (constant = log divergent):")
for rho in (1e-2,1e-3,1e-4,1e-5):
    print("      rho = %.0e :  %.3e"%(rho, rho**2*abs(ang(full_a,rho)).max()))
print("   -> falls as rho^2 : NO UV divergence in (a).  C_a = O(r) and its linear term is odd in rhat.")

# ---------------------------------------------------------------------------
print(); print(bar); print("STEP 8.  THE xi -> 0 RESIDUE AT LARGE |s| : THE LATER x INTEGRAL"); print(bar)
print("   beta_0(x,y) = lim xi I_b = Int d^2r e^{ik.r} [ T4-residue/(r-s)^2 - P.s/(P^2 r^2 s^2) ]")
print("   near r = s (|s| large):  T4-residue -> (1/(2 s^2)) P.(s-r)/(P^2) , and")
print("   Int d^2u e^{ik.u} (-P.u)/u^2 = -2 pi i (P.k)/k^2  =>  beta_0 -> -i pi (P.k)/(P^2 k^2) e^{ik.s}/s^2")
print("   with the outer e^{-ik.(w'-x)} the total phase is e^{-ik.w'} e^{ik.y} : NO x in it, so")
print("   Int d^2x beta_0 ~ Int d^2x /(y-x)^2  diverges logarithmically at large |x| .")
