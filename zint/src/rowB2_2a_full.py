"""Part 2a: is there a UV term, and how does the p+ integration really go?

Variables:  eta = p+/k+ ,  etabar = 1-eta .   In B_2(k+-p+,-z ; p+,-x) the line at z
carries etabar k+ and the line at x carries eta k+ .  Total phase
    exp[-i k.(w' - etabar z - eta x)] ,
so the z-oscillation frequency is etabar k and the x-oscillation frequency is eta k .
"""
import numpy as np
from scipy.linalg import expm
import mpmath as mp
mp.mp.dps = 20
rng = np.random.default_rng(20260927)
Nc = 3
l = np.zeros((8,3,3), dtype=complex)
l[0][0,1]=l[0][1,0]=1; l[1][0,1]=-1j; l[1][1,0]=1j; l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1; l[4][0,2]=-1j; l[4][2,0]=1j; l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j; l[7]=np.diag([1,1,-2])/np.sqrt(3)
t = l/2
f = np.real(-2j*np.einsum('aij,bjk,cki->abc',t,t,t)+2j*np.einsum('bij,ajk,cki->abc',t,t,t))
def adj(V): return np.real(2*np.einsum('aij,jk,bkl,li->ab',t,V,t,V.conj().T))
def randH():
    H=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)); H=(H+H.conj().T)/2
    return H-np.trace(H)*np.eye(3)/3
H1,H2 = randH(), randH()
def Ufield(v, s=0.5, R=1.6):
    """a smooth adjoint field that -> 1 at large |v| (the target has finite support)."""
    g = np.exp(-(v@v)/(2*R*R))
    return adj(expm(1j*s*g*(v[0]*H1 + v[1]*H2)))
def C2a(Ux, Uz): return -np.einsum('adc,fbe,db,ce->af', f, f, Ux, Uz-Ux)

def sq(v): return float(v@v)
def Bracket(D,x,y,z,kp,pp):
    I=np.eye(D); q=kp-pp
    T1=I*D/(2*sq(x-z))*(sq(y-x)-sq(y-z))*(pp*q/kp**2)
    T2=-(pp/q)*( np.outer(y-x,x-z)/sq(x-z) - np.outer(y-x,y-z)/(2*sq(y-z)) )
    T3=(q/kp)*I*( (y-z)@(x-z)/sq(x-z) + (y-x)@(y-z)/(2*sq(y-x)) - (sq(y-x)-sq(y-z))/(2*sq(x-z)) )
    T4=-(q/pp)*( np.outer(y-z,x-z)/sq(x-z) + np.outer(y-z,y-x)/(2*sq(y-x)) )
    T5=(pp/kp)*I*( (y-x)@(x-z)/sq(x-z) - (y-x)@(y-z)/(2*sq(y-z)) - (sq(y-x)-sq(y-z))/(2*sq(x-z)) )
    T6=-( np.outer(x-z,y-x)/sq(x-z) - np.outer(y-z,y-x)/(2*sq(y-z))
          + np.outer(x-z,y-z)/sq(x-z) + np.outer(y-x,y-z)/(2*sq(y-x)) )
    return T1+T2+T3+T4+T5+T6
def scal(x,y,z,xp,wp,kvec,kp,pp):
    """transverse x phase, colour stripped (the braces of (1.15) incl. 1/k+)."""
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x); e=pp/kp
    ph=np.exp(-1j*(kvec@(wp-(1-e)*z-e*x)))
    return ph*float(V@Bracket(2,x,y,z,kp,pp)@K)/(pp*sq(y-x)+(kp-pp)*sq(y-z))/kp
def ang(fn,r,m=8192):
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    return np.mean([fn(r*np.array([np.cos(a),np.sin(a)])) for a in ph], axis=0)

x0=np.array([0.30,-0.20]); y0=np.array([-0.7,1.1]); xp0=np.array([1.6,0.4]); wp0=np.array([-0.5,-1.3])
kvec=np.array([0.5,-0.3]); kp=1.0
V0=(xp0-wp0)/sq(xp0-wp0); R0=y0-x0; W0=(V0@R0)/sq(R0)
Ux0=Ufield(x0)

print("="*78)
print("1.  DOES 2a HAVE A UV DIVERGENCE AT z -> x ?   NO.")
print("="*78)
print("   Reason, before any numerics.  The exact colour factor is")
print("       C_2a = -f^{adc} f^{fbe} U^{db}(x) [ U^{ce}(z) - U^{ce}(x) ]  =  O(|z-x|) ,")
print("   and to leading order it is (N_c/2)[U(x) - U(z)] .  Two powers are gained:")
print("     (i)  one from that linear vanishing, against the kernel's 1/(z-x)^2 and d^2z ~ |u| d|u| ;")
print("     (ii) one more because the linear piece is ODD in (z-x)hat while the UV residue")
print("          C_UV(eta) W /(z-x)^2 is EVEN, so the angular average kills it.")
print()
print("   u^2 x < integrand >  as u = z-x -> 0   (constant = log divergent, -> 0 = finite):")
print("   %-10s %18s %18s"%("|z-x|","2a (with colour)","colour stripped"))
for r in (1e-2,1e-3,1e-4,1e-5):
    a=r**2*abs(ang(lambda d: scal(x0,y0,x0+d,xp0,wp0,kvec,kp,0.4)*C2a(Ux0,Ufield(x0+d)),r)).max()
    b=r**2*abs(ang(lambda d: scal(x0,y0,x0+d,xp0,wp0,kvec,kp,0.4),r))
    print("   %-10.0e %18.3e %18.8f"%(r,a,b))
print("   -> falls as |z-x|^2 .  NOTHING TO REGULARISE: no dimensional regularisation, no")
print("      MSbar, no counterterm.  Keep d_perp = 2 throughout 2a.")

print()
print("="*78)
print("2.  WHAT 2a DOES HAVE: A 1/z^2 TAIL, CUT ONLY BY THE PHASE")
print("="*78)
print("   At large |z| , U(z) -> 1 (not U(x)), so C_2a -> -f^{adc}f^{fbc}U^{db}(x) + N_c U^{af}(x) != 0")
print("   and the colour-stripped integrand has the tail")
print("       V(eta)/z^2 ,   V(eta) = (W/2k+) [ etabar/(d_perp eta) + k+/(k+-p+) ] ,")
print("   verified earlier.  Both ENDPOINTS make V(eta) blow up:")
import sympy as sp
e=sp.Symbol('eta',positive=True); dP=sp.Symbol('d_perp',positive=True)
Vf=( (1-e)/(dP*e) + 1/(1-e) )/2
print("       V(eta) x 2k+/W =", sp.simplify(Vf*2))
print("       eta -> 0   :", sp.limit(sp.simplify(Vf*2)*e, e, 0), "/ eta     -> simple 1/eta")
print("       eta -> 1   :", sp.simplify(sp.limit(sp.simplify(Vf*2)*(1-e), e, 1)), "/ etabar -> simple 1/etabar")
print()
print("   The z-oscillation frequency is  etabar k , so the tail integral is")
print("       Int_{|z|>R0} d^2z e^{i etabar k.z}/z^2 = -2 pi log[ etabar k_perp R0 e^{gamma}/2 ]")
print("                                             = 2 pi log(1/etabar) + const .")
print("   It therefore produces an EXTRA log(1/etabar) at the eta -> 1 end, and NOTHING at the")
print("   eta -> 0 end (where etabar -> 1 and the oscillation is fully on).")
for q in (1e-1,1e-2,1e-3):
    num=2*mp.pi*(mp.quad(lambda r: mp.besselj(0,q*r)/r,[1.0,1/q,60/q])
                 + mp.quadosc(lambda r: mp.besselj(0,q*r)/r,[60/q,mp.inf],period=2*mp.pi/q))
    ana=-2*mp.pi*mp.log(q*mp.e**mp.euler/2)
    print("      etabar k_perp = %-8.0e  tail = %14s   2pi log(1/.) form = %14s"
          %(q, mp.nstr(num,9), mp.nstr(ana,9)))

print()
print("="*78)
print("3.  AND A MIRROR 1/x^2 TAIL, CUT BY THE OTHER OSCILLATION  exp[i eta k.x]")
print("="*78)
print("   NOTE on the diagnostic: with the phase IN, the angular average of a 1/v^2 tail is")
print("   ~ J_0(q|v|)/v^2 , which oscillates and decays, so '|v|^2 x <f> -> constant' is the")
print("   wrong test.  Measure the tail COEFFICIENT with the oscillation stripped, then put the")
print("   oscillation back analytically as Int_{R0} d|v|/|v| J_0(q|v|) = log(1/q) + const .")
def scal0(x,y,z,xp,wp,kp,pp):
    """same integrand with the phase set to 1: exposes the tail coefficient."""
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x)
    return float(V@Bracket(2,x,y,z,kp,pp)@K)/(pp*sq(y-x)+(kp-pp)*sq(y-z))/kp
zf=np.array([0.9,0.35])
print()
print("   tail coefficients, phase stripped   (|v|^2 x < f > -> constant)")
print("   %-10s %20s %20s"%("|v|","z -> inf, eta = 0.4","x -> inf, eta = 0.6"))
for X in (1e2,1e3,1e4,1e5):
    a_=X**2*ang(lambda d: scal0(x0,y0,x0+d,xp0,wp0,kp,0.4),X)
    b_=X**2*ang(lambda d: scal0(zf+d,y0,zf,xp0,wp0,kp,0.6),X)
    print("   %-10.0e %20.8f %20.8f"%(X,a_,b_))
print("   -> both are genuine 1/v^2 tails.  (The mirror symmetry (x <-> z, eta <-> etabar) of")
print("      the colour-stripped integrand, verified earlier to 1.6e-13, guarantees this.)")
print()
print("   the z-tail coefficient, checked against the closed form")
print("      V_z(eta) = (W / 2 k+) [ etabar/(d_perp eta) + k+/(k+-p+) ] :")
print("   %-8s %20s %20s"%("eta","measured","closed form"))
for et in (0.2,0.4,0.6,0.8,0.95):
    m=1e4**2*ang(lambda d: scal0(x0,y0,x0+d,xp0,wp0,kp,et),1e4)
    cf=(W0/(2*kp))*((kp/(et*kp)-1)/2 + 1/(1-et))
    print("   %-8.2f %20.8f %20.8f"%(et,m,cf))
print("   -> V_z ~ (W/2k+)/etabar as eta -> 1 ;  by mirror symmetry V_x ~ (W'/2k+)/eta as eta -> 0 .")
print()
print("   and 2a's colour survives in BOTH tails:")
Uzr=adj(expm(1j*0.5*randH())); Uxr=adj(expm(1j*0.5*randH()))
print("      || C_2a(U(z) -> 1) || = %.4f     || C_2a(U(x) -> 1) || = %.4f"
      %(abs(C2a(Uxr,np.eye(8))).max(), abs(C2a(np.eye(8),Uzr)).max()))
print("   so neither tail is switched off by the colour structure.")

print("="*78)
print("4.  SO THE p+ INTEGRAND HAS log/POLE, NOT JUST 1/POLE")
print("="*78)
print("   Writing G(eta) = Int_{x,y,z,x'} (2a integrand) , the two endpoint behaviours are")
print()
print("       G(eta) ->  [ A_0 log(1/eta)    + B_0 ] / eta        as eta -> 0 ,")
print("       G(eta) ->  [ A_1 log(1/etabar) + B_1 ] / etabar     as eta -> 1 .")
print()
print("   A SINGLE plus prescription is NOT enough.  g(endpoint) for the naive [1/eta]_+ is")
print("   log-divergent (that is exactly the statement A_0 != 0), and the plus prescription")
print("   requires g(endpoint) FINITE.  You need the two standard distributions at each end:")
print()
print("     1/eta            = [1/eta]_+            + delta(eta) log(k+/Lambda)")
print("     log(1/eta)/eta   = [log(1/eta)/eta]_+   + (1/2) delta(eta) log^2(k+/Lambda)")
print()
print("   Check of both, against the direct cut-off integral Int_lam^1 :")
def d1(g,lam): return mp.quad(lambda v: g(v)/v, [lam, 1])
def d1p(g,lam):
    return mp.quad(lambda v: (g(v)-g(0))/v, [0,1]) + g(0)*mp.log(1/lam)
def d2(g,lam): return mp.quad(lambda v: g(v)*mp.log(1/v)/v, [lam, 1])
def d2p(g,lam):
    return mp.quad(lambda v: (g(v)-g(0))*mp.log(1/v)/v, [0,1]) + g(0)*mp.log(1/lam)**2/2
g=lambda v: mp.e**(-1.3*v)*mp.cos(0.7*v)
print("   %-10s %18s %18s %18s %18s"%("lam","1/eta direct","1/eta plus","log/eta direct","log/eta plus"))
for lam in (1e-3,1e-6,1e-9):
    print("   %-10.0e %18.10f %18.10f %18.10f %18.10f"
          %(lam, d1(g,lam), d1p(g,lam), d2(g,lam), d2p(g,lam)))
print("   -> both identities hold, with O(lam) error.  The second one is where log^2(k+/Lambda)")
print("      comes from, and it is the reason this row has a DOUBLE rapidity logarithm.")

print()
print("="*78)
print("5.  THE ASSEMBLED RESULT FOR 2a")
print("="*78)
print("""   Write the z- and x-integrals first, at d_perp = 2 (2a needs no regulator), and call

       G(eta) = Int_{x,y,z,x'} e^{-i k.(w' - etabar z - eta x)} x (braces of (1.15)) x (colour) .

   Its endpoint behaviour, established above, is

       G -> [ A_0 log(1/eta)    + B_0 ] / eta      (eta -> 0 ; from the 1/x^2 tail)
       G -> [ A_1 log(1/etabar) + B_1 ] / etabar   (eta -> 1 ; from the 1/z^2 tail)

       A_1 = (pi/k+) x (x'-w').(y-x)/[(x'-w')^2 (y-x)^2] x [ -f^{adc}f^{fbc}U^{db}(x)
                                                             + N_c U^{af}(x) ]        (z -> oo)
       A_0 = the same with x <-> z and U(x) -> 1 in the colour factor.

   Subtract them and integrate:

     Int_Lambda^{k+-Lambda} dp+/2pi  G  =  (k+/2pi) {
          (A_0 + A_1)/2 x log^2(k+/Lambda)
        + (B_0 + B_1)   x log(k+/Lambda)
        + Int_0^1 d eta [ G(eta) - (A_0 log(1/eta)+B_0)/eta
                                 - (A_1 log(1/etabar)+B_1)/etabar ]  }  + O(Lambda)

   In plus-prescription language the same statement is

       1/eta            = [1/eta]_+          + delta(eta) log(k+/Lambda)
       log(1/eta)/eta   = [log(1/eta)/eta]_+ + (1/2) delta(eta) log^2(k+/Lambda)

   and likewise in etabar.  FOUR distributions in total, two per endpoint.
""")
print("   Consistency with 2b: 2b's counterterm moment I1 = -L^2 + 67/18 - pi^2/3 also carried a")
print("   log^2 .  Both halves of the row have one, and their sum is what has to match the")
print("   double logarithm of the JIMWLK/Sudakov structure.")
