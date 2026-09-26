"""B_2 x C_1 row, part 2a: can the plus prescription do the p+ integration?

eta = p+/k+ , etabar = 1-eta , A = (y-x)^2 , B = (y-z)^2 , kappa = k.(z-x) , lam = Lambda/k+ .
The six blocks of eq (3) carry the p+ weights
   T1: eta etabar      T2: -eta/etabar     T3: etabar
   T4: -etabar/eta     T5: eta             T6: -1
and the ONLY other p+ dependence is the denominator  k+[eta A + etabar B]  and the phase
e^{-i eta kappa} .  The six tensors are p+-independent.
"""
import numpy as np, mpmath as mp, sympy as sp
mp.mp.dps = 25

print("="*78)
print("1.  WHERE THE p+ INTEGRAND IS SINGULAR")
print("="*78)
eta = sp.Symbol('eta', positive=True); eb = 1 - eta
cs = {'T1': eta*eb, 'T2': -eta/eb, 'T3': eb, 'T4': -eb/eta, 'T5': eta, 'T6': sp.Integer(-1)}
print("   %-5s %-14s %-16s %-16s"%("block","weight","eta -> 0","eta -> 1"))
for n,c in cs.items():
    l0 = sp.limit(c, eta, 0); l1 = sp.limit(c, eta, 1)
    print("   %-5s %-14s %-16s %-16s"%(n, c, l0, l1))
print()
print("   -> TWO endpoint poles, not one:  T4 ~ -1/eta at p+ -> 0 ,  T2 ~ -1/etabar at p+ -> k+ .")
print("   This is the signature of a LOOP (self-energy) topology: C_1 splits k+ into p+ and")
print("   k+-p+ and B_2 re-merges them, so both soft ends are singular.  Row IV had only one,")
print("   because there the measured gluon kept k+ and p+ ran up to V-k+ .")

print()
print("="*78)
print("2.  IS THE PLUS PRESCRIPTION LEGITIMATE?  (need g(endpoint) FINITE)")
print("="*78)
print("   The residue of each pole is (tensor) x h(eta) , with")
print("       h(eta) = e^{-i eta kappa} / [ eta A + etabar B ] ,")
print("   and h is smooth on [0,1] because A = (y-x)^2 > 0 and B = (y-z)^2 > 0 , so")
print("       eta A + etabar B  >=  min(A,B)  >  0   : NO PINCH.")
print("       h(0) = 1/B = 1/(y-z)^2          finite")
print("       h(1) = e^{-i kappa}/A = e^{-i k.(z-x)}/(y-x)^2   finite")
print("   => BOTH plus prescriptions are valid.  Answer: YES, and you need it at both ends.")

print()
print("="*78)
print("3.  THE TWO-SIDED PLUS PRESCRIPTION, CHECKED NUMERICALLY")
print("="*78)
def build(A, B, kappa, T):
    """T = dict of the six contracted tensor values (p+-independent numbers)."""
    def h(e): return mp.e**(-1j*e*kappa)/(e*A + (1-e)*B)
    def integ(e):
        w = {'T1': e*(1-e), 'T2': -e/(1-e), 'T3': 1-e, 'T4': -(1-e)/e, 'T5': e, 'T6': -1.0}
        return h(e)*sum(w[n]*T[n] for n in T)
    return h, integ

def direct(A,B,kappa,T,lam):
    h, integ = build(A,B,kappa,T)
    return mp.quad(integ, [lam, 0.5, 1-lam])

def plusform(A,B,kappa,T,lam):
    h, integ = build(A,B,kappa,T)
    L = mp.log((1-lam)/lam)
    # delta pieces: -1/eta x (etabar h T4) at eta=0 ; -1/etabar x (eta h T2) at eta=1
    delta = -L*h(0)*T['T4'] - L*h(1)*T['T2']
    def reg(e):
        g0 = (1-e)*h(e)*T['T4']; g0_0 = h(0)*T['T4']          # companion of 1/eta
        g1 = e*h(e)*T['T2'];     g1_1 = h(1)*T['T2']          # companion of 1/etabar
        rest = h(e)*( e*(1-e)*T['T1'] + (1-e)*T['T3'] + e*T['T5'] - T['T6'] )
        # note: T6 weight is -1, so -T6 above reproduces it with the sign below
        rest = h(e)*( e*(1-e)*T['T1'] + (1-e)*T['T3'] + e*T['T5'] - T['T6'] )
        return -(g0-g0_0)/e - (g1-g1_1)/(1-e) + rest
    return delta + mp.quad(reg, [0, 0.5, 1])

rng = np.random.default_rng(3)
print("   %-12s %26s %26s %12s"%("Lambda/k+","direct integral","two-sided plus form","rel. diff"))
A, B = 1.7, 0.9; kappa = 0.8
T = {n: complex(rng.normal(), rng.normal()) for n in ('T1','T2','T3','T4','T5','T6')}
for lam in (1e-3,1e-5,1e-7,1e-9):
    d = direct(A,B,kappa,T,lam); p = plusform(A,B,kappa,T,lam)
    print("   %-12.0e %26s %26s %12.2e"
          %(lam, mp.nstr(d,10), mp.nstr(p,10), abs(d-p)/abs(d)))
print("   -> agreement improves as Lambda -> 0 ; the residual is the O(Lambda) tail from")
print("      extending the regular integrals down to 0 and up to 1.")

print()
print("="*78)
print("4.  ARE THE TWO POLES MIRROR IMAGES?  (x <-> z together with p+ <-> k+-p+)")
print("="*78)
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
def scalar(x,y,z,xp,wp,kvec,kp,pp):
    """transverse x phase, colour stripped.  Total phase exp[-i k.(w' - etabar z - eta x)]."""
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x); e=pp/kp
    ph=np.exp(-1j*(kvec@(wp - (1-e)*z - e*x)))
    return ph*float(V@Bracket(2,x,y,z,kp,pp)@K)/(pp*sq(y-x)+(kp-pp)*sq(y-z))/kp
rng2=np.random.default_rng(99)
worst=0.0; scale=0.0
for _ in range(200):
    x,y,z,xp,wp,kvec = rng2.normal(size=(6,2))
    kp=1.0; pp=rng2.uniform(0.15,0.85)
    a=scalar(x,y,z,xp,wp,kvec,kp,pp)
    b=scalar(z,y,x,xp,wp,kvec,kp,kp-pp)          # x <-> z  AND  p+ -> k+-p+
    worst=max(worst,abs(a-b)); scale=max(scale,abs(a))
print("   max | F(x,y,z;p+) - F(z,y,x;k+-p+) | = %.3e      (for scale, max|F| = %.3f)"%(worst,scale))
print("   -> the transverse-times-phase factor IS invariant.  The total phase is")
print("        exp[-i k.(w' - etabar z - eta x)] , which is manifestly symmetric.")
print("   So the two endpoint poles are MIRROR IMAGES, exchanged by x <-> z .")
print()
print("   The COLOUR factors are NOT symmetric under x <-> z , though:")
print("     original  f^{adc} f^{ebf} U^{db}(x) U^{ce}(z) + N_c U^{af}(y)   -> x <-> z swaps the lines")
print("     2a        -f^{adc} f^{fbe} U^{db}(x)[U^{ce}(z) - U^{ce}(x)]")
print("     2b        N_c [ U^{af}(y) - U^{af}(x) ]")
print("   so you may NOT use the symmetry to fold one pole onto the other.  Keep both.")

print()
print("="*78)
print("5.  THE RESULT, AND WHAT NOT TO DO NEXT")
print("="*78)
print("   With  eta = p+/k+ ,  A = (y-x)^2 ,  B = (y-z)^2 ,  kappa = k.(z-x) ,")
print("   h(eta) = e^{-i eta kappa} / [eta A + etabar B] , and Int_Lambda^{k+-Lambda} dp+ = k+ Int d eta :")
print()
print("     k+ Int_0^1 d eta {  - [1/eta]_+  etabar h(eta) T4  - [1/etabar]_+ eta h(eta) T2")
print("                         + h(eta) [ eta etabar T1 + etabar T3 + eta T5 - T6 ]  }")
print("     - k+ log(k+/Lambda) [  T4/(y-z)^2  +  e^{-i k.(z-x)} T2/(y-x)^2  ]")
print()
print("   TWO delta-type terms, so 2 log(k+/Lambda) worth of rapidity logarithm, not one.")
print()
print("   DO NOT change variables to xi = k+/(k+ + p+) here.  In this row p+ <= k+ , so")
print("   xi in [1/2, 1] , and the measured momentum stays k -- the phase is")
print("   exp[-i k.(w' - etabar z - eta x)] , with no k/xi anywhere.  There is no")
print("   fragmentation measure d xi/xi^2 and no hard factor at k/xi , so this row carries")
print("   NO DGLAP convolution: it is a virtual correction, and everything it produces is")
print("   delta(1-xi)-like.  P_gg showed up in it only as the TRANSVERSE UV residue (part 2b).")

print()
print("="*78)
print("6.  DO THE TWO delta COEFFICIENTS STAY UV FINITE IN 2a ?")
print("="*78)
from scipy.linalg import expm
Nc=3
l=np.zeros((8,3,3),dtype=complex)
l[0][0,1]=l[0][1,0]=1; l[1][0,1]=-1j; l[1][1,0]=1j; l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1; l[4][0,2]=-1j; l[4][2,0]=1j; l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j; l[7]=np.diag([1,1,-2])/np.sqrt(3)
t=l/2
f=np.real(-2j*np.einsum('aij,bjk,cki->abc',t,t,t)+2j*np.einsum('bij,ajk,cki->abc',t,t,t))
def adj(V): return np.real(2*np.einsum('aij,jk,bkl,li->ab',t,V,t,V.conj().T))
def randH():
    H=rng2.normal(size=(3,3))+1j*rng2.normal(size=(3,3)); H=(H+H.conj().T)/2
    return H-np.trace(H)*np.eye(3)/3
H1,H2=randH(),randH()
def Uf(v,s=0.5): return adj(expm(1j*s*(v[0]*H1+v[1]*H2)))
def T4res(x,y,z,xp,wp,kp):
    """the eta -> 0 delta coefficient: T4 tensor contracted, divided by (y-z)^2."""
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x)
    M=-( np.outer(y-z,x-z)/sq(x-z) + np.outer(y-z,y-x)/(2*sq(y-x)) )
    return float(V@M@K)/sq(y-z)/kp
def T2res(x,y,z,xp,wp,kvec,kp):
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x)
    M=-( np.outer(y-x,x-z)/sq(x-z) - np.outer(y-x,y-z)/(2*sq(y-z)) )
    return np.exp(-1j*(kvec@(z-x)))*float(V@M@K)/sq(y-x)/kp
x0=np.array([0.3,-0.2]); y0=np.array([-0.7,1.1]); xp0=np.array([1.6,0.4]); wp0=np.array([-0.5,-1.3])
kvec=np.array([0.5,-0.3]); Ux=Uf(x0)
def C2a(Ux,Uz): return -np.einsum('adc,fbe,db,ce->af',f,f,Ux,Uz-Ux)
def C2b(Ux,Uy): return Nc*(Uy-Ux)
Uy=Uf(y0)
def ang(fn,r,m=2048):
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    return np.mean([fn(r*np.array([np.cos(a),np.sin(a)])) for a in ph],axis=0)
print("   u^2 x < delta-coefficient x colour >   (constant = log divergent, -> 0 = finite)")
print("   %-10s %14s %14s %14s %14s"%("|z-x|","2a  eta->0","2a  eta->1","2b  eta->0","2b  eta->1"))
for r in (1e-2,1e-3,1e-4):
    a0=r**2*abs(ang(lambda d: T4res(x0,y0,x0+d,xp0,wp0,1.0)*C2a(Ux,Uf(x0+d)),r)).max()
    a1=r**2*abs(ang(lambda d: T2res(x0,y0,x0+d,xp0,wp0,kvec,1.0)*C2a(Ux,Uf(x0+d)),r)).max()
    b0=r**2*abs(ang(lambda d: T4res(x0,y0,x0+d,xp0,wp0,1.0)*C2b(Ux,Uy),r)).max()
    b1=r**2*abs(ang(lambda d: T2res(x0,y0,x0+d,xp0,wp0,kvec,1.0)*C2b(Ux,Uy),r)).max()
    print("   %-10.0e %14.3e %14.3e %14.8f %14.8f"%(r,a0,a1,b0,b1))
print("   -> in 2a both delta coefficients are UV FINITE: the plus prescription introduces")
print("      no new short-distance problem, because C_2a = O(|z-x|) still does its job.")
print("   -> in 2b they are NOT: there the UV log and the rapidity log MULTIPLY, giving")
print("      log_UV x log(k+/Lambda) -- consistent with C_UV(eta) = eta etabar + eta/etabar")
print("      + etabar/eta having poles at BOTH endpoints.  Handle 2b's UV first, then p+ .")
