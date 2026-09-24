"""New row (B_2 x C_1): check eq(2) and eq(3) against the quoted vertex definitions.

Notation used here (all vectors in D transverse dimensions, d_perp = D):
   B_{jk}   = the bracket of B_2^{be}_{2jk}(k+-p+, -z ; p+, -x) after the substitutions
   E_{im}   = C_1's tensor contracted with B:  delta_im trB - (k+/p+) B_mi - (k+/(k+-p+)) B_im
   Br_{im}  = the bracket the user WROTE in eq (2)
"""
import numpy as np
rng = np.random.default_rng(20260924)

def sq(v): return float(v @ v)

def Bmat(D, x, y, z, kp, pp):
    """B_2 bracket with (x_def,z_def,w_def,k+_def,i_def) -> (y,x,z,k+-p+,k)."""
    I = np.eye(D)
    kq = kp - pp                       # the 'k+' slot of the definition
    tot = pp + kq                      # = k+   (the (p+ + k+) of the definition)
    B = I/(2*sq(x-z))*(sq(y-x) - sq(y-z))
    B += (tot/kq)*( np.outer(y-x, x-z)/sq(x-z) - np.outer(y-x, y-z)/(2*sq(y-z)) )
    B += (tot/pp)*( np.outer(x-z, y-z)/sq(x-z) + np.outer(y-x, y-z)/(2*sq(y-x)) )
    return B                            # B[j,k]

def Econtract(D, x, y, z, kp, pp):
    I = np.eye(D)
    B = Bmat(D, x, y, z, kp, pp)
    return I*np.trace(B) - (kp/pp)*B.T - (kp/(kp-pp))*B     # [i,m]

def Bracket_user(D, x, y, z, kp, pp):
    """exactly the six blocks of eq (2), as typed."""
    I = np.eye(D); q = kp - pp
    T1 = I*D/(2*sq(x-z))*(sq(y-x)-sq(y-z))
    T2 = -(kp**2/q**2)*( np.outer(y-x, x-z)/sq(x-z) - np.outer(y-x, y-z)/(2*sq(y-z)) )
    T3 = (kp/pp)*I*( (y-z)@(x-z)/sq(x-z) + (y-x)@(y-z)/(2*sq(y-x))
                     - (sq(y-x)-sq(y-z))/(2*sq(x-z)) )
    T4 = -(kp**2/pp**2)*( np.outer(y-z, x-z)/sq(x-z) + np.outer(y-z, y-x)/(2*sq(y-x)) )
    T5 = (kp/q)*I*( (y-x)@(x-z)/sq(x-z) - (y-x)@(y-z)/(2*sq(y-z))
                    - (sq(y-x)-sq(y-z))/(2*sq(x-z)) )
    T6 = -(kp**2/(pp*q))*( np.outer(x-z, y-x)/sq(x-z) - np.outer(y-z, y-x)/(2*sq(y-z))
                           + np.outer(x-z, y-z)/sq(x-z) + np.outer(y-x, y-z)/(2*sq(y-x)) )
    return T1+T2+T3+T4+T5+T6            # [i,m]

print("="*78)
print("1.  THE TENSOR BRACKET OF EQ (2), FROM C_1 x B_2")
print("="*78)
print("   %-4s %14s %14s %12s"%("D","max|derived|","max|quoted|","max|diff|"))
for D in (2,3,4,5):
    md=mq=mdiff=0.0
    for _ in range(300):
        x,y,z = rng.normal(size=(3,D))
        kp = 1.0; pp = rng.uniform(0.15,0.85)
        E = Econtract(D,x,y,z,kp,pp); Bu = Bracket_user(D,x,y,z,kp,pp)
        md=max(md,abs(E).max()); mq=max(mq,abs(Bu).max()); mdiff=max(mdiff,abs(E-Bu).max())
    print("   %-4d %14.4e %14.4e %12.3e"%(D,md,mq,mdiff))
print("   -> the six blocks of eq (2) ARE the C_1 x B_2 contraction. Exact.")

import sympy as sp
print()
print("="*78)
print("2.  THE PREFACTOR OF EQ (2)")
print("="*78)
g,pi,kp,pp = sp.symbols('g pi kplus pplus', positive=True)
I = sp.I
A_dag  = -I*g/(sp.sqrt(2)*pi*sp.sqrt(kp))                       # (A^(1))^dagger magnitude
C_dag  = -g/(2*pi*sp.sqrt(2*kp))*sp.sqrt(pp*(kp-pp))/kp         # C_1^dagger  (real)
B_mag  =  I*g**2*sp.sqrt(pp*(kp-pp))/(4*pi**2*kp)               # B_2, with (p+ + k+) -> k+
prod = sp.simplify(A_dag*C_dag*B_mag)
print("   A^dag        =", A_dag)
print("   C_1^dag      =", C_dag)
print("   B_2          =", B_mag, "   / [p+(y-x)^2 + (k+-p+)(y-z)^2]")
print("   product      =", prod)
tot = sp.simplify(2/(2*pi)**3 * prod)            # eq (1) carries an overall 2/(2pi)^3
print("   x 2/(2pi)^3  =", tot)
quoted = sp.simplify(1/(2*pi)**3 * g**4/(8*pi**4) * pp*(kp-pp)/kp**3)
print("   eq (2) says  =", quoted)
print("   ratio quoted/derived =", sp.simplify(quoted/tot))
print("   -> magnitude g^4/(8 pi^4), 1/k+^3 and p+(k+-p+) all EXACT; the sign is -1,")
print("      which is restored by the first bracket:  -Abar^dag + A^dag U  ~  -[U(x')-U(w')].")

print()
print("="*78)
print("3.  THE COLOUR ALGEBRA")
print("="*78)
Nc = 3
# Gell-Mann matrices -> t^a = lambda^a/2 -> f^{abc} = -2i tr([t^a,t^b] t^c)
l = np.zeros((8,3,3), dtype=complex)
l[0][0,1]=l[0][1,0]=1
l[1][0,1]=-1j; l[1][1,0]=1j
l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1
l[4][0,2]=-1j; l[4][2,0]=1j
l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j
l[7]=np.diag([1,1,-2])/np.sqrt(3)
t = l/2
f = np.real(-2j*np.einsum('aij,bjk,cki->abc', t, t, t)
            +2j*np.einsum('bij,ajk,cki->abc', t, t, t))
print("   f built from Gell-Mann; antisymmetry check:", abs(f + f.transpose(1,0,2)).max())
print("   f^{abc} f^{dbc} = N_c delta^{ad} ?  max dev =",
      abs(np.einsum('abc,dbc->ad', f, f) - Nc*np.eye(8)).max())
print("   f^{abc} f^{cbd} = -N_c delta^{ad} ? max dev =",
      abs(np.einsum('abc,cbd->ad', f, f) + Nc*np.eye(8)).max())
print("   so  - f^{abc} f^{cbd} U^{df}(y)  =  + N_c U^{af}(y) .   CORRECT as quoted.")
print()
from scipy.linalg import expm
def adjU(rng, s=0.6):
    """a genuine adjoint SU(3) Wilson line: U^{ab} = 2 tr(t^a V t^b V^dag)."""
    H = rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3))
    H = (H + H.conj().T)/2; H = H - np.trace(H)*np.eye(3)/3
    V = expm(1j*s*H)
    return np.real(2*np.einsum('aij,jk,bkl,li->ab', t, V, t, V.conj().T))
Ux, Uz = adjU(rng), adjU(rng)
print("   adjointness check  f^{def} U^{da} U^{eb} U^{fc} = f^{abc} :  max dev =",
      abs(np.einsum('def,da,eb,fc->abc', f, Ux, Ux, Ux) - f).max())
M_gen  = np.einsum('adc,ebf,db,ce->af', f, f, Ux, Uz)
M_same = np.einsum('adc,ebf,db,ce->af', f, f, Ux, Ux)
print()
print("   f^{adc} f^{ebf} U^{db}(x) U^{ce}(z) , compared with -N_c U(x):")
print("     U(z) = U(x) :  ||M + N_c U^{af}(x)|| = %.3e    ||M + N_c U^{fa}(x)|| = %.3e"
      %(abs(M_same + Nc*Ux).max(), abs(M_same + Nc*Ux.T).max()))
print("     U(z) free   :  ||M + N_c U^{af}(x)|| = %.3e    ||M + N_c U^{fa}(x)|| = %.3e"
      %(abs(M_gen + Nc*Ux).max(), abs(M_gen + Nc*Ux.T).max()))
print("     scale for comparison: ||M|| = %.3f , N_c||U|| = %.3f"%(abs(M_gen).max(), Nc*abs(Ux).max()))

print()
print("="*78)
print("4.  EQ (2) -> EQ (3): the w integration and the p+ redistribution")
print("="*78)
kps, pps = sp.symbols('kp pp', positive=True); q = kps - pps
pre2 = pps*q/kps**3                                  # eq (2) prefactor weight
c2 = [sp.Integer(1), -kps**2/q**2, kps/pps, -kps**2/pps**2, kps/q, -kps**2/(pps*q)]
c3 = [pps*q/kps**2, -pps/q, q/kps, -q/pps, pps/kps, sp.Integer(-1)]   # eq (3), with 1/k+ out front
print("   %-6s %-22s %-22s %s"%("block","eq(2) coeff x p+(k+-p+)/k+^3","eq(3) coeff x 1/k+","match"))
for n,(a,b) in enumerate(zip(c2,c3),1):
    lhs = sp.simplify(pre2*a); rhs = sp.simplify(b/kps)
    print("   T%-5d %-22s %-22s %s"%(n, lhs, rhs, sp.simplify(lhs-rhs)==0))
print("   -> all six p+ weights transfer correctly; eq(3) bracket = eq(2) bracket x p+(k+-p+)/k+^2 .")
print()
print("   the delta^(2)[w - z + (p+/k+)(z-x)] sets  w = z - (p+/k+)(z-x) , so")
print("      -k.(w'-w) = -k.(w'-z) - (p+/k+) k.(z-x)   -> the two phases of eq (3).  CORRECT.")
print("   BUT eq (3) still writes Int_{w,w'} : the w integral has been used up by the delta.")

print()
print("="*78)
print("5.  LOCATING THE UV:  z -> x")
print("="*78)
def integrand(D, x, y, z, xp, wp, kp, pp):
    """the full eq(2) integrand, i and m contracted."""
    V = (xp-wp)/sq(xp-wp)
    K = (z-x)/sq(z-x)
    Br = Bracket_user(D, x, y, z, kp, pp)
    den = pp*sq(y-x) + (kp-pp)*sq(y-z)
    return float(V @ Br @ K)/den * (pp*(kp-pp)/kp**3)

def ang_avg(D, x, y, xp, wp, kp, pp, rad, n=4096):
    tot = 0.0
    for _ in range(n):
        d = rng.normal(size=D); d /= np.linalg.norm(d)
        tot += integrand(D, x, y, x + rad*d, xp, wp, kp, pp)
    return tot/n

def predicted(D, x, y, xp, wp, kp, pp):
    """u^2 x (angular average).  Derived below; note (x-z) = -u, and the 2/d_perp
    pieces of T1,T3,T5,T6 cancel identically."""
    V = (xp-wp)/sq(xp-wp); R = y-x; eta = pp/kp; eb = 1-eta
    Cuv = eta*eb + eta/eb + eb/eta          # = P_gg(eta)/(2 N_c)
    return (V@R)/(kp*sq(R)) * Cuv

print("   |u| -> 0 with x, y, x', w' fixed;  u^2 x (angular average of the integrand)")
print("   %-4s %-8s %12s %12s %12s %14s"%("D","eta","|u|=1e-2","1e-3","1e-4","predicted"))
for D in (2,3,4):
    for eta in (0.3, 0.65):
        x,y,xp,wp = rng.normal(size=(4,D))*1.0
        kp, pp = 1.0, eta
        vals = [rad**2*ang_avg(D,x,y,xp,wp,kp,pp,rad) for rad in (1e-2,1e-3,1e-4)]
        print("   %-4d %-8.2f %12.6f %12.6f %12.6f %14.6f"
              %(D,eta,vals[0],vals[1],vals[2],predicted(D,x,y,xp,wp,kp,pp)))
print()
print()
print("   => the integrand behaves as  C / |z-x|^2  with")
print("        C = [ (V.R)/(k+ R^2) ] x [ eta*etabar + eta/etabar + etabar/eta ] ,")
print("        eta = p+/k+ ,  R = y-x ,  V = (x'-w')/(x'-w')^2 ,")
print("      i.e. the UV residue is exactly  C_UV(eta) = P_gg(eta)/(2 N_c) .")
print("   d^2z = |u| d|u| dphi  =>  Int d|u|/|u| : a LOGARITHMIC UV DIVERGENCE at z -> x.")
print()
print("   The term-by-term singular parts (u = z-x, R = y-x, note (x-z) = -u):")
print("     T1  + d_perp delta^{im} (R.u)/u^2                    -> + (V.R)/(d_perp u^2) x d_perp")
print("     T2  + (k+/q)^2  R^i u^m /u^2                         -> + (k+/q)^2 (V.R)/u^2")
print("     T3  - 2 (k+/p+) delta_{im} (R.u)/u^2                 -> - 2(k+/p+)(V.R)/(d_perp u^2)")
print("     T4  + (k+/p+)^2 R^i u^m /u^2                         -> + (k+/p+)^2 (V.R)/u^2")
print("     T5  - 2 (k+/q) delta_{im} (R.u)/u^2                  -> - 2(k+/q)(V.R)/(d_perp u^2)")
print("     T6  + 2 (k+^2/(p+ q)) R^m u^i /u^2                   -> + 2(k+^2/(p+q))(V.R)/(d_perp u^2)")
print("   and the four 1/d_perp pieces cancel because  k+^2/(p+ q) - k+/p+ - k+/q = 0 .")
print("   What is left is  1 + (k+/q)^2 + (k+/p+)^2 , and the p+(k+-p+)/k+^2 weight of eq (3)")
print("   turns that into  eta*etabar + eta/etabar + etabar/eta  exactly.")

print()
print("="*78)
print("6.  ARE THERE ANY OTHER DANGEROUS REGIONS?")
print("="*78)
def scan(D, label, mover, fixed, kp=1.0, pp=0.4):
    print("   %-28s"%label, end="")
    for rad in (1e-1,1e-2,1e-3,1e-4):
        tot=0.0; n=2048
        for _ in range(n):
            d = rng.normal(size=D); d/=np.linalg.norm(d)
            tot += integrand(D, *mover(rad,d,*fixed), kp, pp)
        print(" %12.4e"%(rad**2*tot/n), end="")
    print()
D=2
x0,y0,xp0,wp0 = rng.normal(size=(4,D))
print("   u^2 x <integrand>  as the named separation -> 0   (a constant means log divergent)")
print("   %-28s %12s %12s %12s %12s"%("region","|u|=1e-1","1e-2","1e-3","1e-4"))
scan(D,"z -> x", lambda r,d: (x0, y0, x0+r*d, xp0, wp0), ())
scan(D,"y -> x  (z fixed)", lambda r,d: (x0, x0+r*d, x0+np.array([1.,0.]), xp0, wp0), ())
scan(D,"y -> z  (z fixed)", lambda r,d: (x0, x0+np.array([1.,0.])+r*d, x0+np.array([1.,0.]), xp0, wp0), ())
scan(D,"x' -> w'", lambda r,d: (x0, y0, x0+np.array([1.,0.]), wp0+r*d, wp0), ())
def joint(r,d):
    d2 = rng.normal(size=2); d2/=np.linalg.norm(d2)
    return (x0, x0+r*d2, x0+r*d, xp0, wp0)
print()
print("   Only  z -> x  gives a constant.  y -> x and y -> z give u^2 x <.> -> 0 (integrable);")
print("   so does x' -> w' , the LO Weizsacker-Williams kernel 1/|x'-w'| against |x'-w'| d|x'-w'|.")
print()
print("   LARGE |z| :  NOT convergent.  Two blocks GROW linearly in |z| :")
print("      T3's  (y-x).(y-z)/[2(y-x)^2]  ->  -(R.Z)/(2R^2)   and")
print("      T4's  (y-x)^m (y-z)^i/[2(y-x)^2] , T6's (y-x)^i (y-z)^m/[2(y-x)^2]  ->  ~ |Z| .")
print("   Against kernel ~1/|Z| and denominator ~ (k+-p+) Z^2 that leaves 1/Z^2 , and")
print("   d^2z ~ |z| d|z| gives  Int d|z|/|z| :  a SECOND, INFRARED logarithm.")
def largez(D,x,y,xp,wp,kp,pp,Z,n=20000):
    t=0.0
    for _ in range(n):
        d=rng.normal(size=D); d/=np.linalg.norm(d)
        t+=integrand(D,x,y,x+Z*d,xp,wp,kp,pp)
    return Z**2*t/n
def largez_pred(D,x,y,xp,wp,kp,pp):
    V=(xp-wp)/sq(xp-wp); R=y-x; q=kp-pp
    return (V@R)/(2*kp*sq(R))*( (kp/pp - 1.0)/D + kp/q )
print()
print("   %-4s %-7s %12s %12s %12s %14s"%("D","eta","|z|=1e3","1e4","1e5","predicted"))
for D in (2,3,4):
    for eta in (0.4,0.7):
        x1,y1,xp1,wp1 = rng.normal(size=(4,D))
        v=[largez(D,x1,y1,xp1,wp1,1.0,eta,Z) for Z in (1e3,1e4,1e5)]
        print("   %-4d %-7.2f %12.6f %12.6f %12.6f %14.6f"
              %(D,eta,v[0],v[1],v[2],largez_pred(D,x1,y1,xp1,wp1,1.0,eta)))
print()
print("   residue :  (V.R)/(2 k+ R^2) x [ (1/d_perp)(k+-p+)/p+ + k+/(k+-p+) ] .")
print("   Note it does NOT reduce to C_UV(eta) and its d_perp does NOT cancel -- a different")
print("   divergence from the UV one, with the colour structure U(z) -> 1 at large |z| .")
print()
print("   LARGE |y| :  T2, T4 and T6 all grow like |y| while the denominator grows like k+ y^2,")
print("   so the raw integrand falls only as 1/|y| .  The angular average is what matters, and")
print("   Monte-Carlo cannot resolve it here (the mean is far below the sampling noise), so use")
print("   exact trapezoidal quadrature in the angle -- spectrally accurate for a periodic integrand.")
def ang2d(fn, m=8192):
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    return float(np.mean([fn(np.array([np.cos(t),np.sin(t)])) for t in ph]))
x1,xp1,wp1=rng.normal(size=(3,2)); z1=x1+np.array([1.3,-0.7])
print()
print("   %-10s %16s %16s %16s"%("|y|","<integrand>","|y| x <.>","|y|^2 x <.>"))
for Y in (1e2,1e3,1e4,1e5,1e6):
    a_=ang2d(lambda d: integrand(2,x1,x1+Y*d,z1,xp1,wp1,1.0,0.4))
    print("   %-10.0e %16.6e %16.6e %16.6e"%(Y,a_,Y*a_,Y**2*a_))
print("   -> |y|^2 x <.> settles to a constant : the angular average falls as 1/y^2 (the 1/|y|")
print("      piece is odd in yhat and averages away), so Int d^2y ~ Int d|y|/|y| :")
print("      a THIRD logarithm, infrared, at large |y| .")
print()
print("   THE JOINT REGION z -> x AND y -> x is the SAME singularity, not a new one.")
print("     put y = x + rho1 a , z = x + rho2 b with rho1, rho2 ~ r :")
print("       kernel (z-x)/(z-x)^2 ~ 1/r , bracket ~ O(1) , denominator ~ r^2")
print("       => integrand ~ 1/r^3 , measure d^2y d^2z ~ r^4 dr/r  =>  r^0 dr/r : log")
print("     but that log is already the z -> x one: its residue (V.R)/(k+ R^2) behaves as 1/|R|")
print("     as R = y-x -> 0 , which is integrable against d^2y = |R| d|R| dphi .")
print("     So the y integral of the residue converges and there is NO log^2 from the transverse side.")
print()
print("   The p+ ENDPOINTS are separate: the residue eta*etabar + eta/etabar + etabar/eta")
print("   has 1/etabar and 1/eta poles at p+ -> k+ and p+ -> 0 , which is exactly why")
print("   eq (2) carries the cutoff  Int_Lambda^{k+-Lambda} dp+ .  Those are RAPIDITY")
print("   divergences, not transverse ones.")

print()
print("="*78)
print("7.  WHICH OF THE THREE LOGS NEEDS A REGULATOR")
print("="*78)
print("   z is the UNMEASURED GLUON's transverse position: it carries U(z) and is integrated")
print("     over all space.  Both of its logs are real:")
print("       z -> x        ULTRAVIOLET , residue C_UV(eta) = P_gg(eta)/(2 N_c) , d_perp-free")
print("       |z| -> oo     INFRARED    , residue (1/d_perp)(k+-p+)/p+ + k+/(k+-p+)")
print("                                   -- must cancel against the Bbar_2 (virtual) term,")
print("                                      where U(z) -> 1 : the same BK mechanism as Row IV.")
print()
print("   y and x' are SOURCE positions: they sit under rho^f(y) and rho^{b'}(x') .")
print("     The 1/y^2 tail is therefore controlled by the target correlator <rho rho>, not by")
print("     a regulator.  It is not a divergence of the loop integral.")
print()
print("   p+ -> 0 and p+ -> k+ are RAPIDITY divergences, cut by Lambda : Int_Lambda^{k+-Lambda}.")
print()
print("   ==> THE UV TERM IS  z -> x , and its coefficient is")
print()
print("       [1/(2pi)^3] (g^4/8 pi^4) (1/k+) Int dp+/2pi  P_gg(p+/k+)/(2 N_c)")
print("          x  (x'-w').(y-x) / [ (x'-w')^2 (y-x)^2 ]  x  Int_UV d^2z / (z-x)^2")
print("          x  (colour + sources) ,")
print()
print("   i.e. the UV residue is  P_gg  times the product of the two LO Weizsacker-Williams")
print("   kernels -- exactly the DGLAP form, and with no d_perp left over.")
