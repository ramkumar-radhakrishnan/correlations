"""Complete part 2b: extract the UV at z -> x, regularise it, and do the p+ integral."""
import sympy as sp, mpmath as mp, numpy as np
mp.mp.dps = 25
eta, lam = sp.symbols('eta lamb', positive=True)
eb = 1 - eta
C = eta*eb + eta/eb + eb/eta                 # = P_gg(eta)/(2 N_c)

print("="*78)
print("1.  THE MASTER z-INTEGRAL: the large-|z| end is cut by the PHASE, not divergent")
print("="*78)
print("   In 2b the colour factor N_c[U(y)-U(x)] is z-independent, so the only z dependence")
print("   left is the kernel and the phase.  Writing the total phase as")
print("        exp[-i k.(w' - etabar z - eta x)]  =  e^{-i k.w'} e^{i etabar k.z} e^{i eta k.x} ,")
print("   the z integral of the UV residue is")
print("        Int d^2z  e^{i q.(z-x)} / (z-x)^2 ,        q = etabar k ,")
print("   which converges at large |z| by oscillation and is log divergent only at z -> x .")
print("   With a short-distance cutoff |z-x| > rho :  = -pi log[ q^2 rho^2 e^{2 gamma_E}/4 ] .")
def cutoff(q, rho, T=60.0):
    q = mp.mpf(q); rho = mp.mpf(rho)
    a = mp.quad(lambda r: mp.besselj(0,q*r)/r, [rho, 1/q, T/q])
    b = mp.quadosc(lambda r: mp.besselj(0,q*r)/r, [T/q, mp.inf], period=2*mp.pi/q)
    return 2*mp.pi*(a+b)
print()
print("   %-10s %-6s %20s %24s"%("rho","q","numerical","-pi log[q^2 rho^2 e^{2g}/4]"))
for rho,q in ((1e-4,1.3),(1e-6,1.3),(1e-6,0.4)):
    num = cutoff(q,rho); ana = -mp.pi*mp.log(mp.mpf(q)**2*mp.mpf(rho)**2*mp.e**(2*mp.euler)/4)
    print("   %-10.0e %-6.2f %20s %24s"%(rho,q,mp.nstr(num,10),mp.nstr(ana,10)))
print("   -> in d = 2-2eps with MSbar:  -pi [ 1/eps + log(q^2/mubar^2) ] , mubar^2 = 4 pi mu^2 e^{-gamma_E}")
print("      The pole is purely ULTRAVIOLET (the phase already handled large |z|).")

print()
print("="*78)
print("2.  THE p+ INTEGRAL OF THE UV RESIDUE:  two moments of C_UV")
print("="*78)
I0 = sp.simplify(sp.integrate(C, (eta, lam, 1-lam)))
print("   I0 = Int_lam^{1-lam} C_UV d eta = -2 log(lam) - 11/6 + O(lam)")
print("        sympy:", sp.re(sp.expand(sp.series(I0, lam, 0, 1).removeO())))
print("        decomposition:  Int eta etabar = 1/6 ,  each of Int eta/etabar and Int etabar/eta")
print("                        = log(1/lam) - 1  =>  2 log(1/lam) - 2 + 1/6 = 2 log(1/lam) - 11/6")
print()
parts = {'eta*etabar': eta*eb, 'eta/etabar': eta/eb, 'etabar/eta': eb/eta}
tot = 0
for nm, tm in parts.items():
    if nm == 'etabar/eta':
        # the lam -> 0 limit is finite; do it in closed form
        v = 2*sp.integrate(eb/eta*sp.log(eb), (eta, 0, 1))
    else:
        v = sp.series(sp.simplify(sp.integrate(tm*sp.log(eb**2), (eta, lam, 1-lam))), lam, 0, 1).removeO()
    v = sp.simplify(sp.expand(v))
    print("   Int %-12s log(etabar^2) d eta = %s"%(nm, v))
    tot += v
I1 = sp.simplify(sp.expand(tot))
print("   I1 = ", I1)
print()
L = sp.Symbol('L', positive=True)     # L = log(k+/Lambda) = -log(lam)
I0L = sp.simplify(sp.re(sp.expand(sp.series(I0, lam, 0, 1).removeO())).subs(sp.log(lam), -L))
I1L = sp.simplify(I1.subs(sp.log(lam), -L))
print("   in terms of L = log(k+/Lambda) :")
print("       I0 =", I0L)
print("       I1 =", sp.expand(I1L))

print()
print("="*78)
print("3.  NUMERICAL CHECK OF I0 AND I1")
print("="*78)
Cf = lambda e: e*(1-e) + e/(1-e) + (1-e)/e
print("   %-10s %18s %18s %18s %18s"%("lam","I0 quad","I0 formula","I1 quad","I1 formula"))
for Lm in (1e-3,1e-5,1e-8):
    q0 = mp.quad(Cf, [Lm, 0.5, 1-Lm])
    q1 = mp.quad(lambda e: Cf(e)*mp.log((1-e)**2), [Lm, 0.5, 1-Lm])
    a0 = -2*mp.log(Lm) - mp.mpf(11)/6
    a1 = -mp.log(Lm)**2 + mp.mpf(67)/18 - mp.pi**2/3
    print("   %-10.0e %18.8f %18.8f %18.8f %18.8f"%(Lm,q0,a0,q1,a1))

print()
print("="*78)
print("4.  THE SUBTRACTED REMAINDER IS UV CONVERGENT")
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
def full(x,y,z,xp,wp,kvec,kp,pp):
    """the braces of (1.16) including the 1/k+ and the total phase."""
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x); e=pp/kp
    ph=np.exp(-1j*(kvec@(wp-(1-e)*z-e*x)))
    return ph*float(V@Bracket(2,x,y,z,kp,pp)@K)/(pp*sq(y-x)+(kp-pp)*sq(y-z))/kp
def subtr(x,y,z,xp,wp,kvec,kp,pp):
    """the UV counterterm: C_UV(eta) x (LO WW product) x e^{i q.(z-x)}/(z-x)^2 x total phase at z=x."""
    V=(xp-wp)/sq(xp-wp); R=y-x; e=pp/kp
    Cuv=e*(1-e)+e/(1-e)+(1-e)/e
    ph=np.exp(-1j*(kvec@(wp-x)))                       # the z -> x phase
    osc=np.exp(1j*(1-e)*(kvec@(z-x)))                  # keeps the large-|z| oscillation
    return ph*osc*Cuv*(V@R)/(kp*sq(R))/sq(z-x)
x0=np.array([0.30,-0.20]); y0=np.array([-0.7,1.1]); xp0=np.array([1.6,0.4]); wp0=np.array([-0.5,-1.3])
kvec=np.array([0.5,-0.3]); kp=1.0; pp=0.4
def ang(fn,r,m=4096):
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    return np.mean([fn(r*np.array([np.cos(a),np.sin(a)])) for a in ph])
print("   u^2 x < integrand >  as u = z-x -> 0   (constant = log divergent, -> 0 = finite)")
print("   %-10s %22s %22s %22s"%("|z-x|","full (1.16) braces","counterterm","difference"))
for r in (1e-2,1e-3,1e-4,1e-5):
    F=r**2*ang(lambda d: full(x0,y0,x0+d,xp0,wp0,kvec,kp,pp),r)
    S=r**2*ang(lambda d: subtr(x0,y0,x0+d,xp0,wp0,kvec,kp,pp),r)
    D=r**2*ang(lambda d: full(x0,y0,x0+d,xp0,wp0,kvec,kp,pp)-subtr(x0,y0,x0+d,xp0,wp0,kvec,kp,pp),r)
    print("   %-10.0e %22s %22s %22s"%(r,np.format_float_scientific(F.real,5),
          np.format_float_scientific(S.real,5), np.format_float_scientific(D.real,5)))
print("   -> the counterterm reproduces the full UV residue; the difference vanishes, so")
print("      [ (1.16) - counterterm ] can be integrated at d_perp = 2 with no regulator.")

print()
print("="*78)
print("5.  THE COMPLETED PART 2b")
print("="*78)
print("""
   d3N_2b/d3k  =  UV  +  finite , with

   UV = - [1/(2pi)^3] (g^4 / 16 pi^4)
        { [ 2L - 11/6 ] [ 1/eps + log(k_perp^2/mubar^2) ] - L^2 + 67/18 - pi^2/3 }
        x Int_{x,y,x',w'} e^{-i k.(w'-x)} (x'-w').(y-x) / [(x'-w')^2 (y-x)^2]
        x [U^{ab'}(x') - U^{ab'}(w')] N_c [U^{af}(y) - U^{af}(x)] rho^{b'}(x') rho^f(y)

   finite = [1/(2pi)^3] (g^4/8 pi^4) (1/k+) Int_Lambda^{k+-Lambda} dp+/2pi Int_{x,y,z,x'}
            { (1.16) braces  -  C_UV(eta) e^{i etabar k.(z-x)}/(z-x)^2
                                x (x'-w').(y-x)/[(x'-w')^2 (y-x)^2] / k+ }
            x (same phase, same colour, d_perp = 2)

   where L = log(k+/Lambda) , mubar^2 = 4 pi mu^2 e^{-gamma_E} ,
   and the overall (-pi) x (1/2pi) x (g^4/8pi^4) = -(g^4/16 pi^4) with the k+ cancelling.
""")
