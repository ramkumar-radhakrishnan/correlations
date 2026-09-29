"""Part (b): UV subtraction in transverse coordinates.
Counterterm  S(r) = T_UV(r)/(r^2+s^2)   instead of   e^{iq.r} T_UV(r)/s^2 ."""
import numpy as np, sympy as sp, mpmath as mp
mp.mp.dps=20
bar="="*78
e=sp.Symbol('epsilon'); mu,s2=sp.symbols('mu s2',positive=True)
print(bar); print("1. SCALAR MASTER  J0 = mu^{-2eps} Int d^d r 1/[r^2 (r^2+s^2)] ,  d = 2-2eps"); print(bar)
u,a,c=sp.symbols('u a c',positive=True)
print("   radial: Int d^d r f(r^2) = pi^{d/2}/Gamma(d/2) Int_0^oo du u^{d/2-1} f(u)")
print("   Beta:   Int_0^oo du u^{a-1}/(u+c) = c^{a-1} Gamma(a) Gamma(1-a) ,  here a = d/2-1 = -eps")
J0 = sp.pi**(1-e)/sp.gamma(1-e)*sp.gamma(-e)*sp.gamma(1+e)*(mu**2*s2)**(-e)/s2
ser = sp.expand(sp.series(J0*s2/sp.pi, e, 0, 1).removeO())
print("   J0 * s^2/pi =", sp.simplify(ser))
print("   => J0 = (pi/s^2)[ -1/eps + gamma_E + log(pi mu^2 s^2) ]")
print("        = (pi/s^2)[ 1/eps_UV + log( mubar^2 s^2 e^{2 gamma_E}/4 ) ] ,  mubar^2 = 4 pi mu^2 e^{-gamma_E}")
mb=sp.Symbol('mubar',positive=True)
chk=sp.simplify(sp.EulerGamma+sp.log(sp.pi*mu**2*s2) - sp.log(mb**2*s2*sp.exp(2*sp.EulerGamma)/4)).subs(mb, sp.sqrt(4*sp.pi*mu**2*sp.exp(-sp.EulerGamma)))
print("   check of the MSbar rewriting:", sp.simplify(chk))
print("   cutoff |r|>rho : Int d^2r 1/[r^2(r^2+s^2)] = (pi/s^2) log(1+s^2/rho^2) -> (pi/s^2) log(s^2/rho^2)")

print(); print(bar); print("2. TENSOR PART IS EXACTLY ISOTROPIC"); print(bar)
print("   S depends on r only through |r| apart from the numerator, so")
print("   Int d^d r r^i r^m/[r^4(r^2+s^2)] = (delta^{im}/d) x (scalar J0)   EXACTLY, and with d_perp = d:")
print("   d_perp xi xibar (P^i s^m/P^2) delta^{im}/d  =  xi xibar P.s/P^2   -- no leftover finite term.")
print("   => mu^{-2eps} Int d^d r S(r) = pi W C_UV(xi) [ 1/eps_UV + log(mubar^2 (y-x)^2 e^{2gamma_E}/4) ]")

print(); print(bar); print("3. THE REMAINDER IS UV FINITE"); print(bar)
def sq(v): return float(v@v)
def T_of(r,s,P,xi):
    xb=1-xi; P2=sq(P); r2=sq(r); s2_=sq(s); smr=s-r; smr2=sq(smr); Pr=P@r; Ps=P@s; sr=s@r
    return (2*xi*xb*(Pr*sr/(P2*r2**2)-Pr/(2*P2*r2)) + (xi/xb)*(Ps/(P2*r2))*(1+(r@smr)/(2*smr2))
            + xb*(Pr/(P2*r2))*(2-2*sr/r2-sr/(2*s2_)) + (xb/xi)*((P@smr)/(P2*r2))*(1-sr/(2*s2_))
            + xi*(Pr/(P2*r2))*(0.5-2*sr/r2-(s@smr)/(2*smr2))
            + Pr*(2*sr-r2)/(P2*r2**2) + (P@smr)*sr/(2*P2*r2*smr2) - Ps*(r@smr)/(2*P2*r2*s2_))
def Tuv(r,s,P,xi):
    xb=1-xi; return 2*xi*xb*(P@r)*(s@r)/(sq(P)*sq(r)**2) + (xi/xb+xb/xi)*(P@s)/(sq(P)*sq(r))
s0=np.array([0.7,-0.4]); P0=np.array([-0.3,0.9]); k0=np.array([1.3,0.5])
def ang(fn,rho,m=4096):
    ph=(np.arange(m)+0.5)*2*np.pi/m
    return np.mean([fn(rho*np.array([np.cos(t),np.sin(t)])) for t in ph])
for xi in (0.3,0.7):
    xb=1-xi; M2=xi/xb*sq(s0); q=xb*k0
    F=lambda r: np.exp(1j*(q@r))*T_of(r,s0,P0,xi)/(xb*(sq(r-s0)+M2)) - Tuv(r,s0,P0,xi)/(sq(r)+sq(s0))
    print("   xi=%.1f  rho^2<remainder>:"%xi, "  ".join("%.1e"%abs(rho**2*ang(F,rho)) for rho in (1e-2,1e-3,1e-4)))
print("   -> falls as rho^2 : UV finite.  Large |r|: counterterm ~ 1/r^4 , the rest oscillates.")

print(); print(bar); print("4. EXACT RELATION TO THE MOMENTUM-SPACE SUBTRACTION (d = 2, no regulator needed)"); print(bar)
print("   Int d^2r T_UV [ e^{iq.r}/s^2 - 1/(r^2+s^2) ] : angular integrals")
print("     Int dphi e^{iq.r} = 2pi J0(q r) ;  Int dphi e^{iq.r} rhat^i rhat^m = pi delta J0 - pi(2 qhat qhat - delta) J2")
for qs in (0.3,1.0,3.0):
    I1=mp.quad(lambda t: (mp.besselj(0,qs*t)-1/(1+t**2))/t,[0,1,10])+mp.quadosc(lambda t:(mp.besselj(0,qs*t)-1/(1+t**2))/t,[10,mp.inf],period=2*mp.pi/qs)
    print("     Int_0^oo dx/x [J0(q s x) - 1/(1+x^2)] at qs=%.1f : %s   vs  -1/2 log(q^2 s^2 e^{2g}/4) = %s"
          %(qs, mp.nstr(I1,12), mp.nstr(-0.5*mp.log(qs**2*mp.e**(2*mp.euler)/4),12)))
print("     Int_0^oo dx J2(x)/x =", mp.nstr(mp.quadosc(lambda t: mp.besselj(2,t)/t,[0,mp.inf],period=2*mp.pi),12), " (= 1/2)")
print("   => [old counterterm] - [new counterterm] = -pi W C_UV log(q^2 s^2 e^{2gamma}/4) + pi xi xibar T")
print("      consistent: MSbar(old) - MSbar(new) = pi W C_UV [log(mu^2/q^2) - log(mu^2 s^2 e^{2g}/4)] + pi xi xibar T")

print(); print(bar); print("5. xi INTEGRALS"); print(bar)
print("   UV term: log has no xi  =>  Int C_UV d xi = 2L - 11/6 only.  No L^2 from the UV term.")
print("   remainder: Itilde_b = I_b + pi xi xibar T - pi W C_UV [ log xibar^2 + ell ] , ell = log(k^2 s^2 e^{2g}/4)")
print("     xi -> 0 :  I_b -> beta0/xi , extra -pi W ell / xi          => simple pole, beta0~ = beta0 - pi W ell")
print("     xi -> 1 :  I_b -> (-pi W log(1/xb) + beta1)/xb , extra (+2 pi W log(1/xb) - pi W ell)/xb")
print("               => (+pi W log(1/xb) + beta1 - pi W ell)/xb ,  beta1~ = beta1 - pi W ell")
print("   Int Itilde_b = (beta0~ + beta1~) L + (pi W/2) L^2 + Ftilde_b")
print("   consistency with the momentum-space result: the extra 2 W ell L from the UV log is cancelled by")
print("   (beta0~+beta1~-beta0-beta1)/pi L = -2 W ell L .  Same total, different split.")
