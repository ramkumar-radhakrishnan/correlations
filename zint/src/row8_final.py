"""Row 8: assemble the eps-expanded z-integral, then do the p+ integrals."""
import sympy as sp, numpy as np, mpmath as mp

eps = sp.symbols('epsilon')
P,Kp = sp.symbols('pplus kplus', positive=True)
X,R  = sp.symbols('X R')        # X = (x-y).(x'-y')/[..], R = rhat.(x-y) rhat.(x'-y')/[..]
L    = sp.symbols('Lbar')       # 1/eps - log(r^2 mubar^2)
S    = Kp+P

Pcal = 1/P + P/Kp**2 + P/S**2
Fcal = P*X/S**2 + 2*P*R/S**2 + 2*(P/Kp**2 + 1/P)*X
tot  = sp.expand((L+2)*Pcal*X - Fcal)

print("="*74); print("1.  eps-expanded z-integral, divided by pi"); print("="*74)
print("   total/pi =", sp.simplify(tot))
c1 = sp.simplify(sp.expand(tot*P).subs(P,0))                     # coefficient of 1/p+
c2 = sp.simplify(sp.diff(sp.expand(tot - c1/P - sp.together(0)),P).subs(P,0))
rest = sp.simplify(sp.expand(tot - c1/P))
print()
print("   coefficient of 1/p+            :", sp.factor(c1), "   <- the +2 cancelled")
print("   remainder                      :", sp.factor(rest))
print("   ... regrouped:")
A2 = sp.simplify(sp.expand(rest*Kp**2/P).subs(1/S**2,0))
print("       coeff of p+/k+^2           :", sp.factor(sp.simplify(L*X)), "  (check below)")
targ = c1/P + L*X*P/Kp**2 + ((L+1)*X - 2*R)*P/S**2
print("   PROPOSED closed form total/pi  :  L*X/p+ + L*X p+/k+^2 + [(L+1)X - 2R] p+/(k+ + p+)^2")
print("   difference from the direct expansion :", sp.simplify(sp.expand(tot - targ)))

print(); print("="*74); print("2.  THE THREE p+ INTEGRALS  (kappa = k.(y'-y)/k+, Pv = V - k+)"); print("="*74)
p,kap,Pv,Lam = sp.symbols('p kappa P_vee Lambda', positive=True)
I2 = sp.simplify(sp.integrate(p*sp.exp(-sp.I*kap*p),(p,0,Pv)))
print("   I2 = int_0^Pv dp p e^{-i kappa p}      =", sp.simplify(sp.factor(I2)))
q,V = sp.symbols('q V', positive=True)
I3 = sp.exp(sp.I*kap*Kp)*sp.integrate((q-Kp)/q**2*sp.exp(-sp.I*kap*q),(q,Kp,V))
print("   I3 = int_0^Pv dp p e^{-i kappa p}/(k+ + p)^2 :")
print("        = e^{i kappa k+} [ (1 + i kappa k+)(Ei(-i kappa V) - Ei(-i kappa k+))")
print("                           + k+ ( e^{-i kappa V}/V - e^{-i kappa k+}/k+ ) ]")
# numerical verification of I3 closed form
def I3_num(kk, kp, Vv):
    f = lambda t: t*mp.e**(-1j*kk*t)/(kp+t)**2
    return mp.quad(f, [0, Vv-kp])
def I3_cf(kk, kp, Vv):
    Ei = lambda u: mp.ei(u)
    return mp.e**(1j*kk*kp)*((1+1j*kk*kp)*(Ei(-1j*kk*Vv)-Ei(-1j*kk*kp))
                             + kp*(mp.e**(-1j*kk*Vv)/Vv - mp.e**(-1j*kk*kp)/kp))
for kk,kp,Vv in [(0.7,1.0,20.0),(2.3,0.4,50.0),(0.13,2.0,9.0)]:
    a,b = I3_num(kk,kp,Vv), I3_cf(kk,kp,Vv)
    print("      kappa=%-5.2f k+=%-4.1f V=%-5.1f  quad=%s  closed=%s  diff=%.2e"
          % (kk,kp,Vv,mp.nstr(a,8),mp.nstr(b,8),abs(a-b)))
print()
print("   I1 = int_Lam^Pv dp e^{-i kappa p}/p  (PLUS PRESCRIPTION)")
print("      = log(Pv/Lam) + int_0^Pv dp/p [e^{-i kappa p} - 1]")
print("      = log(Pv/Lam) - Cin(|kappa| Pv) - i sgn(kappa) Si(|kappa| Pv)     + O(kappa Lam)")
for kk,Pvv,Lm in [(0.7,20.0,1e-6),(2.3,50.0,1e-7)]:
    num = mp.quad(lambda t: mp.e**(-1j*kk*t)/t, [Lm, Pvv])
    cf  = mp.log(Pvv/Lm) - (mp.euler+mp.log(kk*Pvv)-mp.ci(kk*Pvv)) - 1j*mp.si(kk*Pvv)
    print("      kappa=%-5.2f Pv=%-5.1f Lam=%.0e  quad=%s closed=%s diff=%.2e"
          % (kk,Pvv,Lm,mp.nstr(num,8),mp.nstr(cf,8),abs(num-cf)))
