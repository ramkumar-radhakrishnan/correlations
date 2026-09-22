"""Group II Row IV: check the change of variable p+ -> xi = k+/(k+ + p+) block by block."""
import numpy as np, mpmath as mp
mp.mp.dps=20
K=1.0; Lam=1e-3; Vee=40.0; kap=0.9; dperp=2.0
XX, Xr, M = 1.7, -0.6, 3.3           # stand-ins for (X.X') , (X.rhat)(X'.rhat) , M
def S(p): return K+p
# --- the p+ form (the verified three blocks; P = 1/(2pi)^3 g^4/4pi^3 Nc/k+ factored out) ---
ph = lambda p: mp.e**(-1j*kap*S(p)/K)
b1 = lambda p: (p/S(p)**2 + p/K**2 + 1/p)*M*XX
b2 = lambda p: dperp*p/S(p)**2*Xr
b3 = lambda p: (p/K**2 + 1/p)*XX
I1 = mp.quad(lambda p: ph(p)*b1(p),[Lam,1,Vee-K])
I2 = mp.quad(lambda p: ph(p)*b2(p),[Lam,1,Vee-K])
I3 = mp.quad(lambda p: ph(p)*b3(p),[Lam,1,Vee-K])
tot_p = I1 - I2 - I3                       # signs of the three blocks: + , - , -
# --- the xi form ---
xL = K/(K+Lam); xV = K/Vee                  # xi at p+ = Lambda  and at p+ = V - k+
phx = lambda x: mp.e**(-1j*kap/x)
C   = lambda x: x*(1-x) + x/(1-x) + (1-x)/x
# correct: each block picks up  dp+ = -(k+/xi^2) dxi , limits xL -> xV
J1 = -mp.quad(lambda x: phx(x)/x**2*C(x)*M*XX,[xL,0.5,xV])
J2 = +mp.quad(lambda x: phx(x)/x**2*dperp*x*(1-x)*Xr,[xL,0.5,xV])
J3 = +mp.quad(lambda x: phx(x)/x**2*(x/(1-x)+(1-x)/x)*XX,[xL,0.5,xV])
print("="*74); print("BLOCK-BY-BLOCK: p+ form vs xi form"); print("="*74)
print("   block 1 :  p+ %s   xi %s"%(mp.nstr(+I1,10), mp.nstr(J1,10)))
print("   block 2 :  p+ %s   xi %s"%(mp.nstr(-I2,10), mp.nstr(J2,10)))
print("   block 3 :  p+ %s   xi %s"%(mp.nstr(-I3,10), mp.nstr(J3,10)))
print("   total   :  p+ %s   xi %s"%(mp.nstr(tot_p,10), mp.nstr(J1+J2+J3,10)))
print()
print("   => with the limits written  xi : k+/(k++Lambda) -> k+/V  (descending),")
print("      block 1 carries an overall MINUS and blocks 2,3 carry a PLUS,")
print("      and ALL THREE carry the 1/xi^2 of the Jacobian.")
print()
yours = -mp.quad(lambda x: phx(x)/x**2*C(x)*M*XX,[xL,0.5,xV]) \
        -mp.quad(lambda x: phx(x)*dperp*x*(1-x)*Xr,[xL,0.5,xV]) \
        -mp.quad(lambda x: phx(x)*(x/(1-x)+(1-x)/x)*XX,[xL,0.5,xV])
print("   what you wrote (all three minus, no 1/xi^2 in blocks 2,3):")
print("      %s     vs correct  %s"%(mp.nstr(yours,10), mp.nstr(J1+J2+J3,10)))
