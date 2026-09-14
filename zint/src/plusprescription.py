"""(1) the w' delta-function step, (2) whether the + prescription is legitimate here."""
import numpy as np, sympy as sp, mpmath as mp
mp.mp.dps = 30
Ei, g = mp.ei, mp.euler

print("A. THE w' STEP  (checks the 1/k+ bookkeeping)")
pp, kp = sp.symbols('p k', positive=True)                 # p+, k+
yp, z = sp.symbols("y' z")                                # 1-d stand-ins for the transverse points
wsol = ((pp+kp)*yp - pp*z)/kp
print(f"   w' = {sp.simplify(wsol)}")
print(f"   w' - z = {sp.factor(sp.simplify(wsol - z))}        (= (p+ + k+)/k+ * (y'-z))")
print(f"   so  (w'-z)^m/(w'-z)^2  ->  k+/(p+ + k+) * (y'-z)^m/(y'-z)^2")
print(f"   w' - y' = {sp.factor(sp.simplify(wsol - yp))}       (phase: e^-ik(w'-w) = e^-ik(y'-w) e^-i(p+/k+)k(y'-z))")
pref2 = 1/kp**2
brack2 = sp.Matrix([1, -(pp+kp)/pp, -(pp+kp)/kp])          # eq (2) bracket coefficients
brack3 = sp.Matrix([1/(pp+kp), -1/pp, -1/kp])              # eq (3) bracket coefficients
lhs = sp.simplify(pref2*kp/(pp+kp)*brack2 - (1/kp)*brack3)
print(f"   (1/k+^2)(k+/(p++k+)) x [eq2 bracket]  -  (1/k+) x [eq3 bracket]  =  {list(lhs)}")
print("   -> the w' integration and the 1/k+ out front of eq (3) are CORRECT.")
print()

print("B. CORRECTION TO WHAT I SAID LAST TIME ABOUT THE OVERALL 1/k+")
xi = sp.symbols('xi', positive=True)
for nm, term in [("1/(p+ + k+)", 1/(pp+kp)), ("1/p+", 1/pp), ("1/k+", 1/kp)]:
    sub = term.subs(pp, kp*xi)*kp                          # dp+ = k+ dxi
    print(f"   int dp+ [ {nm:12s} ] = int dxi [ {sp.simplify(sub)} ]")
print("   every term carries exactly ONE inverse power of momentum, so the k+ of")
print("   dp+ = k+ dxi cancels INSIDE the bracket.  int dp+ [bracket] = int dxi [A/(1+xi) - B/xi - C]")
print("   and the 1/k+ out front SURVIVES.  My earlier 'used twice' flag was wrong.")
print()

print("C. IS THE + PRESCRIPTION CORRECT?   identity vs the Ei form")
print("   int_lam^Xi dxi e^{-i kap xi}/xi  =  int_0^Xi dxi [e^{-i kap xi} - 1]/xi  +  log(Xi/lam)")
for kap in [1.7, -0.4, 5.0]:
    for lam in ['1e-6','1e-9']:
        k_, l_, X_ = mp.mpf(kap), mp.mpf(lam), mp.mpf('37.0')
        exact = Ei(-1j*k_*X_) - Ei(-1j*k_*l_)
        plus  = mp.quad(lambda u: (mp.e**(-1j*k_*u)-1)/u, [0, 1/abs(k_), X_]) + mp.log(X_/l_)
        print(f"   kap={kap:+5.1f} lam={lam}:  Ei form {mp.nstr(exact,12):>32}   + form {mp.nstr(plus,12):>32}"
              f"   diff {mp.nstr(abs(exact-plus),3)}")
print("   identical.  The subtracted piece is Ei(-i kap Xi) - gamma - log(i kap Xi), finite.")
print()

print("D. BUT: IS THE SUBTRACTION TERM g(0) TRANSVERSELY FINITE?")
print("   At p+ = 0 the phase e^{-i(p+/k+)k.(y'-z)} -> 1.  The delta_jm delta_ik' structure")
print("   then leaves  (y'-z).(x-z)/[(y'-z)^2 (x-z)^2]  in the z-integration -- the BK REAL kernel.")
yy = np.array([-0.5, 0.8]); xx = np.array([0.3, -0.2])
th = np.linspace(0, 2*np.pi, 8192, endpoint=False)
u = np.stack([np.cos(th), np.sin(th)], -1)
print(f"   {'rho':>8} {'rho^2 <K_real>':>16} {'rho^2 <K_BK full>':>20}")
for rho in (1e1, 1e2, 1e3, 1e4):
    zz = u*rho
    a = 1/np.sum((xx-zz)**2, -1); b = 1/np.sum((yy-zz)**2, -1)
    K = np.sum((xx-zz)*(yy-zz), -1)*a*b
    print(f"   {rho:8.0e} {np.mean(-2*K)*rho**2:16.6f} {np.mean(a+b-2*K)*rho**2:20.6f}")
print("   -> -2K alone tends to -2/rho^2 : int rho drho/rho^2 is LOG DIVERGENT.")
print("      the full BK kernel M = (x-y')^2/[(x-z)^2 (y'-z)^2] tends to 0 x 1/rho^2 (falls as 1/rho^4).")
print()

print("E. WHICH LOG DOES THE + PRESCRIPTION ISOLATE?")
print("   log(Xi/lam) = log((V-k+)/Lambda),  NOT l_k = log(k+/Lambda).")
print("   difference log((V-k+)/k+) is a finite scheme choice; it shifts what you call")
print("   the impact factor.  Numbers for V = 100 k+:")
for lam in (1e-6, 1e-9):
    lk = np.log(1/lam); X = 99.0
    print(f"   Lambda/k+={lam:.0e}:  log(Xi/lam) = {np.log(X/lam):10.5f}   l_k = {lk:10.5f}"
          f"   difference = {np.log(X):8.5f}")
