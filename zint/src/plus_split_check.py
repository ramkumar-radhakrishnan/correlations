"""Check the +-prescription split of the four-kernel row, and hunt for divergences
in the first (the "finite") term."""
import numpy as np, sympy as sp, mpmath as mp
mp.mp.dps = 25
Ei, g = mp.ei, mp.euler

print("A. THE SPLIT, TERM BY TERM")
pp, kp, P, L = sp.symbols('p k P Lam', positive=True)
print("   bracket middle term is  -delta_jm delta_ik' / p+ , and")
print("     -int_L^P dp+/p+ g  =  -int_L^P dp+ [g(p+)-g(0)]/p+  -  g(0) log(P/L)")
print("   with the overall minus in front of the row, the log piece comes out with a PLUS:")
print("     (-1) x (-1) x g(0) log(P/L)  =  + g(0) log(P/L)   <- matches your second term's sign")
print(f"   prefactor: (1/4pi^4) x (1/2pi) = 1/(8 pi^5)  ->  {sp.nsimplify(1/(4*sp.pi**4)/(2*sp.pi))}")
print("   so the 1/(8 pi^5) and the + sign of your second term are CORRECT.")
print()
print("   delta_jm delta_ik' contracted on (y'-z)^m (x'-y')^k' (y-w)^i (x-z)^j :")
print("      m=j : (y'-z).(x-z)     i=k' : (x'-y').(y-w)")
print("   which is exactly the numerator you wrote in the second term -> CORRECT.")
print()
print("   ==> BUT the second term still carries the whole bracket")
print("       [ d_k'm d_ij/(p++k+) - d_jm d_ik'/p+ - d_im d_jk'/k+ ].")
print("       That is a leftover: the contraction has ALREADY been done (it is in the")
print("       numerator), the 1/p+ has ALREADY been integrated (it is the log), and there")
print("       is no dp+ left for a free p+ to live in.  DELETE the bracket.")
print()

print("B. LARGE-|z| BEHAVIOUR OF THE THREE STRUCTURES  (rho^2 x angular average)")
rng = np.random.default_rng(3)
yp, xx, xpp, yy, ww = rng.normal(size=(5,2))
th = np.linspace(0, 2*np.pi, 16384, endpoint=False)
u  = np.stack([np.cos(th), np.sin(th)], -1)
a, b = xpp-yp, yy-ww                       # the spectator vectors
def zparts(z):
    A = yp - z; B = xx - z
    S1 = (A@a if A.ndim==1 else A@a)/np.sum(A*A,-1) * (B@b)/np.sum(B*B,-1)
    S2 = np.sum(A*B,-1)/(np.sum(A*A,-1)*np.sum(B*B,-1))
    S3 = (A@b)/np.sum(A*A,-1) * (B@a)/np.sum(B*B,-1)
    return S1, S2, S3
print(f"   {'rho':>8} {'S1 (d_k'+chr(39)+'m d_ij)':>22} {'S2 (d_jm d_ik'+chr(39)+')':>22} {'S3 (d_im d_jk'+chr(39)+')':>22}")
for rho in (1e2, 1e3, 1e4):
    S1,S2,S3 = zparts(u*rho)
    print(f"   {rho:8.0e} {np.mean(S1)*rho**2:22.6f} {np.mean(S2)*rho**2:22.6f} {np.mean(S3)*rho**2:22.6f}")
print(f"   predicted:  S1 -> (x'-y').(y-w)/2 = {a@b/2:.6f}   S2 -> 1   S3 -> (y-w).(x'-y')/2 = {a@b/2:.6f}")
print("   ALL THREE fall only as 1/rho^2, i.e. every one is log divergent WITHOUT a phase.")
print()

print("C. WHICH ONES ARE RESCUED BY THE PHASE?   phase = e^{-i xi k.(y'-z)}")
kk = np.array([0.9, 0.6])
def avg_with(phasefun, S, rho):
    z = u*rho
    return np.mean(phasefun(z)*S)
for xi in (0.3, 0.05):
    print(f"   xi = {xi}")
    for rho in (1e2, 1e3, 1e4):
        z = u*rho
        ph = np.exp(-1j*xi*((yp-z)@kk))
        S1,S2,S3 = zparts(z)
        m1 = np.mean(ph*S1)*rho**2; m2 = np.mean(ph*S2)*rho**2; m3 = np.mean(ph*S3)*rho**2
        m2m = np.mean((ph-1)*S2)*rho**2
        print(f"      rho={rho:.0e}:  |S1 x ph| {abs(m1):9.2e}  |S2 x ph| {abs(m2):9.2e}"
              f"  |S3 x ph| {abs(m3):9.2e}   |S2 x (ph-1)| {abs(m2m):9.6f}")
print("   with the full phase all three averages decay (the oscillation kills them):")
print("   the z-integral is finite for every xi > 0.  But the SUBTRACTED combination")
print("   (phase - 1) x S2 keeps the phase-less -1, whose average stays at -1/rho^2.")
print()

print("D. SO: DOES THE FIRST TERM DIVERGE?")
print("   [1/p+]_+ acts as  [g(p+) - g(0)]/p+ .  g(0) is the integrand with phase = 1.")
print("   Its z-integrand is S2 = (y'-z).(x-z)/[(y'-z)^2(x-z)^2], the BK REAL kernel,")
print("   with NO phase left to cut it.  int rho drho (-1/rho^2) = -log : DIVERGENT,")
print("   for every value of p+.  The same divergence sits in your second term")
print("   (it IS g(0)).  They cancel only in the sum, i.e. only if you never split.")
print()

print("E. A SECOND, SEPARATE PROBLEM: THE UPPER LIMIT")
print("   int_0^Xi dxi [e^{-i kap xi} - 1]/xi = Ei(-i kap Xi) - gamma - log(i kap Xi)")
print("   -> -log(Xi) as Xi -> infinity.  The unsplit integral has a finite V -> infinity limit,")
print("      so the split introduces a spurious log(V) into EACH piece:")
kap = mp.mpf('1.7'); lam = mp.mpf('1e-9')
for X in ['1e2','1e4','1e6']:
    Xi = mp.mpf(X)
    unsplit = Ei(-1j*kap*Xi) - Ei(-1j*kap*lam)
    plus    = Ei(-1j*kap*Xi) - g - mp.log(1j*kap*Xi)
    logpc   = mp.log(Xi/lam)
    print(f"   Xi={X}:  unsplit {mp.nstr(unsplit,10):>28}   +piece {mp.nstr(plus,10):>28}"
          f"   log(Xi/lam) {float(logpc):9.4f}")
print("   the unsplit column is V-stable; the + piece slides by -log(Xi) and the log piece by +log(Xi).")
print()
print("F. THE FIX FOR (E): subtract only below k+  (theta(k+ - p+) in the subtraction)")
print("   [1/p+]_+^(k+) g = [ g(p+) - theta(k+-p+) g(0) ] / p+ ,  isolating l_k = log(k+/Lambda)")
for X in ['1e2','1e4','1e6']:
    Xi = mp.mpf(X)
    plus_k = (Ei(-1j*kap*Xi) - Ei(-1j*kap*1)) + (Ei(-1j*kap*1) - g - mp.log(1j*kap))
    print(f"   Xi={X}:  + piece with theta cut = {mp.nstr(plus_k,10):>28}    (l_k = log(1/lam) = {float(mp.log(1/lam)):.4f})")
print("   now the + piece has a finite V -> infinity limit and the log is l_k.")
