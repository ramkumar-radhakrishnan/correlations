"""Clean limits Lambda -> 0, V -> infinity, and what multiplies the rapidity log."""
import mpmath as mp, numpy as np
mp.mp.dps = 30
Ei, g = mp.ei, mp.euler

print("A. Ei at large imaginary argument is a CONSTANT, not zero:")
for X in [1e2, 1e4, 1e6]:
    print(f"   Ei(-i*{X:.0e}) = {mp.nstr(Ei(-1j*mp.mpf(X)),10)}      -i*pi = {mp.nstr(-1j*mp.pi,10)}")
for X in [1e2, 1e4, 1e6]:
    print(f"   Ei(+i*{X:.0e}) = {mp.nstr(Ei(1j*mp.mpf(X)),10)}      +i*pi = {mp.nstr(1j*mp.pi,10)}")
print("   -> Ei(-i kap Xi) -> -i pi sgn(kap)  as Xi -> infinity.")
print()

print("B. the -1/xi bracket in the double limit  lam -> 0, Xi -> infinity")
print("   claim:  -[Ei(-i kap Xi) - Ei(-i kap lam)]  ->  gamma + log( i kap lam )")
for kap in [1.7, -1.7, 0.3, -4.2]:
    k = mp.mpf(kap); lam = mp.mpf('1e-9'); Xi = mp.mpf('1e9')
    ex = -(Ei(-1j*k*Xi) - Ei(-1j*k*lam))
    cl = g + mp.log(1j*k*lam)
    print(f"   kap={kap:+5.1f}:  exact {mp.nstr(ex,12):>34}   gamma+log(i kap lam) {mp.nstr(cl,12):>34}"
          f"   diff {mp.nstr(abs(ex-cl),3)}")
print("   with lam = Lambda/k+ and l_k = log(k+/Lambda):")
print("      =  -l_k + gamma_E + log( i k.(y'-z) )        <- THE rapidity log, coefficient -1")
print()

print("C. the other two brackets in the same limit: finite, no l_k")
for kap in [1.7, -1.7]:
    k = mp.mpf(kap); lam = mp.mpf('1e-9'); Xi = mp.mpf('1e9')
    A = mp.e**(1j*k)*(Ei(-1j*k*(1+Xi)) - Ei(-1j*k*(1+lam)))
    Acl = mp.e**(1j*k)*(-1j*mp.pi*mp.sign(k) - Ei(-1j*k))
    C = (mp.e**(-1j*k*Xi) - mp.e**(-1j*k*lam))/(1j*k)
    Ccl = -1/(1j*k)
    print(f"   kap={kap:+5.1f}  1/(1+xi): exact {mp.nstr(A,10):>30}  limit {mp.nstr(Acl,10):>30}")
    print(f"             -1     : exact {mp.nstr(C,10):>30}  limit {mp.nstr(Ccl,10):>30}  (+ oscillatory)")
print()

print("D. WHAT MULTIPLIES l_k  --  it is the BK/JIMWLK dipole kernel")
print("   the -1/xi structure is  (y'-z).(x-z) / [ (y'-z)^2 (x-z)^2 ]  x  (x'-y').(y-w)/[...]")
print("   the first factor is exactly  K(x,y';z) = (x-z).(y'-z)/[(x-z)^2 (y'-z)^2],")
print("   the real-emission piece of the BK kernel, with z the emission point.")
print()
print("E. the transverse IR log CANCELS once the virtual partners are added.")
print("   BK: M(x,y';z) = (x-y')^2/[(x-z)^2 (y'-z)^2]")
print("                 = 1/(x-z)^2 + 1/(y'-z)^2 - 2 (x-z).(y'-z)/[(x-z)^2 (y'-z)^2]")
x = np.array([0.3, -0.2]); yp = np.array([-0.5, 0.8])
th = np.linspace(0, 2*np.pi, 8192, endpoint=False)
print("   rho^2 x <angular average> at large rho :")
print(f"   {'rho':>8} {'1/(x-z)^2':>12} {'1/(y-z)^2':>12} {'-2 K':>12} {'sum = M':>12}")
for rho in (1e1, 1e2, 1e3, 1e4):
    z = np.stack([np.cos(th), np.sin(th)], -1)*rho
    a = 1/np.sum((x-z)**2, -1); b = 1/np.sum((yp-z)**2, -1)
    K = np.sum((x-z)*(yp-z), -1)*a*b
    print(f"   {rho:8.0e} {np.mean(a)*rho**2:12.6f} {np.mean(b)*rho**2:12.6f} "
          f"{np.mean(-2*K)*rho**2:12.6f} {np.mean(a+b-2*K)*rho**2:12.6f}")
print("   -> the -2K piece alone gives -2/rho^2 (log divergent, this is the log I flagged);")
print("      the full BK combination gives 0 x 1/rho^2, i.e. it falls as 1/rho^4: CONVERGENT.")
