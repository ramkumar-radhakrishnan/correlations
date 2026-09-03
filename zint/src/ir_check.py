"""Where (if anywhere) does the four-kernel row diverge in the TRANSVERSE integrals?

  master:  B^{mj}(Q,c) = int d^2r e^{-iQ.r} (r^m/r^2)((c-r)^j/(c-r)^2)
  z-integral : r = z-y',  c = x-y',  Q = -xi k
  y'-integral: r = y'-z,  c = x'-z,  Q = (1+xi) k
"""
import numpy as np
from scipy.integrate import quad
from nlo2_master import Bten
g = 0.5772156649015329
c = np.array([0.7, -0.4]); cn = np.hypot(*c)

print("1. SHORT DISTANCE.  angular average of the integrand near each singular point")
print("   (a log needs a nonzero angular average of the 1/rho^2 piece)")
def ang(rho, about):          # about=0 -> r->0 ; about=1 -> r->c
    th = np.linspace(0, 2*np.pi, 4096, endpoint=False)
    r = (np.stack([np.cos(th), np.sin(th)], -1)*rho) + (0 if about == 0 else c)
    den = np.sum(r*r, -1)*np.sum((c-r)**2, -1)
    M = np.einsum('am,aj->mj', r/den[:, None], (c-r))/len(th)
    return M
for about, nm in [(0, "r -> 0  (y'->z  /  z->y')"), (1, "r -> c  (x'->y' /  x->z)")]:
    for rho in (1e-1, 1e-2, 1e-3):
        M = ang(rho, about)
        print(f"   {nm:26s} rho={rho:7.0e}  rho^2*<integrand> max = {np.abs(M).max()*rho**2:.3e}")
print("   -> vanishes like rho^2 * (1/rho) * rho = rho^2 : integrand ~ 1/rho, and")
print("      int rho drho / rho is finite.  NO short-distance (UV) divergence.")

print()
print("2. LARGE DISTANCE.  angular average of the integrand at large rho")
th = np.linspace(0, 2*np.pi, 8192, endpoint=False)
for rho in (1e1, 1e2, 1e3, 1e4):
    r = np.stack([np.cos(th), np.sin(th)], -1)*rho
    den = np.sum(r*r, -1)*np.sum((c-r)**2, -1)
    M = np.einsum('am,aj->mj', r/den[:, None], (c-r))/len(th)
    print(f"   rho={rho:8.0e}  rho^2*<integrand> = [[{M[0,0]*rho**2:+.6f} {M[0,1]*rho**2:+.6f}]"
          f" [{M[1,0]*rho**2:+.6f} {M[1,1]*rho**2:+.6f}]]")
print("   -> tends to -delta^{mj}/2 : integrand ~ -delta^{mj}/(2 rho^2), and")
print("      int d^2r (-delta/(2 r^2)) = -pi delta^{mj} int drho/rho  :  LOG DIVERGENT.")
print("      the divergence is pure trace (traceless part angular-averages to zero).")

print()
print("3. the same statement from the exact algebraic identity")
print("   r.(c-r)/(r^2 (c-r)^2) = -1/2 [ 1/r^2 + 1/(c-r)^2 - c^2/(r^2 (c-r)^2) ]")
rng = np.random.default_rng(0); r = rng.normal(size=(5, 2))
lhs = np.sum(r*(c-r), -1)/(np.sum(r*r, -1)*np.sum((c-r)**2, -1))
rhs = -0.5*(1/np.sum(r*r, -1) + 1/np.sum((c-r)**2, -1)
            - (c@c)/(np.sum(r*r, -1)*np.sum((c-r)**2, -1)))
print(f"   max |lhs-rhs| = {np.abs(lhs-rhs).max():.2e}")
print("   the last term converges at large r; the first two give -1/r^2 : the log.")

print()
print("4. CUT-OFF FORM.  set Q=0 and cut the radius at R:")
print("   delta_mj int_{|r|<R} d^2r (r^m/r^2)((c-r)^j/(c-r)^2)  =  -2 pi log(R) + const")
def trunc(R):
    def f(rho):
        r = np.stack([np.cos(th), np.sin(th)], -1)*rho
        d = np.sum(r*r, -1)*np.sum((c-r)**2, -1)
        return rho*np.mean(np.sum(r*(c-r), -1)/d)*2*np.pi
    v = 0.0
    pts = [0.0, cn*0.999, cn*1.001, 2*cn] + list(np.geomspace(4*cn, R, 60))
    for a, b in zip(pts[:-1], pts[1:]):
        if b > R: b = R
        if b <= a: continue
        v += quad(f, a, b, limit=300)[0]
    return v
prev = None
for R in (1e2, 1e3, 1e4, 1e5):
    v = trunc(R)
    s = "" if prev is None else f"   d/dlogR = {(v-prev)/np.log(10):+.6f}  (-2pi = {-2*np.pi:+.6f})"
    print(f"   R={R:8.0e}   I(R) = {v:+.6f}{s}")
    prev = v

print()
print("5. WHAT REGULATES IT.  with Q != 0 the phase cuts the log at r ~ 2 e^-gamma/|Q|:")
print("   trace B(Q,c) = -2 pi log( 2 e^-gamma / (|Q||c|) )")
for sc in (1e-1, 1e-2, 1e-3, 1e-4):
    Q = sc*np.array([0.9, 0.6]); Qn = np.hypot(*Q)
    t = np.trace(Bten(Q, c)).real
    print(f"   |Q||c| = {Qn*cn:9.2e}   trace = {t:+.6f}   "
          f"-2pi log(2e^-g/(|Q||c|)) = {-2*np.pi*np.log(2*np.exp(-g)/(Qn*cn)):+.6f}")
print("   -> the transverse integral DIVERGES logarithmically as |Q| -> 0.")

print()
print("6. WHERE |Q| -> 0 HAPPENS")
print("   y'-integration: Q = (1+xi) k  ->  |Q| >= k  never zero  ->  finite")
print("   z -integration: Q = -xi k     ->  |Q| = xi k -> 0 at the soft end xi -> 0")
print("   so   log(1/(xi k |c|))^  multiplies the 1/xi of the T2 term:")
print("        int_lambda^Xi dxi/xi * 2 log(1/(xi k|c|))  =  log^2 ... : DOUBLE LOG")
for lam in (1e-2, 1e-3, 1e-4):
    I = quad(lambda u: 2*np.log(1/u)/u, lam, 1.0)[0]
    print(f"   lambda={lam:7.0e}:  int_lambda^1 dxi/xi 2log(1/xi) = {I:12.4f}   log^2(1/lambda) = {np.log(1/lam)**2:12.4f}")
