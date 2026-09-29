"""Row V part 2a: check the user's plus-prescription rewrite, term by term.

Everything is compared at FIXED transverse positions (x, y, z, x', w'), with the common factors
(1/(2pi)^3) (g^4/8pi^4)(1/k+)(1/2pi) and the colour factor stripped (they are p+-independent).
k+ = 1, so p+ = xi.
"""
import numpy as np
from scipy import integrate
from scipy.linalg import expm

rng = np.random.default_rng(7)
sq = lambda v: float(v @ v)
bar = "=" * 78

def blocks(x, y, z, xp, wp, xi):
    """the six blocks of the user's 2a bracket, contracted with (x'-w')^i (z-x)^m/[(x'-w')^2 (z-x)^2]; d_perp = 2."""
    I = np.eye(2); xb = 1 - xi
    B = [I*2/(2*sq(x-z))*(sq(y-x)-sq(y-z))*xi*xb,
         -(xi/xb)*(np.outer(y-x, x-z)/sq(x-z) - np.outer(y-x, y-z)/(2*sq(y-z))),
         xb*I*((y-z)@(x-z)/sq(x-z) + (y-x)@(y-z)/(2*sq(y-x)) - (sq(y-x)-sq(y-z))/(2*sq(x-z))),
         -(xb/xi)*(np.outer(y-z, x-z)/sq(x-z) + np.outer(y-z, y-x)/(2*sq(y-x))),
         xi*I*((y-x)@(x-z)/sq(x-z) - (y-x)@(y-z)/(2*sq(y-z)) - (sq(y-x)-sq(y-z))/(2*sq(x-z))),
         -(np.outer(x-z, y-x)/sq(x-z) - np.outer(y-z, y-x)/(2*sq(y-z))
           + np.outer(x-z, y-z)/sq(x-z) + np.outer(y-x, y-z)/(2*sq(y-x)))]
    Pv = (xp-wp)/sq(xp-wp); K = (z-x)/sq(z-x)
    return np.array([Pv@b@K for b in B])

def B2c(x, y, z, xp, wp):   # bracket of T2 without -(xi/xibar)
    Pv = (xp-wp)/sq(xp-wp); K = (z-x)/sq(z-x)
    return Pv@(np.outer(y-x, x-z)/sq(x-z) - np.outer(y-x, y-z)/(2*sq(y-z)))@K
def B4c(x, y, z, xp, wp):   # bracket of T4 without -(xibar/xi)
    Pv = (xp-wp)/sq(xp-wp); K = (z-x)/sq(z-x)
    return Pv@(np.outer(y-z, x-z)/sq(x-z) + np.outer(y-z, y-x)/(2*sq(y-x)))@K

def phase(k, x, z, wp, xi):
    return np.exp(-1j*(k@(wp-z)) - 1j*xi*(k@(z-x)))
def den(x, y, z, xi):
    return xi*sq(y-x) + (1-xi)*sq(y-z)

def cq(f, a, b, pts=None):
    g = lambda t: f(t)
    re = integrate.quad(lambda t: g(t).real, a, b, points=pts, limit=1000, epsabs=1e-13, epsrel=1e-12)[0]
    im = integrate.quad(lambda t: g(t).imag, a, b, points=pts, limit=1000, epsabs=1e-13, epsrel=1e-12)[0]
    return re + 1j*im

print(bar); print("1. POINTWISE CHECK OF THE REWRITE (fixed x,y,z,x',w')"); print(bar)
print("   LHS  = Int_Lam^{1-Lam} dxi phase/den * sum(T1..T6)")
print("   user = your 4 terms, exactly as written (T2 kept in the first part AND moved to the q+ part)")
print("   fix  = the same with T2 removed from the first part")
for trial in range(3):
    x, y, z, xp, wp = rng.normal(size=(5, 2)); k = rng.normal(size=2)
    F = lambda xi: phase(k, x, z, wp, xi)/den(x, y, z, xi)*blocks(x, y, z, xp, wp, xi).sum()
    g = lambda xi: phase(k, x, z, wp, xi)/den(x, y, z, xi)*B4c(x, y, z, xp, wp)          # T4 residue fn
    h = lambda q: phase(k, x, z, wp, 1-q)/den(x, y, z, 1-q)*B2c(x, y, z, xp, wp)        # T2 residue fn (q = 1-xi)
    def part1(xi, keepT2):
        b = blocks(x, y, z, xp, wp, xi)
        reg = b[0] + b[2] + b[4] + b[5] + (b[1] if keepT2 else 0) + B4c(x, y, z, xp, wp)   # +[B4] from -(1/xi - 1)
        return phase(k, x, z, wp, xi)/den(x, y, z, xi)*reg - (g(xi) - g(0))/xi
    part2 = -cq(lambda q: (h(q) - h(0))/q - h(q), 0, 1)       # Int_{k+}^0 dq {[1/q]_+ B2 - B2}
    for Lam in (1e-3, 1e-5):
        L = np.log((1 - Lam)/Lam)
        lhs = cq(F, Lam, 1 - Lam, pts=[10*Lam, 1 - 10*Lam])
        fix = cq(lambda t: part1(t, False), 0, 1) + part2 - L*g(0) - L*h(0)
        user = cq(lambda t: part1(t, True), 0, 1 - Lam) + part2 - L*g(0) - L*h(0)   # T2 kept: needs a cutoff at 1
        print(f"   config {trial}, Lambda={Lam:.0e}:  LHS {lhs:.6f}   fix {fix:.6f} (diff {abs(lhs-fix):.1e})"
              f"   user {user:.6f} (diff {abs(lhs-user):.1e})")
print("   -> the rewrite is right once T2 is removed from the first part; as written T2 is counted twice")
print("      and Int_0^{k+} of the leftover -(p+/(k+-p+))[...] diverges at p+ -> k+ (the 'user' column grows like log(1/Lambda)).")

# ----------------------------------------------------------------------------------------------
print(); print(bar); print("2. UV (z -> x): each new term, times the colour factor, is finite"); print(bar)
Nc = 3
l = np.zeros((8, 3, 3), dtype=complex)
l[0][0,1]=l[0][1,0]=1; l[1][0,1]=-1j; l[1][1,0]=1j; l[2][0,0]=1; l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1; l[4][0,2]=-1j; l[4][2,0]=1j; l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j; l[7]=np.diag([1,1,-2])/np.sqrt(3)
t = l/2
f = np.real(-2j*np.einsum('aij,bjk,cki->abc', t, t, t) + 2j*np.einsum('bij,ajk,cki->abc', t, t, t))
def adj(V): return np.real(2*np.einsum('aij,jk,bkl,li->ab', t, V, t, V.conj().T))
def randH():
    H = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)); H = (H + H.conj().T)/2
    return H - np.trace(H)*np.eye(3)/3
H1, H2 = randH(), randH()
def U(v, s_=0.5, R=1.6):
    gg = np.exp(-(v@v)/(2*R*R)); return adj(expm(1j*s_*gg*(v[0]*H1 + v[1]*H2)))
Wt = rng.normal(size=(8, 8))
def col(xv, zv):   # one fixed projection of -f f U(x) U(z) + Nc U(x)
    Ux, Uz = U(xv), U(zv)
    return float(np.sum(Wt*(-np.einsum('adc,fbe,db,ce->af', f, f, Ux, Uz) + Nc*Ux)))
x = np.array([0.3, -0.2]); y = np.array([1.1, 0.4]); xp = np.array([-0.5, 0.9]); wp = np.array([0.7, 1.3])
print("   rho^2 < term * colour >_angles  as z -> x  (a constant would mean a log UV divergence):")
for name, fn in (("delta(p+) term  B4/(y-z)^2", lambda zv: B4c(x, y, zv, xp, wp)/sq(y-zv)),
                 ("delta(q+) term  B2/(y-x)^2", lambda zv: B2c(x, y, zv, xp, wp)/sq(y-x))):
    out = []
    for rho in (1e-2, 1e-3, 1e-4):
        phs = np.linspace(0, 2*np.pi, 256, endpoint=False)
        vals = [fn(x + rho*np.array([np.cos(p), np.sin(p)]))*col(x, x + rho*np.array([np.cos(p), np.sin(p)])) for p in phs]
        out.append(rho**2*np.mean(vals))
    print(f"     {name}:  " + "  ".join(f"{v:.2e}" for v in out))
print("   -> falls with rho: no UV divergence (the colour factor vanishes at z = x).")

# ----------------------------------------------------------------------------------------------
print(); print(bar); print("3. THE delta(q+) (p+ -> k+) LOG TERM: ITS z INTEGRAL DIVERGES AT LARGE |z|"); print(bar)
print("   its phase e^{-ik.(w'-x)} has no z, and the colour factor -> C_inf(x) != 0 as U(z) -> 1.")
W = (xp-wp)@(y-x)/(sq(xp-wp)*sq(y-x))
Cinf = float(np.sum(Wt*(-np.einsum('adc,fbe,db,ce->af', f, f, U(x), np.eye(8)) + Nc*U(x))))
for R in (1e1, 1e2, 1e3):
    phs = np.linspace(0, 2*np.pi, 512, endpoint=False)
    vals = [B2c(x, y, x + R*np.array([np.cos(p), np.sin(p)]), xp, wp)/sq(y-x) for p in phs]
    print(f"     |z-x| = {R:.0e}:  |z-x|^2 < B2/(y-x)^2 > = {R**2*np.mean(vals):+.6f}    predicted -W/2 = {-W/2:+.6f}")
print(f"   colour at large z: C_inf = {Cinf:+.4f} (not 0).  So Int d^2z ~ -pi W C_inf log(R): log divergent.")

# ----------------------------------------------------------------------------------------------
print(); print(bar); print("4. THE delta(p+) (p+ -> 0) LOG TERM: z INTEGRAL FINITE, x INTEGRAL LOG DIVERGENT"); print(bar)
z = np.array([0.2, 0.5])
print("   large |z| : phase e^{ik.z} oscillates and the term falls like 1/|z|^2 -> finite.")
print("   large |x| : phase e^{-ik.(w'-z)} has no x; colour -> -f f 1 U(z) + Nc 1 != 0 :")
for R in (1e1, 1e2, 1e3):
    phs = np.linspace(0, 2*np.pi, 512, endpoint=False)
    vals = [B4c(R*np.array([np.cos(p), np.sin(p)]), y, z, xp, wp)/sq(y-z) for p in phs]
    pred = -0.5*((xp-wp)@(y-z))/(sq(xp-wp)*sq(y-z))
    print(f"     |x| = {R:.0e}:  |x|^2 < B4/(y-z)^2 > = {R**2*np.mean(vals):+.6f}    predicted -P.(y-z)/(2 P^2 (y-z)^2) = {pred:+.6f}")
Cx = float(np.sum(Wt*(-np.einsum('adc,fbe,db,ce->af', f, f, np.eye(8), U(z)) + Nc*np.eye(8))))
print(f"   colour at large x: {Cx:+.4f} (not 0).  So Int d^2x of this term is log divergent at large |x|,")
print("   unless it is combined with the other diagrams (as in the derivation of the JIMWLK/BK kernel).")
