"""Checks for the note on Group I, Row V (RowV_UV.pdf).

 1. the contraction C1-vertex x B2-bracket = the bracket of Group_V.tex (d_perp = 2 and 3)
 2. the prefactor g^4/(8 pi^4) p+(k+-p+)/(k+^3 D)
 3. colour identities (general adjoint U, SU(3))
 4. small-r behaviour of the six brackets, angular average in d = 2 and d = 3
 5. full kernel V at small r: r^2 <V> -> Phihat(xi) (y-z)/(y-z)^2
 6. part 2a: colour structure vanishes linearly at x = z
 7. xi integrals (sympy): int Phihat, plus-prescribed remainder, int c(xi) log xi
 8. dimensional-regularization integrals (sympy, MS-bar)
 9. 2D identity between the phase and the cut regulator
10. part 2b: after the plus prescription and the counterterm, the x-integrand is UV finite
-> rowV_checks.json"""
import json
import os
import sys
import numpy as np
import sympy as sp
from scipy import integrate, special, linalg

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'three_rho'))
import numerics as nm  # noqa: E402

out = {}
rng = np.random.default_rng(2026)


def sq(v):
    return float(v @ v)


# ------------------------------------------------------------------ brackets of Group_V.tex
def brackets(x, y, z, dperp):
    """B1..B6 of the notes (index order i, m), vectors of any dimension"""
    n = len(x)
    I = np.eye(n)
    yx, yz, xz = y - x, y - z, x - z
    B1 = I * dperp / (2 * sq(xz)) * (sq(yx) - sq(yz))
    B2 = np.outer(yx, xz) / sq(xz) - np.outer(yx, yz) / (2 * sq(yz))
    B3 = I * (yz @ xz / sq(xz) + yx @ yz / (2 * sq(yx)) - (sq(yx) - sq(yz)) / (2 * sq(xz)))
    B4 = np.outer(yz, xz) / sq(xz) + np.outer(yz, yx) / (2 * sq(yx))
    B5 = I * (yx @ xz / sq(xz) - yx @ yz / (2 * sq(yz)) - (sq(yx) - sq(yz)) / (2 * sq(xz)))
    B6 = (np.outer(xz, yx) / sq(xz) - np.outer(yz, yx) / (2 * sq(yz)) + np.outer(xz, yz) / sq(xz)
          + np.outer(yx, yz) / (2 * sq(yx)))
    return [B1, B2, B3, B4, B5, B6]


def coeffs(xi):
    return [xi * (1 - xi), -xi / (1 - xi), 1 - xi, -(1 - xi) / xi, xi, -1.0]


def line56(x, y, z, kp, pp, dperp):
    """the bracket of Eq. (line 56) of Group_V.tex (before the factor p+(k+-p+)/k+^2)"""
    B = brackets(x, y, z, dperp)
    q = kp - pp
    return (B[0] - kp**2 / q**2 * B[1] + kp / pp * B[2] - kp**2 / pp**2 * B[3] + kp / q * B[4] - kp**2 / (pp * q) * B[5])


# 1. contraction
res = []
for n in (2, 3):
    for trial in range(3):
        x, y, z = rng.normal(size=(3, n))
        kp = 1.7
        pp = rng.uniform(0.1, 1.6)
        q = kp - pp
        I = np.eye(n)
        yx, yz, xz = y - x, y - z, x - z
        # B2 bracket with J <-> gluon p+ at x, K <-> gluon k+-p+ at z
        Bb = (I * (sq(yx) - sq(yz)) / (2 * sq(xz)) + kp / q * (np.outer(yx, xz) / sq(xz) - np.outer(yx, yz) / (2 * sq(yz)))
              + kp / pp * (np.outer(xz, yz) / sq(xz) + np.outer(yx, yz) / (2 * sq(yx))))
        # vertex C_{i J K; m} = d_im d_JK - (k+/p+) d_Jm d_iK - (k+/(k+-p+)) d_Km d_iJ
        T = np.einsum('im,jk,jk->im', I, I, Bb) - kp / pp * np.einsum('jm,ik,jk->im', I, I, Bb) - kp / q * np.einsum('km,ij,jk->im', I, I, Bb)
        res.append(float(np.abs(T - line56(x, y, z, kp, pp, n)).max() / np.abs(T).max()))
out['contraction_rel_diff'] = max(res)

# 2. prefactor
g, kp_, pp_, D = sp.symbols('g k_p p_p D', positive=True)
A = g / (sp.sqrt(2) * sp.pi * sp.sqrt(kp_))                       # |A^(1)| without the kernel
C = g / (2 * sp.pi * sp.sqrt(2 * kp_)) * sp.sqrt(pp_ * (kp_ - pp_)) / kp_
B = g**2 * sp.sqrt(pp_ * (kp_ - pp_)) / (4 * sp.pi**2 * D * kp_)     # B2 with k+ -> k+-p+, p+ -> p+, (p+ + k+) -> k+
pref = sp.simplify(2 * A * C * B / (g**4 / (8 * sp.pi**4) * pp_ * (kp_ - pp_) / (kp_**3 * D)))
out['prefactor_ratio'] = str(pref)

# 3. colour
ts, f = nm.structure_constants(3)
Nc = 3
U = nm.random_adjoint(3, ts, rng, False)
lhs = np.einsum('adc,fbe,db,ce->af', f, f, U, U)
out['ffUU_same_point'] = float(np.abs(lhs - Nc * U).max())
out['ff_contraction'] = float(np.abs(np.einsum('abc,cbd->ad', f, f) + Nc * np.eye(8)).max())
out['f_antisym'] = float(np.abs(np.einsum('ebf->fbe', f) + f).max())

# 4./5. small-r behaviour


def sphere(n, m):
    """quadrature nodes on the unit circle/sphere (equal weights; exact for low-order polynomials)"""
    if n == 2:
        th = 2 * np.pi * (np.arange(m) + 0.5) / m
        return np.stack([np.cos(th), np.sin(th)], 1), np.full(m, 1.0 / m)
    ct, wt = np.polynomial.legendre.leggauss(m)
    ph = 2 * np.pi * (np.arange(2 * m) + 0.5) / (2 * m)
    pts, wts = [], []
    for c, w in zip(ct, wt):
        s_ = np.sqrt(1 - c * c)
        for p in ph:
            pts.append([s_ * np.cos(p), s_ * np.sin(p), c])
            wts.append(w / 2 / (2 * m))
    return np.array(pts), np.array(wts)


def ang_avg_fast(n, r, Y, xi, ndir):
    dirs, wts = sphere(n, ndir)
    z = np.zeros(n)
    terms = np.zeros((6, n))
    fullV = np.zeros(n)
    cs = coeffs(xi)
    for u, wq in zip(dirs, wts):
        xv = u * r
        Bs = brackets(xv, Y, z, n)
        w = -xv / (xv @ xv)
        Dn = xi * sq(Y - xv) + (1 - xi) * sq(Y)
        for k in range(6):
            t = Bs[k] @ w
            terms[k] += wq * t
            fullV += wq * cs[k] * t / Dn
    return terms * r * r, fullV * r * r


expect = {2: [1, -1, -1, -1, -1, -1], 3: [1, -1, -2 / 3, -1, -2 / 3, -2 / 3]}
small = {}
for n in (2, 3):
    Y = rng.normal(size=n)
    xi = 0.37
    rows = []
    for r in (1e-2, 1e-3, 1e-4):
        terms, fullV = ang_avg_fast(n, r, Y, xi, 64)
        coef = [float(t @ Y / (Y @ Y)) for t in terms]
        rows.append(dict(r=r, coef=coef))
    small[n] = rows
    Phat = xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi
    out['small_r_d%d' % n] = dict(rows=rows, expected=expect[n])
# full kernel V (d = 2), several xi
fullrows = []
Y = np.array([0.8, -0.5])
for xi in (0.1, 0.37, 0.8):
    Phat = xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi
    for r in (1e-3, 1e-4):
        _, fv = ang_avg_fast(2, r, Y, xi, 64)
        fullrows.append(dict(xi=xi, r=r, ratio=[float(v) for v in fv * (Y @ Y) / Y], Phat=Phat))
out['full_kernel_small_r'] = fullrows

# 6. part 2a: colour structure vanishes linearly at x = z
gens = [f[a] for a in range(8)]       # real antisymmetric generators f^{abc}: U = exp(theta^a f^a) is an adjoint rotation
Agen = [sum(rng.normal() * gens[a] for a in range(8)) for _ in range(2)]


def Ufield(p):
    return linalg.expm(p[0] * Agen[0] + p[1] * Agen[1])


zz = rng.normal(size=2)
o2a = []
for r in (1e-1, 1e-2, 1e-3):
    xx = zz + r * np.array([0.6, 0.8])
    O = -np.einsum('adc,fbe,db,ce->af', f, f, Ufield(xx), Ufield(zz)) + Nc * Ufield(xx)
    o2a.append(dict(r=r, norm=float(np.abs(O).max())))
out['O2a_small_r'] = o2a

# 7. xi integrals
xi_, L = sp.symbols('xi Lambda', positive=True)
Phat_s = xi_ * (1 - xi_) + xi_ / (1 - xi_) + (1 - xi_) / xi_
intPhat = sp.series(sp.simplify(sp.integrate(Phat_s, (xi_, L, 1 - L))), L, 0, 1).removeO()
out['int_Phat'] = str(sp.simplify(intPhat))
c_s = xi_ * (1 - xi_) - 2
out['int_plus_remainder'] = str(sp.integrate(sp.simplify(Phat_s - 1 / xi_ - 1 / (1 - xi_)), (xi_, 0, 1)))
out['int_c_logxi'] = str(sp.integrate(c_s * sp.log(xi_), (xi_, 0, 1)))
out['int_c_log_xi1mxi'] = str(sp.simplify(sp.integrate(c_s * sp.log(xi_ * (1 - xi_)), (xi_, 0, 1))))

# 8. dim reg (d = 2 - eps, MS-bar factor (mu^2 e^gamma/(4 pi))^(eps/2))
eps, R, q, mu = sp.symbols('epsilon R q mu', positive=True)
d = 2 - eps
cd2 = sp.gamma(d / 2)**2 * sp.pi**(2 - d)
Omega = 2 * sp.pi**(d / 2) / sp.gamma(d / 2)
msb = (mu**2 * sp.exp(sp.EulerGamma) / (4 * sp.pi))**(eps / 2)
I_cut = cd2 * Omega * R**eps / eps * msb                              # int_{r<R} d^d r c_d^2 (r^2)^(1-d)
a = 1 - eps
I_phase = cd2 * sp.pi**(d / 2) * 2**(d - 2 * a) * sp.gamma(d / 2 - a) / sp.gamma(a) * (q**2)**(a - d / 2) * msb
ser_cut = sp.series(I_cut, eps, 0, 1).removeO()
ser_phase = sp.series(I_phase, eps, 0, 1).removeO()
target_cut = sp.pi * (2 / eps + sp.log(R**2 * mu**2 * sp.exp(2 * sp.EulerGamma) / 4))
target_phase = sp.pi * (2 / eps + sp.log(mu**2 / q**2))
out['dimreg_cut_minus_target'] = str(sp.simplify(sp.expand_log(ser_cut - target_cut, force=True)))
out['dimreg_phase_minus_target'] = str(sp.simplify(sp.expand_log(ser_phase - target_phase, force=True)))

# 9. 2D identity: int d^2 r [e^{iq.r} - theta(R-r)]/r^2 = -pi log(q^2 R^2 e^{2 gamma}/4)
ident = []
for qv, Rv in ((1.0, 0.7), (3.0, 2.0), (0.2, 5.0)):
    f1 = integrate.quad(lambda r: (special.j0(qv * r) - 1.0) / r, 0, Rv, limit=400)[0]
    f2 = integrate.quad(lambda r: special.j0(qv * r) / r, Rv, np.inf, limit=4000)[0] if False else \
        -special.j0(qv * Rv) * 0  # placeholder replaced below
    # int_R^inf J0(qr)/r dr = int_{qR}^inf J0(t)/t dt
    tail = integrate.quad(lambda t: special.j0(t) / t, qv * Rv, qv * Rv + 2000, limit=20000)[0]
    lhs = 2 * np.pi * (f1 + tail)
    rhs = -np.pi * np.log(qv**2 * Rv**2 * np.exp(2 * np.euler_gamma) / 4)
    ident.append(dict(q=qv, R=Rv, lhs=lhs, rhs=float(rhs)))
out['phase_vs_cut_identity'] = ident


# 10. part 2b: plus prescription (xi -> 0, 1) at fixed x, then the counterterm c(xi) e^{i xi k.r} K(y-z)/r^2
def Vvec(xi, xv, Y):
    z = np.zeros(2)
    Bs = brackets(xv, Y, z, 2)
    w = -xv / (xv @ xv)
    Dn = xi * sq(Y - xv) + (1 - xi) * sq(Y)
    cs = coeffs(xi)
    return sum(c * (B @ w) for c, B in zip(cs, Bs)) / Dn


def plus_integral(xv, Y, kv):
    z = np.zeros(2)
    Bs = brackets(xv, Y, z, 2)
    w = -xv / (xv @ xv)
    V0 = -(Bs[3] @ w) / sq(Y)              # lim xi->0 of xi V
    V1 = -(Bs[1] @ w) / sq(Y - xv)         # lim xi->1 of (1-xi) V
    ph = lambda xi: np.exp(1j * xi * (kv @ xv))
    out_ = np.zeros(2, complex)
    for comp in range(2):
        f = lambda xi: Vvec(xi, xv, Y)[comp] * ph(xi) - V0[comp] / xi - V1[comp] * ph(1.0) / (1 - xi)
        re = integrate.quad(lambda t: f(t).real, 0, 1, limit=400, epsabs=1e-11, epsrel=1e-11)[0]
        im = integrate.quad(lambda t: f(t).imag, 0, 1, limit=400, epsabs=1e-11, epsrel=1e-11)[0]
        out_[comp] = re + 1j * im
    return out_


def ct_integral(xv, Y, kv):
    cfun = lambda xi: (xi * (1 - xi) - 2) * np.exp(1j * xi * (kv @ xv))
    re = integrate.quad(lambda t: cfun(t).real, 0, 1)[0]
    im = integrate.quad(lambda t: cfun(t).imag, 0, 1)[0]
    return (re + 1j * im) * Y / sq(Y) / sq(xv)


Y = np.array([0.8, -0.5]); kv = np.array([1.3, 0.4])
uv = []
for r in (1e-1, 1e-2, 1e-3):
    dirs, wts = sphere(2, 32)
    avgP = np.zeros(2, complex); avgD = np.zeros(2, complex); mxD = 0.0
    for u, wq in zip(dirs, wts):
        xv = u * r
        P = plus_integral(xv, Y, kv)
        Dd = P - ct_integral(xv, Y, kv)
        avgP += wq * P; avgD += wq * Dd
        mxD = max(mxD, float(np.abs(Dd).max()))
    uv.append(dict(r=r, r2_avg_plus=[float(v) for v in (r * r * avgP * sq(Y) / Y).real],
                   r2_avg_after_ct=float(np.abs(r * r * avgD).max()), r_max_after_ct=r * mxD))
out['plus_remainder_uv'] = uv

# 11. plus prescription at fixed x: direct integral over (lam, 1-lam) = plus integral + endpoint logs + O(lam)
pp_rows = []
xv = np.array([0.35, 0.2]); Y = np.array([0.8, -0.5]); kv = np.array([1.3, 0.4])
z = np.zeros(2)
Bs = brackets(xv, Y, z, 2); w = -xv / (xv @ xv)
V0 = -(Bs[3] @ w) / sq(Y); V1 = -(Bs[1] @ w) / sq(Y - xv)
P = plus_integral(xv, Y, kv)
for lam in (1e-2, 1e-3, 1e-4):
    direct = np.zeros(2, complex)
    for comp in range(2):
        fr = lambda t: (Vvec(t, xv, Y)[comp] * np.exp(1j * t * (kv @ xv)))
        direct[comp] = (integrate.quad(lambda t: fr(t).real, lam, 1 - lam, limit=400, points=[lam * 10, 1 - lam * 10])[0]
                        + 1j * integrate.quad(lambda t: fr(t).imag, lam, 1 - lam, limit=400, points=[lam * 10, 1 - lam * 10])[0])
    L_ = np.log((1 - lam) / lam)
    pred = P + V0 * L_ + V1 * np.exp(1j * (kv @ xv)) * L_
    pp_rows.append(dict(lam=lam, diff=float(np.abs(direct - pred).max()), size=float(np.abs(direct).max())))
out['plus_prescription_decomposition'] = pp_rows

json.dump(out, open(os.path.join(HERE, 'rowV_checks.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k in ('plus_prescription_decomposition',)}, indent=0))
