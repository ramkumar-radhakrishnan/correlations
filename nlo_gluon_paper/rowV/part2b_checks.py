"""Checks for Part2b_MSbar.pdf (Group I, Row V, part 2b in dimensional regularization, MS-bar).

 1. short distance x -> z: r^2 <V> -> [xi(1-xi) + xi/(1-xi) + (1-xi)/xi] (y-z)/(y-z)^2, in d = 2 and d = 3 (d_perp = d)
 2. the counterterm kernel S = K(z,z,x) - K(y,z,x): its angular average is theta(|y-z|-|x-z|)/(x-z)^2 (2 dims)
 3. the rest (original minus counterterm, smooth test colour factor, phase kept): r^2 <rest> -> 0
 4. the x integral of the rest converges at small |x-z|; the original grows like 2 pi log(1/r0)
 5. the rest falls like 1/|x-z|^3 at large |x-z|
 6. dimensional regularization (sympy): theta form and K form of the counterterm, MS-bar
 7. exact p+ integral (sympy)
 8. prefactors (sympy): UV term = (alpha_s N_c / 2 pi) x LO; beta_0 part cancels the one of G1-II(1)
-> part2b_checks.json"""
import json
import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}


def sq(v):
    return float(v @ v)


# ------------------------------------------------------------------ the bracket of part 2b (as in rowV_checks.py)
def brackets(x, y, z, dperp):
    """B1..B6 of the six-term bracket (index order i, m), vectors of any dimension"""
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
    """p+(k+-p+)/k+^2, -p+/(k+-p+), (k+-p+)/k+, -(k+-p+)/p+, p+/k+, -1"""
    return [xi * (1 - xi), -xi / (1 - xi), 1 - xi, -(1 - xi) / xi, xi, -1.0]


def phihat(xi):
    return xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi


def V(xi, x, y, z, dperp):
    """(z-x)^m Bracket^{im} / [xi (y-x)^2 + (1-xi)(y-z)^2]  (= k+ x the integrand without prefactor, colour, phase)"""
    Bs = brackets(x, y, z, dperp)
    w = (z - x) / sq(z - x)
    return sum(c * (B @ w) for c, B in zip(coeffs(xi), Bs)) / (xi * sq(y - x) + (1 - xi) * sq(y - z))


def S(x, y, z):
    """K(z,z,x) - K(y,z,x) = (y-x).(y-z) / [(y-x)^2 (x-z)^2]"""
    return (y - x) @ (y - z) / (sq(y - x) * sq(x - z))


def Vvec(xi, X, Yv, dperp=2):
    """V for many points X (N,2), z = 0, y = Yv; same as V()"""
    r2 = np.einsum('ij,ij->i', X, X)
    W = -X / r2[:, None]
    YX = Yv[None, :] - X
    yx2 = np.einsum('ij,ij->i', YX, YX)
    Y2 = Yv @ Yv
    dot = lambda A, B: np.einsum('ij,ij->i', A, B)
    Yb = np.broadcast_to(Yv, X.shape)
    xw, Yw, yxw = dot(X, W), W @ Yv, dot(YX, W)
    B1 = (dperp / (2 * r2) * (yx2 - Y2))[:, None] * W
    B2 = YX * (xw / r2)[:, None] - YX * (Yw / (2 * Y2))[:, None]
    B3 = (X @ Yv / r2 + YX @ Yv / (2 * yx2) - (yx2 - Y2) / (2 * r2))[:, None] * W
    B4 = Yb * (xw / r2)[:, None] + Yb * (yxw / (2 * yx2))[:, None]
    B5 = (dot(YX, X) / r2 - YX @ Yv / (2 * Y2) - (yx2 - Y2) / (2 * r2))[:, None] * W
    B6 = X * (yxw / r2)[:, None] - Yb * (yxw / (2 * Y2))[:, None] + X * (Yw / r2)[:, None] + YX * (Yw / (2 * yx2))[:, None]
    c = coeffs(xi)
    tot = c[0] * B1 + c[1] * B2 + c[2] * B3 + c[3] * B4 + c[4] * B5 + c[5] * B6
    return tot / (xi * yx2 + (1 - xi) * Y2)[:, None]


def Svec(X, Yv):
    YX = Yv[None, :] - X
    return (YX @ Yv) / (np.einsum('ij,ij->i', YX, YX) * np.einsum('ij,ij->i', X, X))


def sphere(n, m):
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


# 1. short distance, d = 2 and d = 3
rng = np.random.default_rng(7)
small = []
for n in (2, 3):
    Y = rng.normal(size=n)
    z = np.zeros(n)
    dirs, wts = sphere(n, 64)
    for xi in (0.1, 0.37, 0.8):
        r = 1e-4
        avg = sum(wq * V(xi, u * r, Y, z, n) for u, wq in zip(dirs, wts))
        ratio = r * r * avg * sq(Y) / Y
        small.append(dict(d=n, xi=xi, ratio=float(np.max(np.abs(ratio - phihat(xi)))) / phihat(xi), phihat=phihat(xi)))
out['short_distance'] = small

# 2. angular average of S
Y = np.array([0.8, -0.5]); z = np.zeros(2); R = np.sqrt(sq(Y))
dirs, wts = sphere(2, 512)
sav = []
for r in (0.01, 0.3, 0.9 * R, 1.1 * R, 3.0, 30.0):
    avg = sum(wq * S(u * r, Y, z) for u, wq in zip(dirs, wts))
    sav.append(dict(r_over_R=r / R, r2_avg=float(r * r * avg), expected=1.0 if r < R else 0.0))
out['S_angular_average'] = sav

# 3./4./5. the rest with a smooth test colour factor F(x) = G(y) - G(x) and the phase e^{i xi k.(x-z)} (z = 0)
kv = np.array([1.3, 0.4])


def Gv(X):
    X = np.atleast_2d(X)
    return np.cos(0.7 * X[:, 0] - 0.4 * X[:, 1]) + 0.5 * np.sin(0.3 * X[:, 0] + 0.9 * X[:, 1])


def Fv(X):
    return Gv(Y)[0] - Gv(X)


def orig(xi, X, phase=True):
    ph = np.exp(1j * xi * (X @ kv)) if phase else np.ones(len(X))
    return (ph * Fv(X))[:, None] * Vvec(xi, X, Y)


def ct(xi, X):
    return phihat(xi) * (Y / sq(Y))[None, :] * (Svec(X, Y) * Fv(z)[0])[:, None]


def rest(xi, X):
    return orig(xi, X) - ct(xi, X)


# vectorized kernel = loop kernel
Xr = rng.normal(size=(20, 2))
out['Vvec_vs_V'] = max(float(np.abs(Vvec(xi, Xr, Y)[j] - V(xi, Xr[j], Y, z, 2)).max() / np.abs(V(xi, Xr[j], Y, z, 2)).max())
                       for xi in (0.1, 0.6) for j in range(20))

dirs, wts = sphere(2, 128)
r2rows = []
for xi in (0.02, 0.37, 0.95):
    for r in (1e-2, 1e-3, 1e-4):
        ao = wts @ orig(xi, dirs * r)
        ar = wts @ rest(xi, dirs * r)
        r2rows.append(dict(xi=xi, r=r, r2_orig=float(np.abs(r * r * ao).max()), r2_rest=float(np.abs(r * r * ar).max())))
out['rest_small_r'] = r2rows

gl_t, gl_w = np.polynomial.legendre.leggauss(24)


def ring_integral(fun, xi, edges, ndir):
    """int d^2x fun over edges[0] < |x-z| < edges[-1]: Gauss-Legendre in r on each panel, trapezoid in angle"""
    dr, wr = sphere(2, ndir)
    tot = np.zeros(2, complex)
    for a_, b_ in zip(edges[:-1], edges[1:]):
        rr = 0.5 * (b_ - a_) * gl_t + 0.5 * (b_ + a_)
        P = (rr[:, None, None] * dr[None, :, :]).reshape(-1, 2)
        vals = fun(xi, P).reshape(len(rr), ndir, 2)
        ang = np.einsum('j,ijk->ik', wr, vals) * 2 * np.pi
        tot += np.einsum('i,i,ik->k', 0.5 * (b_ - a_) * gl_w, rr, ang)
    return tot


# x integral over r0 < |x-z| < 0.3 (inside |y-z|): rest converges, original grows like log(1/r0^2)
conv = []
for xi in (0.2, 0.7):
    for r0 in (1e-2, 1e-3, 1e-4, 1e-5):
        edges = np.geomspace(r0, 0.3, int(3 * np.log10(0.3 / r0)) + 3)
        Ir = ring_integral(rest, xi, edges, 128)
        Io = ring_integral(orig, xi, edges, 128)
        logpart = phihat(xi) * Fv(z)[0] * Y / sq(Y) * np.pi * np.log(0.09 / r0**2)
        conv.append(dict(xi=xi, r0=r0, rest=[float(Ir[0].real), float(Ir[0].imag)], orig=float(Io[0].real),
                         orig_minus_log=float((Io - logpart)[0].real)))
out['rest_x_integral'] = conv

# large |x-z|: counterterm ~ 1/r^3; original (= rest there) ~ 1/r^2; without the phase its angular average is not zero,
# with the phase e^{i xi k.(x-z)} the x integral converges
large = []
dl, wl = sphere(2, 256)
for xi in (0.2, 0.7):
    for r in (1e1, 1e2, 1e3):
        P = dl * r
        large.append(dict(xi=xi, r=r, r3_max_ct=float(np.abs(ct(xi, P)).max() * r**3),
                          r2_max_orig=float(np.abs(orig(xi, P)).max() * r**2),
                          r2_avg_nophase_F1=float(np.abs(wl @ Vvec(xi, P, Y)).max() * r**2)))
out['large_r'] = large
tail = []
for xi in (0.2, 0.7):
    for Rmax in (25.0, 50.0, 100.0, 200.0):
        edges = np.linspace(1.5, Rmax, int(np.ceil((Rmax - 1.5) / 0.5)) + 1)
        T = ring_integral(rest, xi, edges, 2048)
        tail.append(dict(xi=xi, Rmax=Rmax, tail=[float(T[0].real), float(T[0].imag), float(T[1].real), float(T[1].imag)]))
out['rest_tail'] = tail

# 6. dimensional regularization, d = 2 - eps, MS-bar factor (mu^2 e^gamma / (4 pi))^(eps/2) per transverse integral
eps, Rs, mu = sp.symbols('epsilon R mu', positive=True)
d = 2 - eps
msb = (mu**2 * sp.exp(sp.EulerGamma) / (4 * sp.pi))**(eps / 2)
cd2 = sp.gamma(d / 2)**2 * sp.pi**(2 - d)
Omega = 2 * sp.pi**(d / 2) / sp.gamma(d / 2)
target = sp.pi * (2 / eps + 2 * sp.EulerGamma + sp.log(Rs**2 * mu**2 / 4))
# theta form: int_{|x-z|<R} d^dx c_d^2 ((x-z)^2)^(1-d)
I_theta = msb * cd2 * Omega * Rs**eps / eps
# K form: int d^dx [K_d(z,z,x) - K_d(y,z,x)] = 0 - 4 pi^2 int d^dq/(2pi)^d e^{iq(y-z)}/q^2
I_K = -msb * 4 * sp.pi**2 * sp.gamma(d / 2 - 1) / (4 * sp.pi**(d / 2)) * (Rs**2)**(1 - d / 2)
for name, I in (('theta', I_theta), ('K', I_K)):
    ser = sp.series(I, eps, 0, 1).removeO()
    out['dimreg_%s_minus_target' % name] = str(sp.simplify(sp.expand_log(ser - target, force=True)))
# d = 3 value of the K-form formula: int d^3x c_3^2 (y-x).(z-x)/(|y-x|^3|z-x|^3) = pi/R  (from int grad(1/r).grad(1/r') = 4 pi/R)
dd = 3
val3 = 4 * sp.pi**2 * sp.gamma(sp.Rational(dd, 2) - 1) / (4 * sp.pi**sp.Rational(dd, 2)) * (Rs**2)**(1 - sp.Rational(dd, 2))
c3sq = (sp.gamma(sp.Rational(3, 2)) * sp.pi**(1 - sp.Rational(3, 2)))**2
out['K_formula_d3_minus_known'] = str(sp.simplify(val3 - c3sq * 4 * sp.pi / Rs))
# smooth form used before: R^2/(r^2 (r^2+R^2)): radial integral int_0^inf t^(1-d)/(1+t^2) dt = pi/(2 sin(pi eps/2))
I_smooth = msb * cd2 * Omega * Rs**eps * sp.pi / (2 * sp.sin(sp.pi * eps / 2))
out['dimreg_smooth_minus_target'] = str(sp.simplify(sp.expand_log(sp.series(I_smooth, eps, 0, 1).removeO() - target, force=True)))

# 7. exact p+ integral
xi_, lam = sp.symbols('xi lambda', positive=True)
Ph = xi_ * (1 - xi_) + xi_ / (1 - xi_) + (1 - xi_) / xi_
exact = sp.integrate(Ph, (xi_, lam, 1 - lam))
claim = 2 * sp.log((1 - lam) / lam) - sp.Rational(11, 6) + 4 * lam - lam**2 + sp.Rational(2, 3) * lam**3
out['xi_integral_minus_claim'] = max(abs(float(sp.N((exact - claim).subs(lam, lv), 30))) for lv in (sp.Rational(1, 3), sp.Rational(1, 10), sp.Rational(1, 1000)))
lv = 0.013
out['xi_integral_numeric'] = [float(sp.N(exact.subs(lam, lv))), float(sp.N(claim.subs(lam, lv)))]

# 8. prefactors
g, Nc, kp = sp.symbols('g N_c k_p', positive=True)
alpha = g**2 / (4 * sp.pi)
pref_UV = 1 / (2 * sp.pi)**3 * g**4 / (8 * sp.pi**4) / kp / (2 * sp.pi) * sp.pi * Nc   # (g^4/8pi^4)(1/k+)(dp+/2pi)(1/k+) x k+ x pi x N_c
out['pref_UV_over_g4Nc_16pi4'] = str(sp.simplify(pref_UV / (1 / (2 * sp.pi)**3 * g**4 * Nc / (16 * sp.pi**4) / kp)))
P_LO = 1 / (2 * sp.pi)**3 * g**2 / (2 * sp.pi**2 * kp)
out['pref_UV_over_LO'] = str(sp.simplify(pref_UV / P_LO / (alpha * Nc / (2 * sp.pi))))
beta0 = sp.Rational(11, 3) * Nc
out['beta0_part'] = str(sp.simplify(alpha * Nc / (2 * sp.pi) * sp.Rational(-11, 6) + alpha * beta0 / (4 * sp.pi)))
# G1-II(1): (g^4 N_c/32 pi^4)(1/k+) [-(bracket)] (11/3) x colour [U(w)-U(x)] = -colour of LO  ->  +(g^4 N_c/32 pi^4)(11/3) per LO-colour
G1II = g**4 * Nc / (32 * sp.pi**4) * (-1) * sp.Rational(11, 3) * (-1)
V2b = g**4 * Nc / (16 * sp.pi**4) * sp.Rational(-11, 6)
out['beta0_G1II_plus_2b'] = str(sp.simplify(G1II + V2b))

json.dump(out, open(os.path.join(HERE, 'part2b_checks.json'), 'w'), indent=1)
for k in ('dimreg_theta_minus_target', 'dimreg_K_minus_target', 'dimreg_smooth_minus_target', 'K_formula_d3_minus_known',
          'xi_integral_minus_claim', 'xi_integral_numeric', 'pref_UV_over_g4Nc_16pi4', 'pref_UV_over_LO', 'beta0_part', 'beta0_G1II_plus_2b'):
    print(k, out[k])
print('short distance max rel dev', max(r['ratio'] for r in small))
print('Vvec vs V', out['Vvec_vs_V'])
print('S angular', [(round(r['r_over_R'], 3), round(r['r2_avg'], 12)) for r in sav])
for r in r2rows:
    print('xi %.2f r %.0e  r2|<orig>| %.4f  r2|<rest>| %.2e' % (r['xi'], r['r'], r['r2_orig'], r['r2_rest']))
for r in conv:
    print('xi %.1f r0 %.0e  int rest %s   orig %.5f  orig - log %.6f' % (r['xi'], r['r0'], r['rest'], r['orig'], r['orig_minus_log']))
for r in large:
    print(r)
for r in tail:
    print(r)
