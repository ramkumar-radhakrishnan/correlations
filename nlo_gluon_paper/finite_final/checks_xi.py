"""Group III and Group I Row V in the variables of AABKL (arXiv:1411.2869):
  Group III:  xi = p+/(p+ + k+) in [0, 1-k+/vee],  y' = (1-xi) u + xi z   (u: measured gluon, z: unobserved, y': parent)
  Group I V:  xi = p+/k+ in [0,1]
Checks (output checks_xi.json):
 1. the p+ integrals of the plus prescriptions in xi:
    int_0^{vee-k+} dp+ [1/p+]_+^{(k+)} F  =  int_0^{1-x} dxi (g(xi)-g(0))/xi,  g = F(p+(xi))/(1-xi),  x = k+/vee
    ... + F(0) log(1-x)  =  int_0^{1-x} dxi g/(xi)_+     ((xi)_+ the usual plus distribution on [0,1])
    int dp+ [ -(1/k+) F3(p+) + F3(vee) k+/(p+ + k+) ]  =  -int_0^{1-x} dxi (g3(xi)-g3(1))/(1-xi),  g3 = k+ F3/(1-xi)
 2. pointwise, at random points: k+ x [p+ integrand at p+(xi), y'(u)] = [xi integrand at (xi,u)] for every term of
    Rows V (= III(1), VI: same kernel) and III part 2, with an arbitrary colour function C(y');
    the third term of Row V / the p+/k+^2 term of III(2) at xi = 1 equals the integrand of the (B) line."""
import json
import numpy as np
from scipy.integrate import quad

rng = np.random.default_rng(7)
out = {}


def cquad(f, a, b):
    re = quad(lambda t: f(t).real, a, b, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    im = quad(lambda t: f(t).imag, a, b, limit=400, epsabs=1e-13, epsrel=1e-12)[0]
    return re + 1j * im


# ---------------------------------------------------------------- 1. plus prescriptions
kp, V = 1.3, 40.0
x = kp / V
a, b = 0.7, 0.4
F = lambda p: np.exp(-1j * a * p / kp) / (1 + b * p / kp)          # F(p+), F(0) = 1
pof = lambda xi: kp * xi / (1 - xi)
lhs = cquad(lambda p: (F(p) - F(0) * kp / (p + kp)) / p if p > 0 else -1j * a / kp - b / kp + 1 / kp, 0, V - kp)
g = lambda xi: F(pof(xi)) / (1 - xi)
rhs = cquad(lambda xi: (g(xi) - g(0)) / xi if xi > 0 else 0j, 0, 1 - x)
out['[1/p+]_+^(k+) -> [1/xi]_+'] = abs(lhs - rhs)
std = cquad(lambda xi: (g(xi) - g(0)) / xi if xi > 0 else 0j, 0, 1 - x) + g(0) * np.log(1 - x)
out['... + F(0) log(1-x) - int g/(xi)_+ (by definition)'] = 0.0
# F3: integrand of the -(1/k+) term after int_{y'}; its large-p+ limit is F3(vee) k+/(p+ + k+)  -> model:
c3 = 0.8
F3 = lambda p: (kp / (p + kp)) * (1 + c3 * kp / (p + kp))          # -> k+/(p+ + k+) at large p+, g3(1) = 1
lhs3 = quad(lambda p: -F3(p) / kp + 1 / (p + kp), 0, V - kp, limit=400)[0]
g3 = lambda xi: kp * F3(pof(xi)) / kp / (1 - xi)                   # = (1 + c3 (1-xi))
rhs3 = -quad(lambda xi: (g3(xi) - 1.0) / (1 - xi), 0, 1 - x, limit=400)[0]
out['-(1/k+) F3 + (B) line -> -[1/(1-xi)]_+'] = abs(lhs3 - rhs3)
lim = [-quad(lambda xi: (g3(xi) - 1.0) / (1 - xi), 0, 1 - kp / VV, limit=400)[0] for VV in (1e2, 1e4, 1e6)]
out['[1/(1-xi)]_+ term for vee = 1e2, 1e4, 1e6 (finite)'] = lim
alone = [quad(lambda p: -F3(p) / kp, 0, VV - kp, limit=400)[0] for VV in (1e2, 1e4, 1e6)]
out['-(1/k+) term alone for vee = 1e2, 1e4, 1e6 (grows like -log)'] = alone

# ---------------------------------------------------------------- 2. pointwise maps
K = lambda r: r / (r @ r)
k = np.array([0.9, -0.4])
C = lambda yp: np.cos(yp @ np.array([0.3, 1.1])) + 0.5j * np.sin(yp @ np.array([-0.7, 0.2]))   # any colour function


def rowV_p(pp, yp, P, d):
    """the three terms of the Row V p+ integrand (without the common prefactor), [1/p+] term unsubtracted"""
    ph = np.exp(-1j * k @ (yp - P['w'])) * np.exp(-1j * (k @ (yp - P['z'])) * pp / kp)
    Kyz, Kxy, Kyw, Kxz = K(yp - P['z']), K(P['xp'] - yp), K(P['y'] - P['w']), K(P['x'] - P['z'])
    t1 = (Kyz @ Kxy) * (Kyw @ Kxz) / (pp + kp)
    t2 = -(Kyz @ Kxz) * (Kxy @ Kyw) / pp
    t3 = -(Kyz @ Kyw) * (Kxy @ Kxz) / kp
    return ph * C(yp) * np.array([t1, t2, t3])


def rowV_xi(xi, u, P):
    yp = (1 - xi) * u + xi * P['z']
    ph = np.exp(-1j * k @ (u - P['w']))
    Kuz, Kxy, Kyw, Kxz = K(u - P['z']), K(P['xp'] - yp), K(P['y'] - P['w']), K(P['x'] - P['z'])
    return ph * C(yp) * np.array([(Kuz @ Kxy) * (Kyw @ Kxz), -(Kuz @ Kxz) * (Kxy @ Kyw) / xi,
                                  -(Kuz @ Kyw) * (Kxy @ Kxz) / (1 - xi)])


def rowV_B(u, P):
    """integrand of the (B) line of Row V, y' renamed u (coefficient 1/(p+ + k+) -> 1/(1-xi))"""
    ph = np.exp(-1j * k @ (u - P['w']))
    return ph * C(P['z']) * (K(u - P['z']) @ K(P['y'] - P['w'])) * (K(P['xp'] - P['z']) @ K(P['x'] - P['z']))


def row32_p(pp, yp, P, d=2.0):
    ph = np.exp(-1j * k @ (yp - P['w'])) * np.exp(-1j * (k @ (yp - P['z'])) * pp / kp)
    m, kk = K(yp - P['z']), K(P['xp'] - yp)
    xz, xw, zw = P['x'] - P['z'], P['x'] - P['w'], P['z'] - P['w']
    D = pp * (xz @ xz) + kp * (xw @ xw)
    A1 = (m @ kk) * d * pp / (2 * (zw @ zw) * (pp + kp) ** 2) * ((xz @ xz) - (xw @ xw))
    A2 = pp * (m @ kk) / (kp * (pp + kp)) * ((xz @ zw) / (zw @ zw) - (xz @ xw) / (2 * (xw @ xw))
                                               - ((xz @ xz) - (xw @ xw)) / (2 * (zw @ zw)))
    A3 = (m @ kk) / (pp + kp) * ((xw @ zw) / (zw @ zw) + (xz @ xw) / (2 * (xz @ xz)) - ((xz @ xz) - (xw @ xw)) / (2 * (zw @ zw)))
    A4 = -1 / pp * ((kk @ xw) * (m @ zw) / (zw @ zw) + (m @ xz) * (kk @ xw) / (2 * (xz @ xz)))
    A5 = -pp / kp ** 2 * ((kk @ xz) * (m @ zw) / (zw @ zw) - (kk @ xz) * (m @ xw) / (2 * (xw @ xw)))
    A6 = -1 / kp * ((m @ xz) * (kk @ zw) / (zw @ zw) - (m @ xz) * (kk @ xw) / (2 * (xw @ xw))
                    + (m @ xw) * (kk @ zw) / (zw @ zw) + (kk @ xz) * (m @ xw) / (2 * (xz @ xz)))
    return ph * C(yp) / D * np.array([A1, A2, A3, A4, A5, A6])


def row32_xi(xi, u, P, d=2.0):
    yp = (1 - xi) * u + xi * P['z']
    ph = np.exp(-1j * k @ (u - P['w']))
    m, kk = K(u - P['z']), K(P['xp'] - yp)
    xz, xw, zw = P['x'] - P['z'], P['x'] - P['w'], P['z'] - P['w']
    Dx = xi * (xz @ xz) + (1 - xi) * (xw @ xw)
    A1 = (m @ kk) * d * xi * (1 - xi) / (2 * (zw @ zw)) * ((xz @ xz) - (xw @ xw))
    A2 = xi * (m @ kk) * ((xz @ zw) / (zw @ zw) - (xz @ xw) / (2 * (xw @ xw)) - ((xz @ xz) - (xw @ xw)) / (2 * (zw @ zw)))
    A3 = (1 - xi) * (m @ kk) * ((xw @ zw) / (zw @ zw) + (xz @ xw) / (2 * (xz @ xz)) - ((xz @ xz) - (xw @ xw)) / (2 * (zw @ zw)))
    A4 = -(1 - xi) / xi * ((kk @ xw) * (m @ zw) / (zw @ zw) + (m @ xz) * (kk @ xw) / (2 * (xz @ xz)))
    A5 = -xi / (1 - xi) * ((kk @ xz) * (m @ zw) / (zw @ zw) - (kk @ xz) * (m @ xw) / (2 * (xw @ xw)))
    A6 = -((m @ xz) * (kk @ zw) / (zw @ zw) - (m @ xz) * (kk @ xw) / (2 * (xw @ xw))
           + (m @ xw) * (kk @ zw) / (zw @ zw) + (kk @ xz) * (m @ xw) / (2 * (xz @ xz)))
    return ph * C(yp) / (kp * Dx) * np.array([A1, A2, A3, A4, A5, A6])


def row32_B(u, P):
    ph = np.exp(-1j * k @ (u - P['w']))
    xz, xw, zw = P['x'] - P['z'], P['x'] - P['w'], P['z'] - P['w']
    m = K(u - P['z'])
    return -ph * C(P['z']) * (K(P['xp'] - P['z']) @ K(xz)) * ((m @ zw) / (zw @ zw) - (m @ xw) / (2 * (xw @ xw))) / kp


errV, err32 = 0.0, 0.0
eps = (1e-4, 1e-6, 1e-8)
errB = {e: 0.0 for e in eps}
for _ in range(200):
    P = {n: rng.normal(size=2) for n in ('w', 'z', 'y', 'x', 'xp')}
    xi, u = rng.uniform(0.02, 0.98), rng.normal(size=2)
    pp, yp = pof(xi), (1 - xi) * u + xi * P['z']
    a_ = kp * rowV_p(pp, yp, P, 2)          # dp+ d^2y' = k+ dxi d^2u (Jacobians (1-xi)^-2 and (1-xi)^2)
    errV = max(errV, np.max(np.abs(a_ - rowV_xi(xi, u, P))) / np.max(np.abs(a_)))
    b_ = kp * row32_p(pp, yp, P)
    err32 = max(err32, np.max(np.abs(b_ - row32_xi(xi, u, P))) / np.max(np.abs(b_)))
    # xi -> 1: (1-xi) x [third term of V] and (1-xi) x [A5 term of III(2)] tend to minus the (B) integrands
    for e1 in eps:
        t3 = e1 * rowV_xi(1 - e1, u, P)[2]
        t5 = e1 * row32_xi(1 - e1, u, P)[4]
        errB[e1] = max(errB[e1], abs(t3 + rowV_B(u, P)) / abs(rowV_B(u, P)),
                       abs(t5 - row32_B(u, P)) / abs(row32_B(u, P)))
out['Row V (III(1), VI): k+ x p+ integrand = xi integrand, max rel. diff'] = errV
out['Row III(2): k+ x p+ integrand = xi integrand, max rel. diff'] = err32
out['xi -> 1: (1-xi) x term = (B) integrand, max rel. diff at 1-xi = 1e-4, 1e-6, 1e-8'] = [errB[e] for e in eps]

if __name__ == '__main__':
    for kk_, vv in out.items():
        print(kk_, vv)
    json.dump(out, open('checks_xi.json', 'w'), indent=1, default=float)
