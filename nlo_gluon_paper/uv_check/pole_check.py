r"""Group I Row III: the p+ integrands of your original row (Group_I.tex, before the p+ integral) have the denominator
p+ (x-z)^2 - k+ (x-w)^2, which vanishes at p0 = k+ (x-w)^2/(x-z)^2.  (1) residues of the five p+ integrands at p0
and their sum over the integrals whose range contains p0; (2) principal-value integral of the third integrand vs the
closed form with log|.|.   python3 pole_check.py -> pole_check.json"""
import json
import os
import numpy as np
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__))
K = 1.0
VEE = 50.0


def sq(v):
    return float(v @ v)


def pieces(p, P):
    """the five integrands (common factor i g^4 f^{abd}/(8 pi^4) and colour dropped), WITHOUT the pole factor;
    returns list of (h(p), sign) with integrand = h(p) / (sign*(p A - k B))"""
    x, y, z, w, xp, wp = (P[n] for n in ('x', 'y', 'z', 'w', 'xp', 'wp'))
    A, B = sq(x - z), sq(x - w)
    X1, Y, ZW = xp - wp, y - z, z - w
    k1 = (X1 @ Y) / (sq(X1) * sq(Y))
    h12 = k1 * (p + K) / (K - p)**2                                           # T1, T2: / (pA - kB)
    h3 = -(X1 @ (x - w)) * ((x - z) @ Y) / (sq(X1) * B * A * sq(Y)) * (K * B + p * A) / (p * K)   # / (pA - kB)
    V = K * (x - w) - p * (x - z)
    I2 = np.eye(2)

    def br(i, j, k, m):
        return I2[i, k] * I2[j, m] - (K / p) * I2[j, k] * I2[i, m] - (K / (K - p)) * I2[i, j] * I2[k, m]
    s4 = s5 = 0.0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for m in range(2):
                    c = X1[i] * Y[j] * ZW[k] * V[m] * br(i, j, k, m)
                    s4 += c / (sq(X1) * sq(Y) * sq(ZW))
                    s5 += c / (sq(X1) * sq(Y) * sq(V))
    h4 = -2.0 * (1 / K) * s4 / (K - p)          # 1/(4 pi^4) = 2/(8 pi^4);  / (kB - pA) = -/(pA - kB)
    h5 = -2.0 * s5 * p / (K - p)
    return A, B, [(h12, 1, 'T1', (K, VEE)), (h12, 1, 'T2', (0, K)), (h3, 1, 'T3', (0, VEE)),
                  (h4, -1, 'T4', (K, VEE)), (h5, -1, 'T5', (0, K))]


out = {'residue sums': []}
rng = np.random.default_rng(3)
for trial in range(8):
    P = {n: rng.normal(size=2) for n in ('x', 'y', 'z', 'w', 'xp', 'wp')}
    A, B, _ = pieces(0.5, P)
    p0 = K * B / A
    A, B, ps = pieces(p0, P)
    res = {name: sg * h / A for h, sg, name, (lo, hi) in ps if lo < p0 < hi}
    out['residue sums'].append({'p0/k+': p0, 'residues': res, 'sum': sum(res.values())})
    print('p0/k+ = %.3f' % p0, {k: '%.3e' % v for k, v in res.items()}, 'sum = %.3e' % sum(res.values()))

# (2) PV of T3's integral vs closed form
P = {n: rng.normal(size=2) for n in ('x', 'y', 'z', 'w', 'xp', 'wp')}
A, B = sq(P['x'] - P['z']), sq(P['x'] - P['w'])
L = 1e-3
f = lambda p: (K * B + p * A) / ((p * A - K * B) * p * K)
p0 = K * B / A
pv = integrate.quad(lambda p: (K * B + p * A) / (p * K * A), L, VEE, weight='cauchy', wvar=p0)[0]
closed_abs = (1 / K) * (-np.log(VEE / L) + 2 * np.log(abs((VEE * A - K * B) / (L * A - K * B))))
out['PV of third integrand'] = {'p0': p0, 'numerical PV': pv, 'closed form with log|.|': closed_abs,
                                'argument of your log (without |.|)': (K * B - VEE * A) / (K * B - L * A)}
print(out['PV of third integrand'])
json.dump(out, open(os.path.join(HERE, 'pole_check.json'), 'w'), indent=1, default=float)
