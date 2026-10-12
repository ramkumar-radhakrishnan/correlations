"""Checks for the Fourier transform of Group I, Row II (Part I of the combined note).
1. partial fractions of A^i/(k^2 p^2)  (sympy)
2. A + B = 0 with the complete bracket, = k^i with the bracket of the message  (sympy)
3. bracket = A log(vee/Lambda) + B log(vee/k+) + R + O(k+/vee)  (numbers)
4. alpha term = [N2, A1]: [V, rho^a(u)] identity with explicit charge matrices (SU(2), 3 sites)
5. main.tex xi-integral vs the two log forms"""
import json
import os
import numpy as np
import sympy as sp
from scipy.integrate import quad

out = {}
k1, k2, p1, p2 = sp.symbols('k1 k2 p1 p2', real=True)
k = sp.Matrix([k1, k2]); p = sp.Matrix([p1, p2]); q = k - p
dot = lambda a, b: (a.T * b)[0]
K2, P2, Q2, KP = dot(k, k), dot(p, p), dot(q, q), dot(k, p)
E = (K2 * p - KP * k) / Q2
C = ((2 * P2 - 2 * KP + K2) * k + K2 * p) / Q2
D = (P2 * k + K2 * p) / Q2
A = -E - 2 * k - D
# 1. partial fractions
lhs = A / (K2 * P2)
rhs = -2 * p / (Q2 * P2) + dot(p, q) * k / (Q2 * K2 * P2) - 2 * k / (K2 * P2)
out['partial_fractions'] = all(sp.simplify(x) == 0 for x in (lhs - rhs))
out['E_plus_D'] = all(sp.simplify(x) == 0 for x in (E + D - (2 * K2 * p - dot(p, q) * k) / Q2))
# 2. coefficient of log vee
B_msg = E + 2 * k + C
B_full = B_msg - k
out['A_plus_B_full'] = [str(sp.simplify(x)) for x in (A + B_full)]
out['A_plus_B_msg'] = [str(sp.simplify(x)) for x in (A + B_msg)]
out['C_minus_D'] = [str(sp.simplify(x)) for x in (C - D)]
# 3. large-vee expansion of the complete bracket
rng = np.random.default_rng(3)
num = []
for trial in range(3):
    kk, pp = rng.normal(size=2), rng.normal(size=2)
    qq = kk - pp
    kk2, pp2, qq2, kp = kk @ kk, pp @ pp, qq @ qq, kk @ pp
    En = (kk2 * pp - kp * kk) / qq2
    Cn = ((2 * pp2 - 2 * kp + kk2) * kk + kk2 * pp) / qq2
    Dn = (pp2 * kk + kk2 * pp) / qq2
    An = -En - 2 * kk - Dn
    Bn = En + kk + Cn
    Rn = -2 * kk * np.log(pp2 / kk2) + Cn * np.log(kk2 / qq2)
    for V in (1e2, 1e4, 1e6):
        kp_, Lam = 1.0, 1e-3
        LL, LK = np.log(V / Lam), np.log(V / kp_)
        full = (-LL + LK) * En - 2 * kk * LL - 2 * kk * np.log(kp_ * pp2 / (V * kk2 + kp_ * pp2)) \
            + Cn * np.log(1 + (V / kp_ - 1) * kk2 / qq2) - Dn * (LL + np.log(1 - kp_ / V)) - kk * LK
        approx = An * LL + Bn * LK + Rn
        num.append(dict(trial=trial, vee=V, diff=float(np.abs(full - approx).max()), size=float(np.abs(full).max())))
out['large_vee'] = num
# 4. commutator identity with explicit charges: sites i, SU(2) spin-1/2 on each site
sx = np.array([[0, 1], [1, 0]], complex) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.array([[1, 0], [0, -1]], complex) / 2
T = [sx, sy, sz]
eps = np.zeros((3, 3, 3))
for a, b, c, s in ((0, 1, 2, 1), (1, 2, 0, 1), (2, 0, 1, 1), (1, 0, 2, -1), (0, 2, 1, -1), (2, 1, 0, -1)):
    eps[a, b, c] = s
ns = 3
I2 = np.eye(2)


def site_op(m, i):
    ops = [I2] * ns
    ops[i] = m
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


rho = [[site_op(T[a], i) for a in range(3)] for i in range(ns)]
# check algebra [rho^a(i), rho^b(j)] = i f^{abc} rho^c(i) delta_ij
alg = max(np.abs(rho[i][a] @ rho[j][b] - rho[j][b] @ rho[i][a] - (1j * sum(eps[a, b, c] * rho[i][c] for c in range(3)) if i == j else 0)).max()
          for i in range(ns) for j in range(ns) for a in range(3) for b in range(3))
Kij = rng.normal(size=(ns, ns)); Kij = Kij + Kij.T          # symmetric kernel K(x,y,z) integrated over z
cK = rng.normal(size=ns)                                     # emission kernel K^i(u-w), one component
V = sum(Kij[i, j] * rho[i][b] @ rho[j][b] for i in range(ns) for j in range(ns) for b in range(3))
worst = 0
for a in range(3):
    Aop = sum(cK[l] * rho[l][a] for l in range(ns))
    lhs = V @ Aop - Aop @ V
    rhs = -1j * sum(Kij[i, j] * cK[j] * eps[a, b, c] * (rho[i][b] @ rho[j][c] + rho[j][c] @ rho[i][b])
                    for i in range(ns) for j in range(ns) for b in range(3) for c in range(3))
    # odd projection of the kernel K_ij cK_j under i<->j
    rhs2 = -0.5j * sum(Kij[i, j] * (cK[j] - cK[i]) * eps[a, b, c] * (rho[i][b] @ rho[j][c] + rho[j][c] @ rho[i][b])
                       for i in range(ns) for j in range(ns) for b in range(3) for c in range(3))
    worst = max(worst, np.abs(lhs - rhs).max(), np.abs(lhs - rhs2).max())
out['commutator'] = dict(algebra=float(alg), identity=float(worst), size=float(np.abs(V).max()))
# 5. main.tex xi integral vs the two log forms
rng = np.random.default_rng(7)
kp_, L = 1.0, 1e-4
xi = []
for V in (1e2, 1e4):
    for trial in range(2):
        kk, pp = rng.normal(size=2), rng.normal(size=2)
        qq = kk - pp
        kk2, pp2, qq2, kp = kk @ kk, pp @ pp, qq @ qq, kk @ pp
        c = 0
        En = (kk2 * pp[c] - kp * kk[c]) / qq2
        Cn = ((2 * pp2 - 2 * kp + kk2) * kk[c] + kk2 * pp[c]) / qq2
        Dn = (pp2 * kk[c] + kk2 * pp[c]) / qq2
        f = lambda x: x / ((x * qq2 + (1 - x) * kk2) * (1 - x)) * ((2 / x * (pp2 - kp) - pp2 + 2 * kp) * kk[c] + kk2 * pp[c] / x)
        I = quad(f, L / kp_, 1 - L / kp_, limit=800, points=[0.5, 0.9, 0.99])[0]
        main = np.log(L / kp_) * En - 2 * kk[c] * np.log(V / L) - 2 * kk[c] * np.log(kp_ * pp2 / (V * kk2 + kp_ * pp2)) - I
        msg = (-np.log(V / L) + np.log(V / kp_)) * En - 2 * kk[c] * np.log(V / L) - 2 * kk[c] * np.log(kp_ * pp2 / (V * kk2 + kp_ * pp2)) \
            + Cn * np.log(1 + (V / kp_ - 1) * kk2 / qq2) - Dn * (np.log(V / L) + np.log(1 - kp_ / V))
        xi.append(dict(vee=V, main=float(main), full=float(msg - kk[c] * np.log(V / kp_)), msg=float(msg)))
out['xi'] = xi
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ft_checks.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
