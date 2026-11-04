r"""Group II Row IV from the amplitude: checks.
1. Polarization sum  sum_{ij} T^m_{ijk} T^{m'}_{ijk'}  (parent k, measured i, unmeasured j), any d,
   T^m_{ijk} = d_{km} d_{ij} - alpha d_{jm} d_{ik} - beta d_{im} d_{jk},  alpha = P/p = 1/(1-xi), beta = P/k = 1/xi
   = (alpha^2+beta^2) S1 + 2 alpha beta S2 + (d - 2 alpha - 2 beta) S3.
2. xi(1-xi) x coefficients: n1 = (1-xi)/xi + xi/(1-xi), n2 = 2, n3 = d xi(1-xi) - 2.
3. Colour: f^{c'da} f^{cda} = N_c delta^{cc'}.
4. Normalization: 4/(2pi)^3 * (1/2pi) * g^4 N_c/(16 pi^4) = (1/(2pi)^3) g^4 N_c/(8 pi^5).
5. Traceless part of the z integral in d=2: J_{mm'} = A delta + B' Xhat Xhat with B' = -pi.
python3 checks.py -> checks.json"""
import json
import itertools
import numpy as np
import sympy as sp
from scipy.integrate import dblquad

out = {}
# 1. polarization sum, numerically for d = 2,3,4
rng = np.random.default_rng(1)
res = []
for d in (2, 3, 4):
    xi = 0.37
    al, be = 1 / (1 - xi), 1 / xi
    I = np.eye(d)
    T = np.zeros((d, d, d, d))           # T[m,i,j,k]
    for m, i, j, k in itertools.product(range(d), repeat=4):
        T[m, i, j, k] = I[k, m] * I[i, j] - al * I[j, m] * I[i, k] - be * I[i, m] * I[j, k]
    lhs = np.einsum('mijk,nijl->mnkl', T, T)            # (m, m', k, k')
    S1 = np.einsum('mn,kl->mnkl', I, I)
    S2 = np.einsum('ml,kn->mnkl', I, I)
    S3 = np.einsum('mk,nl->mnkl', I, I)
    rhs = (al**2 + be**2) * S1 + 2 * al * be * S2 + (d - 2 * al - 2 * be) * S3
    res.append(float(np.abs(lhs - rhs).max()))
out['polsum_maxdiff_d234'] = res
print('1. polarization sum, max |lhs-rhs| for d=2,3,4:', res)

# 2. coefficients times xi(1-xi)
x, d = sp.symbols('xi d', positive=True)
al, be = 1 / (1 - x), 1 / x
n1 = sp.simplify(x * (1 - x) * (al**2 + be**2) - ((1 - x) / x + x / (1 - x)))
n2 = sp.simplify(x * (1 - x) * 2 * al * be - 2)
n3 = sp.simplify(x * (1 - x) * (d - 2 * al - 2 * be) - (d * x * (1 - x) - 2))
out['n_checks'] = [str(n1), str(n2), str(n3)]
print('2. n1, n2, n3 residuals:', n1, n2, n3)

# 3. colour f^{c'da} f^{cda} = N_c delta  (SU(3))
lam = [np.array(m, dtype=complex) for m in (
    [[0, 1, 0], [1, 0, 0], [0, 0, 0]], [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
    [[0, 0, 1], [0, 0, 0], [1, 0, 0]], [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], np.diag([1, 1, -2]) / np.sqrt(3))]
t = [l / 2 for l in lam]
f = np.zeros((8, 8, 8))
for a, b, c in itertools.product(range(8), repeat=3):
    f[a, b, c] = (-2j * np.trace((t[a] @ t[b] - t[b] @ t[a]) @ t[c])).real
col = np.einsum('xda,cda->xc', f, f)
out['colour_maxdiff'] = float(np.abs(col - 3 * np.eye(8)).max())
print('3. f f - N_c delta:', out['colour_maxdiff'])

# 4. normalization
g, Nc = sp.symbols('g N_c', positive=True)
pref = sp.Rational(4) / (2 * sp.pi)**3 * 1 / (2 * sp.pi) * g**4 * Nc / (16 * sp.pi**4)
out['norm_residual'] = str(sp.simplify(pref - g**4 * Nc / (8 * sp.pi**5) / (2 * sp.pi)**3))
print('4. normalization residual:', out['norm_residual'])

# 5. traceless part of J in d=2
X = np.array([1.0, 0.0])
def fz(r, th):
    z = np.array([r * np.cos(th), r * np.sin(th)]); u = z - X
    return r * ((z @ X) * (u @ X) - 0.5 * (z @ u)) / ((z @ z) * (u @ u))
v = sum(dblquad(fz, 0, 2 * np.pi, a, b, epsabs=1e-10)[0]
        for a, b in [(0, 0.5), (0.5, 0.99), (0.99, 1.01), (1.01, 2), (2, 50), (50, 5000)])
out['Bprime'] = 2 * v
print('5. B\' =', 2 * v, ' (-pi =', -np.pi, ')')
json.dump(out, open('checks.json', 'w'), indent=1)
