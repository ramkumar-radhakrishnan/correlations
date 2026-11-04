r"""Negative logs with the i*epsilon prescription.
(a) int_{p1}^{p2} dp A/(p A - k B + i eps)  ->  log|(p2 A - k B)/(p1 A - k B)| - i pi theta(p0 in (p1,p2)),  p0 = k B/A.
(b) colour: Tr(T^a T^b T^c) + Tr(T^c T^b T^a) is the symmetric (d^{abc}) part only:
    = 0 for SU(2); = d^{abc}/2 for SU(3) fundamental (odderon-type), so it vanishes for a Gaussian (C-even) projectile.
python3 check.py -> check.json"""
import json
import itertools
import numpy as np

out = {'pplus': [], 'colour': {}}
k = 1.0
for A, B, p1, p2 in [(0.7, 1.9, 0.0, 10.0), (0.7, 1.9, 1.0, 10.0), (1.3, 0.4, 0.0, 1.0), (0.5, 0.2, 1.0, 10.0), (2.0, 0.3, 1.0, 50.0)]:
    p0 = k * B / A
    for eps in (1e-3, 1e-6):
        # exact antiderivative of A/(pA - kB + i eps): principal branch log is continuous along the path (Im > 0)
        val = np.log(p2 * A - k * B + 1j * eps) - np.log(p1 * A - k * B + 1j * eps)
        pred = np.log(abs((p2 * A - k * B) / (p1 * A - k * B))) - 1j * np.pi * (p1 < p0 < p2)
        out['pplus'].append({'A': A, 'B': B, 'p1': p1, 'p2': p2, 'p0': p0, 'eps': eps,
                             'Re': val.real, 'Im/pi': val.imag / np.pi, 'pred Re': pred.real, 'pred Im/pi': pred.imag / np.pi})
for r in out['pplus']:
    print({k_: (round(v, 10) if isinstance(v, float) else v) for k_, v in r.items()})
# numerical quadrature check of one case with finite eps
from scipy.integrate import quad
A, B, p1, p2, eps = 0.7, 1.9, 0.0, 10.0, 1e-4
re = quad(lambda p: A * (p * A - k * B) / ((p * A - k * B)**2 + eps**2), p1, p2, points=[k * B / A], limit=2000)[0]
im = quad(lambda p: -A * eps / ((p * A - k * B)**2 + eps**2), p1, p2, points=[k * B / A], limit=2000)[0]
out['quad'] = {'Re': re, 'Im/pi': im / np.pi}
print('quadrature (eps=1e-4):', re, im / np.pi)

# (b) colour
s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
T2 = [m / 2 for m in s]
lam = [np.array(m, dtype=complex) for m in (
    [[0, 1, 0], [1, 0, 0], [0, 0, 0]], [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
    [[0, 0, 1], [0, 0, 0], [1, 0, 0]], [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
    [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], np.diag([1, 1, -2]) / np.sqrt(3))]
T3 = [m / 2 for m in lam]
def sym_rev(T):
    n = len(T)
    return max(abs(np.trace(T[a] @ T[b] @ T[c] + T[c] @ T[b] @ T[a])) for a, b, c in itertools.product(range(n), repeat=3))
def antisym_rev(T):
    n = len(T)
    return max(abs(np.trace(T[a] @ T[b] @ T[c] - T[c] @ T[b] @ T[a])) for a, b, c in itertools.product(range(n), repeat=3))
d3 = np.zeros((8, 8, 8))
for a, b, c in itertools.product(range(8), repeat=3):
    d3[a, b, c] = 2 * np.trace((T3[a] @ T3[b] + T3[b] @ T3[a]) @ T3[c]).real
dev = max(abs(np.trace(T3[a] @ T3[b] @ T3[c] + T3[c] @ T3[b] @ T3[a]) - d3[a, b, c] / 2)
          for a, b, c in itertools.product(range(8), repeat=3))
out['colour'] = {'SU2 max|Tr(abc)+Tr(cba)|': float(sym_rev(T2)), 'SU2 max|Tr(abc)-Tr(cba)|': float(antisym_rev(T2)),
                 'SU3 max|Tr(abc)+Tr(cba)|': float(sym_rev(T3)), 'SU3 max|Tr(abc)+Tr(cba) - d/2|': float(dev)}
print(out['colour'])
json.dump(out, open('check.json', 'w'), indent=1, default=float)
