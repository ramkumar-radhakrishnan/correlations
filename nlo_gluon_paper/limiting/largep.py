r"""Large-p+ behaviour of Group III Rows III, V, VI (upper limit vee - k+ -> infinity).
At large p+ the y' integral is a Fourier transform at q = k (1 + p+/k+) of (y'-z)^m/(y'-z)^2 g(y'), where g is
smooth once x' is integrated (charges and Wilson lines are smooth).  Line 1 has it with the weights 1/(p+ + k+) and
[k+/(p+(p+ + k+))]_+  ->  integrand ~ (1/p+) (1/q) ~ 1/p+^2.  Line 2 (weight 1/k+) has g(y') - g(z) (the x' -> z, y' -> z
subtraction), whose transform is O(1/q^2)  ->  integrand ~ 1/p+^2.  Numerical check with a smooth model g:
q |I_unsub| -> 2 pi |g(0)| and q^2 |I_sub| -> const.   python3 largep.py -> largep.json"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N, L = 4096, 12.0
dx = L / N
x = (np.arange(N) - N // 2 + 0.5) * dx          # offset grid: r = 0 is not a node
X, Y = np.meshgrid(x, x, indexing='ij')
R2 = X**2 + Y**2
c = np.array([0.4, -0.3])
g = np.exp(-((X - c[0])**2 + (Y - c[1])**2)) * (1 + 0.3 * X)
g0 = np.exp(-(c**2).sum())
qgrid = 2 * np.pi * np.fft.fftfreq(N, d=dx)
out = {}
for name, h in (('unsub', g), ('sub', g - g0)):
    f = X / R2 * h * np.exp(-R2 / 6.0)            # m = 1 component; wide Gaussian damps the far region of g0
    F = np.fft.fft2(np.fft.ifftshift(f)) * dx * dx
    rows = []
    for qv in (5, 10, 20, 40, 80):
        i = int(round(qv / (qgrid[1])))            # q along the x axis
        q = qgrid[i]
        val = abs(F[i, 0])
        rows.append({'q': float(q), '|I|': float(val), 'q|I|': float(q * val), 'q^2|I|': float(q * q * val)})
    out[name] = rows
    print(name, ['q=%.0f q|I|=%.4f q^2|I|=%.3f' % (r['q'], r['q|I|'], r['q^2|I|']) for r in rows])
print('2 pi g(0) =', 2 * np.pi * g0)
out['2pi g0'] = 2 * np.pi * g0
json.dump(out, open(os.path.join(HERE, 'largep.json'), 'w'), indent=1)
