r"""The small-xi limit that produces log(vee/k+) H4.
At fixed x, x', y the y' integral of the DGLAP integrand is  I(xi) = int d^2y' e^{-ik.(y'-y)/xi} log((y'-y)^2) F(y'),
F smooth (kernel x Wilson lines x charges, x' integrated).  Claim: I(xi) = -(4 pi xi^2/k^2) F(y) + O(xi^4),
so (1-xi)/xi^3 [I(xi) + 4 pi xi^2 F(y)/k^2] is integrable at xi -> 0 (the subtracted DGLAP integrand) while
(1-xi)/xi^3 I(xi) ~ -(4pi/k^2) F(y)/xi gives the log(vee/k+).   python3 smallxi.py -> smallxi.json"""
import json
import numpy as np

N, L = 4096, 7.0
dx = L / N
x = (np.arange(N) - N // 2 + 0.5) * dx
X, Y = np.meshgrid(x, x, indexing='ij')
R2 = X**2 + Y**2
c = np.array([0.5, -0.2])
F = np.exp(-((X - c[0])**2 + (Y - c[1])**2)) * (1 + 0.4 * Y)          # smooth model, y = 0
F0 = np.exp(-(c**2).sum()) * 1.0
f = np.log(R2) * F
FT = np.fft.fft2(np.fft.ifftshift(f)) * dx * dx
qg = 2 * np.pi * np.fft.fftfreq(N, d=dx)
k = 2.0
rows = []
for xi in (0.5, 0.2, 0.1, 0.05, 0.03):
    qv = k / xi
    i = int(round(qv / qg[1]))
    q = qg[i]
    I = FT[i, 0] * np.exp(-1j * q * x[N // 2])          # the grid origin sits at +dx/2: undo that phase
    lead = -4 * np.pi / q**2 * F0
    pred = 1 - 2j * (2 * c[0]) / q                     # next order: FT[log r^2 r_j] = i d_j(-4pi/q^2) -> -2i q.gradF/(q^2 F)
    rows.append({'pred_next_order': [pred.real, pred.imag], 'xi': xi, 'q': float(q), 'I/lead': complex(I / lead).real, 'Im(I/lead)': complex(I / lead).imag,
                 '(I-lead)/lead': abs((I - lead) / lead)})
    print('xi=%.3f  I/leading = %.6f %+.6fi  (1 - i xi: %+.3fi)  |I-leading|/|leading| = %.2e' % (xi, (I / lead).real, (I / lead).imag, pred.imag, abs((I - lead) / lead)))
json.dump(rows, open('smallxi.json', 'w'), indent=1)
