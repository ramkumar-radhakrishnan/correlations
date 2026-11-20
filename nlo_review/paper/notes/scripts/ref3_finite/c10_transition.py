# (C10) transition region: measured-gluon Fourier transform of the two-gluon amplitudes of Rows II.I (2nd part, A_B),
# II.IV (A_C = C1 A1) and of their difference (Row III.III, 2nd part, is -2Re A_C^+ A_B), at fixed t = p+/k+,
#     At(k, z) = int d^2w e^{-i k.w} A(w, z),   z = Z n,
# computed over the WHOLE w plane (all of w near the projectile, w near z, w near -t z, and |w| ~ |w - z| ~ Z, i.e. the
# transition region), for growing Z.  Prediction from the inner regions only:
#     Z At = c_soft + c_joint e^{-ik.z} + c_coll e^{+itk.z} + o(1)
# (c_soft: w near the projectile, A_B only; c_joint: w - z fixed, the joint limit of fix_tail; c_coll: parent at the
# projectile, A_C only, the collinear region).  If the region |w| ~ |w-z| ~ Z contributed at O(1/Z), the residual
# Z At - prediction would not vanish.  Coherent amplitude A_B - A_C: c_joint cancels (the parent pole).
# Amplitudes in coordinate space as in notes/fix_tail.tex (t1_density.py), x = 0 (one charge), colour stripped
# (uniform far target), k = (1,0), k+ = 1.
import sys, numpy as np
t = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
xi = 1 / (1 + t); NB = np.sqrt(t) / (1 + t); NCn = np.sqrt(t) / (1 + t) ** 2
k = np.array([1.0, 0.0])
def tC(kk, m, i, j):
    d = lambda a, b: 1.0 if a == b else 0.0
    return d(kk, m) * d(i, j) - (1 + t) / t * d(j, m) * d(i, kk) - (1 + t) * d(i, m) * d(kk, j)
TC = np.array([[[[tC(kk, m, i, j) for j in range(2)] for i in range(2)] for m in range(2)] for kk in range(2)])
def A_B(W1, W2, z):
    # x = 0: xz = -z, xw = -w, zw = z - w
    xw = (-W1, -W2); xz = (-z[0], -z[1]); zw = (z[0] - W1, z[1] - W2)
    A = W1 ** 2 + W2 ** 2; B = z @ z; C = zw[0] ** 2 + zw[1] ** 2
    D = t * B + A
    out = np.empty((2, 2) + W1.shape)
    for i in range(2):
        for j in range(2):
            T = ((1.0 if i == j else 0.0) * (B - A) / (2 * C)
                 + (1 + t) * (xz[j] * zw[i] / C - xz[j] * xw[i] / (2 * A))
                 + (1 + t) / t * (xw[i] * zw[j] / C + xz[j] * xw[i] / (2 * B)))
            out[i, j] = NB * T / D
    return out
def A_C(W1, W2, z):
    y1 = z[0] + xi * (W1 - z[0]); y2 = z[1] + xi * (W2 - z[1])
    K1 = (-y1 / (y1 ** 2 + y2 ** 2), -y2 / (y1 ** 2 + y2 ** 2))          # K(x - y), x = 0
    a1, a2 = W1 - z[0], W2 - z[1]; K2 = (a1 / (a1 ** 2 + a2 ** 2), a2 / (a1 ** 2 + a2 ** 2))   # K(w - z)
    out = np.zeros((2, 2) + W1.shape)
    for i in range(2):
        for j in range(2):
            for m in range(2):
                for kk in range(2):
                    if TC[kk, m, i, j] != 0: out[i, j] += K2[m] * K1[kk] * TC[kk, m, i, j]
    return NCn * out
# singular parts c^{ij l} (w - s)^l/|w - s|^2 e^{-|w-s|^2/(2 sig^2)} at each singular point s (exact 1/r coefficients)
def sing_parts(z, which):
    Z2 = z @ z; parts = []
    if which in ('B', 'coh'):
        sgn = 1.0
        c = np.zeros((2, 2, 2)); D0 = t * Z2
        for i in range(2):
            for j in range(2): c[i, j, i] += -(1 + t) * z[j] / 2 / D0
        parts.append((np.zeros(2), NB * c * sgn))
        c = np.zeros((2, 2, 2)); D1 = (1 + t) * Z2
        for i in range(2):
            for j in range(2):
                for l in range(2):
                    c[i, j, l] = (-(1.0 if i == j else 0.0) * z[l] + (1 + t) * z[j] * (1.0 if i == l else 0.0)
                                  + (1 + t) / t * z[i] * (1.0 if j == l else 0.0)) / D1
        parts.append((z.copy(), NB * c * sgn))
    if which in ('C', 'coh'):
        sgn = -1.0 if which == 'coh' else 1.0
        c = np.zeros((2, 2, 2))
        for i in range(2):
            for j in range(2):
                for m in range(2):
                    for kk in range(2): c[i, j, m] += -z[kk] / Z2 * TC[kk, m, i, j]
        parts.append((z.copy(), NCn * c * sgn))
        c = np.zeros((2, 2, 2))
        for i in range(2):
            for j in range(2):
                for m in range(2):
                    for kk in range(2): c[i, j, kk] += z[m] / Z2 * TC[kk, m, i, j]
        parts.append((-t * z, NCn * c * sgn))
    return parts
def At_numeric(z, which, h=0.1, Lfac=3.0, sig=1.0):
    Z = np.sqrt(z @ z); L = Lfac * Z; Lg = 1.6 * L
    n = int(2 * Lg / h); off = np.array([0.37, 0.61]) * h
    xs = -Lg + off[0] + h * np.arange(n); ys = -Lg + off[1] + h * np.arange(n)
    parts = sing_parts(z, which)
    res = np.zeros((2, 2), complex)
    for i0 in range(0, n, 400):
        W1, W2 = np.meshgrid(xs, ys[i0:i0 + 400], indexing='xy')
        if which == 'B': Am = A_B(W1, W2, z)
        elif which == 'C': Am = A_C(W1, W2, z)
        else: Am = A_B(W1, W2, z) - A_C(W1, W2, z)
        for s, c in parts:
            d1, d2 = W1 - s[0], W2 - s[1]; r2 = d1 ** 2 + d2 ** 2; g = np.exp(-r2 / (2 * sig ** 2)) / r2
            Am = Am - (c[:, :, 0, None, None] * d1 + c[:, :, 1, None, None] * d2) * g
        win = np.exp(-((W1 ** 2 + W2 ** 2) / L ** 2) ** 2)
        ph = np.exp(-1j * (k[0] * W1 + k[1] * W2)) * win
        res += np.einsum('ijab,ab->ij', Am, ph) * h * h
    # analytic FT of the subtracted singular parts: -2 pi i k^l/k^2 (1 - e^{-k^2 sig^2/2}) e^{-i k.s}
    for s, c in parts:
        res += np.einsum('ijl,l->ij', c, -2j * np.pi * k / (k @ k)) * (1 - np.exp(-(k @ k) * sig ** 2 / 2)) * np.exp(-1j * k @ s)
    return res
def predicted(z, which):
    Z = np.sqrt(z @ z); nv = z / Z; out = np.zeros((2, 2), complex); ft = -2j * np.pi * k / (k @ k)
    cs = np.zeros((2, 2), complex); cj = np.zeros((2, 2), complex); cc = np.zeros((2, 2), complex)
    for i in range(2):
        for j in range(2):
            cs[i, j] = NB * (-(1 + t) * nv[j] / (2 * t)) * ft[i]
            cj[i, j] = NB / (1 + t) * sum((-(1.0 if i == j else 0.0) * nv[l] + (1 + t) * nv[j] * (1.0 if i == l else 0.0)
                                            + (1 + t) / t * nv[i] * (1.0 if j == l else 0.0)) * ft[l] for l in range(2))
    cjC = np.zeros((2, 2), complex); ccC = np.zeros((2, 2), complex)
    for i in range(2):
        for j in range(2):
            cjC[i, j] = NCn * sum(-nv[kk] * TC[kk, m, i, j] * ft[m] for m in range(2) for kk in range(2))
            ccC[i, j] = NCn * sum(nv[m] * TC[kk, m, i, j] * ft[kk] for m in range(2) for kk in range(2))
    ej = np.exp(-1j * k @ z); ec = np.exp(1j * t * k @ z)
    if which == 'B': return cs + cj * ej, cj
    if which == 'C': return cjC * ej + ccC * ec, cjC
    return cs + (cj - cjC) * ej - ccC * ec, cj - cjC
if __name__ == '__main__':
    phi = 0.6; nv = np.array([np.cos(phi), np.sin(phi)])
    print('t = %.2f, direction phi = %.2f' % (t, phi))
    for which in ('B', 'C', 'coh'):
        for Z in (10.0, 20.0, 40.0, 80.0):
            z = Z * nv
            num = At_numeric(z, which); pred, cj = predicted(z, which)
            num2 = At_numeric(z, which, h=0.08, Lfac=4.0) if Z <= 20 else None
            print('%-4s Z=%5.0f  max|Z At| = %.4f   |c_joint| = %.4f   max|Z At - prediction| = %.2e%s' % (
                which, Z, np.abs(Z * num).max(), np.abs(cj).max(), np.abs(Z * num - pred).max(),
                '' if num2 is None else '   (finer grid, larger window: %.2e)' % np.abs(Z * num2 - pred).max()), flush=True)
