"""Checks for the Wigner-function section of the notes (hbar = I = 1).

Branch state u(phi) = e^{i(k0 phi - g sin phi)}/sqrt(2pi), c_{k0-j} = J_j(g).
Closed form: W(phi, l) = J_{2(k0-l)}(2 g cos phi) / 2pi,  l in Z/2.
Parity (even) state: W = [J_{2(k0-l)}(2g cos) + J_{2(k0+l)}(2g cos)]/4pi
                        + cos(2 k0 phi - pi l) J_{2l}(2 g sin phi)/2pi.
"""
import numpy as np
from scipy.special import jv

def W_from_coeffs(c, ks, phi):
    """W(phi,l) = (1/2pi) sum_{k+k'=2l} c_k c_k'^* e^{i(k-k')phi}; returns dict 2l -> array."""
    out = {}
    for a, k in enumerate(ks):
        for b, kp in enumerate(ks):
            if abs(c[a]) < 1e-15 or abs(c[b]) < 1e-15:
                continue
            key = k + kp
            out[key] = out.get(key, 0) + c[a] * np.conj(c[b]) * np.exp(1j * (k - kp) * phi) / (2 * np.pi)
    return out

g, k0 = 6.0, 40
phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)
dphi = phi[1] - phi[0]
js = np.arange(-30, 31)
# branch state
ks = k0 - js
c = jv(js, g).astype(complex)
W = W_from_coeffs(c, ks, phi)
err = 0
for key, val in W.items():
    l2 = key                       # 2l
    closed = jv(2 * k0 - l2, 2 * g * np.cos(phi)) / (2 * np.pi)
    err = max(err, np.abs(val - closed).max())
print('branch: max |W(sum) - J_{2(k0-l)}(2g cos phi)/2pi| =', err)
# marginals
mi, mh = 0, 0
for key, val in W.items():
    m = (val.real.sum() * dphi)
    if key % 2 == 0:
        j = k0 - key // 2
        mi = max(mi, abs(m - jv(j, g) ** 2))
    else:
        mh = max(mh, abs(m))
print('branch: integer-l marginal vs J_j^2: %.2e ; half-integer marginal: %.2e' % (mi, mh))
posm = sum(v.real for v in W.values())
print('branch: sum_l W - 1/2pi: %.2e' % np.abs(posm - 1 / (2 * np.pi)).max())
neg = sum(np.abs(v.real).sum() * dphi for v in W.values()) - 1
print('branch: negativity volume (g=%g) = %.3f' % (g, neg))
# Neumann integral J_n^2(z) = (1/pi) int_0^pi J_2n(2 z cos t) dt
t = np.linspace(0, np.pi, 20001)
for n in (0, 3, 7):
    val = np.trapezoid(jv(2 * n, 2 * g * np.cos(t)), t) / np.pi
    print('  Neumann: n=%d  %.10f vs J_n^2 = %.10f' % (n, val, jv(n, g) ** 2))
# phi = 0 cut vs Debye chord formula, large gamma
G = 200.0
for j in (20, 75, 130):
    beta = np.arccos(j / G)
    Phi = 2 * np.sqrt(G * G - j * j) - 2 * j * beta - np.pi / 2
    exact = jv(2 * j, 2 * G) / (2 * np.pi)
    debye = (1 / (2 * np.pi)) * np.sqrt(1 / (np.pi * G * np.sin(beta))) * np.cos(Phi + np.pi / 4)
    print('  W(0, k0-j) at gamma=200, j=%d: exact %.6e  Debye %.6e' % (j, exact, debye))
# parity (even) state
ks2 = np.concatenate([k0 - js, -(k0 - js)])
c2 = np.concatenate([jv(js, g), jv(js, g)]).astype(complex) / np.sqrt(2)
W2 = W_from_coeffs(c2, ks2, phi)
err2 = 0
for key, val in W2.items():
    l = key / 2
    closed = (jv(2 * k0 - key, 2 * g * np.cos(phi)) + jv(2 * k0 + key, 2 * g * np.cos(phi))) / (4 * np.pi) \
        + np.cos(2 * k0 * phi - np.pi * l) * jv(key, 2 * g * np.sin(phi)) / (2 * np.pi)
    err2 = max(err2, np.abs(val - closed).max())
print('even state: max |W(sum) - closed form| =', err2)
# reduced Wigner entropy = S_vN
p = np.concatenate([jv(js, g) ** 2, jv(js, g) ** 2]) / 2
p = p[p > 1e-300]
print('even state: S from reduced Wigner = S_vN = %.6f = ln2 + S_branch = %.6f'
      % (-np.sum(p * np.log(p)), np.log(2) - np.sum(jv(js, g) ** 2 * np.log(jv(js, g) ** 2))))
