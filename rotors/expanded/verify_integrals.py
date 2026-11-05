import numpy as np
from scipy.integrate import quad
from scipy.special import ellipk
from rotorlib import *

print('int_0^pi ln sin u du =', quad(lambda u: np.log(np.sin(u)), 0, np.pi)[0],
      '   -pi ln2 =', -np.pi * np.log(2))
g = quad(lambda x: (1 + np.cos(x)) * np.log(1 + np.cos(x) + 1e-300), 0, 2 * np.pi,
         points=[np.pi], limit=400)[0] / (2 * np.pi)
print('<(1+cos)ln(1+cos)> =', g, '   1 - ln2 =', 1 - np.log(2))

# normalisation integrals (hbar=I=1): int_band dp / sqrt(lam^2-(E-p^2)^2), p>0
for E, lam in ((300., 100.), (60., 100.)):
    pmin, pmax = np.sqrt(max(E - lam, 0)), np.sqrt(E + lam)
    num = quad(lambda p: 1 / np.sqrt(max(lam ** 2 - (E - p * p) ** 2, 1e-300)),
               pmin, pmax, points=[pmin, pmax], limit=400)[0]
    if E > lam:
        k2 = 2 * lam / (E + lam)
        th = 0.5 * (2 / np.sqrt(E + lam)) * ellipk(k2)
    else:
        k2 = (E + lam) / (2 * lam)
        th = 0.5 * np.sqrt(2 / lam) * ellipk(k2)
    print('E=%5.0f lam=%5.0f : band integral numeric %.8f  elliptic %.8f'
          % (E, lam, num, th))

# free-rotation differential entropy vs closed form ln(pi lam/(2 sqrt E))
print('free rotation: differential entropy vs ln[pi lam /(2 sqrt E)]')
for lam in (200.,):
    for r in (5., 20., 80., 320.):
        E = r * lam
        pmin, pmax = np.sqrt(E - lam), np.sqrt(E + lam)
        f = lambda p: 1 / np.sqrt(max(lam ** 2 - (E - p * p) ** 2, 1e-300))
        Z = 2 * quad(f, pmin, pmax, points=[pmin, pmax], limit=400)[0]
        h = lambda p: (f(p) / Z) * np.log(f(p) / Z)
        Sd = -2 * quad(h, pmin, pmax, points=[pmin, pmax], limit=400)[0]
        print('   E/lam=%5.0f : S_diff = %.5f   closed form = %.5f'
              % (r, Sd, np.log(np.pi * lam / (2 * np.sqrt(E)))))

# DKS ground state cross-check: f = C^2 / (4 w W (w+W)^2) vs exact Gaussian
w, W, C = 1.0, 10.0, 0.5
f_pt = C ** 2 / (4 * w * W * (w + W) ** 2)
# exact: covariance of the ground state of H = (px^2+py^2+w^2x^2+W^2y^2+2Cxy)/2
K = np.array([[w * w, C], [C, W * W]])
ev, U = np.linalg.eigh(K)
Om = np.sqrt(ev)
X = U @ np.diag(1 / (2 * Om)) @ U.T
P = U @ np.diag(Om / 2) @ U.T
nu = np.sqrt(X[0, 0] * P[0, 0])
fex = (nu - 0.5) / (nu + 0.5)
print('DKS ground state: f perturbative = %.4e   exact (Gaussian) = %.4e' % (f_pt, fex))
