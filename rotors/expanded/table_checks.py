"""Every analytic claim in the paper vs. exact numerics; prints LaTeX rows."""
import numpy as np
from scipy.special import mathieu_a, mathieu_b, jv, ellipk
from scipy.integrate import quad
from rotorlib import *

rows = []
def row(what, ana, num):
    rows.append((what, ana, num))

# perturbative amplitudes, n = 7, lam = 0.02
lam, (m1, m2) = 0.02, (9, 2); n = 7; a = lam / 2
psi, E, st, idx = exact_state(m1, m2, 40, lam)
psi = psi / psi[idx[(m1, m2)]]
g = lambda d: psi[idx[(m1 + d, m2 - d)]]
row(r'$c_{+}$, Eq.~\eqref{eq:first}', -a / (n + 1), g(1))
row(r'$c_{-}$, Eq.~\eqref{eq:first}', a / (n - 1), g(-1))
row(r'$c_{++}$, Eq.~\eqref{eq:second}', a * a / (2 * (n + 1) * (n + 2)), g(2))
row(r'$c_{--}$, Eq.~\eqref{eq:second}', a * a / (2 * (n - 1) * (n - 2)), g(-2))
row(r'$E^{(2)}$, Eq.~\eqref{eq:E2}', lam ** 2 / (2 * (n * n - 1)), E - 0.5 * (m1 ** 2 + m2 ** 2))
# normalisation coefficient (branch-resolved), lam = 0.4
lam = 0.4; a = lam / 2
psi, E, st, idx = exact_state(m1, m2, 40, lam)
pos = np.array([psi[i] ** 2 for i, (u, v) in enumerate(st) if u > (m1 + m2) / 2])
c0sq = psi[idx[(m1, m2)]] ** 2 / pos.sum()
row(r'$p_{0}=1-\epsp-\epsm$ at $\lambda=0.4$', 1 - a * a / 64 - a * a / 36, c0sq)
# weak coupling S at n=0
lam = 0.02; a = lam / 2
psi, E, st, idx = exact_state(4, 4, 30, lam)
Sx = shannon(rho1_diag(psi, st, 30)); ep = a * a
row(r'$S$ at $n=0$, $\lambda=0.02$, Eq.~\eqref{eq:Sweak}',
    -(1 - 2 * ep) * np.log(1 - 2 * ep) - 2 * ep * np.log(ep), Sx)
# large-n asymptotic, n = 30, lam = 0.05
n, lam = 30, 0.05
psi, E, st, idx = exact_state(30, 0, 80, lam)
row(r'$S$ at $n=30$, $\lambda=0.05$, Eq.~\eqref{eq:Sasym}',
    (lam / n) ** 2 * (np.log(2 * n / lam) + 0.5), shannon(rho1_diag(psi, st, 80)))
# Mathieu
lam = 7.3; q = 2 * lam
H, st, idx = sector_H(0, 80, lam); w0 = np.sort(np.linalg.eigvalsh(H))
H, st, idx = sector_H(1, 80, lam); w1 = np.sort(np.linalg.eigvalsh(H)) - 0.25
row(r'$M=0$ lowest level vs $a_{0}(q)/4$', mathieu_a(0, q) / 4, w0[0])
row(r'$M=0$ third level vs $a_{2}(q)/4$', mathieu_a(2, q) / 4, w0[2])
row(r'$M=1$ lowest level vs $b_{1}(q)/4$', mathieu_b(1, q) / 4, w1[0])
row(r'$M=1$ second level vs $a_{1}(q)/4$', mathieu_a(1, q) / 4, w1[1])
# elliptic band integrals
for E_, l_ in ((300., 100.), (60., 100.)):
    pmin, pmax = np.sqrt(max(E_ - l_, 0)), np.sqrt(E_ + l_)
    num = quad(lambda p: 1 / np.sqrt(max(l_ ** 2 - (E_ - p * p) ** 2, 1e-300)),
               pmin, pmax, points=[pmin, pmax], limit=400)[0]
    ana = (np.sqrt(1 / (E_ + l_)) * ellipk(2 * l_ / (E_ + l_)) if E_ > l_
           else np.sqrt(1 / (2 * l_)) * ellipk((E_ + l_) / (2 * l_)))
    row(r'$B$, $\Erel=%d$, $\lambda=%d$, Eqs.~\eqref{eq:Brot},\eqref{eq:Blib}' % (E_, l_), ana, num)
# Bessel
n, lam = 400, 2000.
N = n // 2 + 80
H, st, idx = sector_H(0, N, lam); w, v = np.linalg.eigh(H)
k = int(np.argmin(np.abs(w - (n // 2) ** 2)))
p = np.zeros(2 * N + 1)
for i, (u, vv) in enumerate(st):
    p[u + N] += v[i, k] ** 2
row(r'$p_{j=0}$ at $n=400$, $\gamma=5$ vs $J_{0}^{2}(5)$', jv(0, 5.) ** 2, p[N + n // 2])
row(r'$p_{j=3}$ at $n=400$, $\gamma=5$ vs $J_{3}^{2}(5)$', jv(3, 5.) ** 2, p[N + n // 2 + 3])
# arcsine entropy
gm = 1e4
row(r'arcsine entropy $\ln(\pi\gamma/2)$, $\gamma=10^{4}$', np.log(np.pi * gm / 2),
    -quad(lambda x: (1 / (np.pi * np.sqrt(gm * gm - x * x))) *
          np.log(1 / (np.pi * np.sqrt(gm * gm - x * x))), -gm, gm, limit=400)[0])
# deficit
jm = int(1e6 + 4000); j = np.arange(-jm, jm + 1)
pj = jv(j, 1e6) ** 2
e = np.clip(np.arange(-jm - 0.5, jm + 1.5), -1e6, 1e6)
wc = (np.arcsin(e[1:] / 1e6) - np.arcsin(e[:-1] / 1e6)) / np.pi
row(r'$\Scl^{\rm bin}-S_{q}$ at $\gamma=10^{6}$ vs $1-\ln2$', 1 - np.log(2), shannon(wc) - shannon(pj))
# deep libration
E, D, P, N = pendulum_parity(10240.)
row(r'$S_{0}$ at $\lambda I/\hbar^{2}=10240$, Eq.~\eqref{eq:S0}',
    0.25 * np.log(10240) + 0.5 * np.log(np.pi * np.e / np.sqrt(2)) - 1 / (8 * np.sqrt(2 * 10240)),
    shannon(D[0]))
E, D, P, N = pendulum_parity(800.)
row(r'$N_{\rm lib}$ at $\lambda I/\hbar^{2}=800$ vs $(4/\pi)\sqrt q$', 4 / np.pi * np.sqrt(1600), (E < 800).sum())

def fmt(x):
    if abs(x) >= 1e-3 and abs(x) < 1e4:
        return '%.7f' % x if abs(x) < 10 else '%.5f' % x
    e = int(np.floor(np.log10(abs(x)))); return r'$%.6f\times10^{%d}$' % (x / 10 ** e, e)
with open('table_checks.tex', 'w') as f:
    for what, ana, num in rows:
        f.write('%s & %s & %s \\\\\n' % (what, fmt(ana), fmt(num)))
for what, ana, num in rows:
    print('%-70s %16.10g %16.10g' % (what[:70], ana, num))
