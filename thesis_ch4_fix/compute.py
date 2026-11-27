"""Corrected numerics for Chapters 3 and 4 (hbar = m = 1).  Writes results.json and figures."""
import json
import numpy as np
from numpy.polynomial.hermite import hermval
from math import factorial, sqrt, pi
from scipy.linalg import eigh

# ------------------------------------------------------------------ helpers
def hfun(m, y):
    c = np.zeros(m + 1); c[m] = 1
    return hermval(y, c) * np.exp(-y * y / 2) / sqrt(2 ** m * factorial(m) * sqrt(pi))

def aho_coeffs(n, g):
    """first-order coefficients for H = p^2/2 + x^2/2 + g x^4 (general n)."""
    c = {n: 1.0}
    c[n + 4] = -g * np.sqrt((n + 1) * (n + 2) * (n + 3) * (n + 4)) / 16
    c[n + 2] = -g * (2 * n + 3) * np.sqrt((n + 1) * (n + 2)) / 4
    if n >= 2:
        c[n - 2] = g * (2 * n - 1) * np.sqrt(n * (n - 1)) / 4
    if n >= 4:
        c[n - 4] = g * np.sqrt(n * (n - 1) * (n - 2) * (n - 3)) / 16
    N = 1 / np.sqrt(sum(v * v for v in c.values()))
    return {m: N * v for m, v in c.items()}

def psi_from(c, x):
    return sum(v * hfun(m, x) for m, v in c.items())

def phi_from(c, k, n):
    """momentum wave function divided by (-i)^n (real)."""
    return sum(v * (-1) ** ((m - n) // 2) * hfun(m, k) for m, v in c.items())

def shannon_1d(rho, d):
    r = rho[rho > 1e-300]
    return float(-np.sum(r * np.log(r)) * d)

def wigner(psi, x):
    """W(x_i, p_j) on the same grid for p, real psi, via sum over y = k dx."""
    d = x[1] - x[0]; N = len(x)
    p = x.copy()
    W = np.zeros((N, N))
    ks = np.arange(-N + 1, N)
    E = np.cos(2 * np.outer(ks * d, p))        # real part suffices (W real; sine terms cancel by symmetry)
    pad = np.concatenate([np.zeros(N), psi, np.zeros(N)])
    for i in range(N):
        f = pad[N + i + ks] * pad[N + i - ks]
        W[i] = (f @ E) * d / np.pi
    return p, W

def husimi(psi, x):
    """Q(x_i,p_j) = |<alpha|psi>|^2 / 2pi with symmetric coherent states (kernel e^{-(x-x')^2-(p-p')^2}/pi on W)."""
    d = x[1] - x[0]
    p = x.copy()
    G = np.exp(-(x[None, :] - x[:, None]) ** 2 / 2) * psi[None, :]
    A = G @ np.exp(-1j * np.outer(x, p)) * d
    return p, np.abs(A) ** 2 / (2 * np.pi * np.sqrt(np.pi))

# ------------------------------------------------------------------ Chapter 3
res = {}
x = np.linspace(-9, 9, 361); dx = x[1] - x[0]
ch3 = {}
lams = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10]
for n in (0, 1, 5):
    rows = []
    for lam in (lams if n < 5 else [1e-4]):
        g = lam / 4
        c = aho_coeffs(n, g)
        psi = psi_from(c, x); phi = phi_from(c, x, n)
        rx, rp = psi ** 2, phi ** 2
        p, W = wigner(psi, x)
        _, Q = husimi(psi, x)
        dA = dx * dx
        Wx, Wp = W.sum(1) * dx, W.sum(0) * dx
        Qx, Qp = Q.sum(1) * dx, Q.sum(0) * dx
        Sx, Sp = shannon_1d(rx, dx), shannon_1d(rp, dx)
        SxH, SpH = shannon_1d(Qx, dx), shannon_1d(Qp, dx)
        SH = shannon_1d(Q.ravel(), dA)
        Fx, Fp = 4 * np.sum(x ** 2 * rp) * dx, 4 * np.sum(x ** 2 * rx) * dx   # F_x = 4<p^2>, F_p = 4<x^2> (real states)
        R2W = -np.log(np.sum(W ** 2) * dA); R4W = -np.log(np.sum(W ** 4) * dA) / 3
        R2H = -np.log(np.sum(Q ** 2) * dA); R4H = -np.log(np.sum(Q ** 4) * dA) / 3
        m = (rx > 1e-300) & (Qx > 1e-300)
        KL = float(np.sum(rx[m] * np.log(rx[m] / Qx[m])) * dx)
        DCS = float(-np.log(np.sum(W * Q) * dA / np.sqrt(np.sum(W ** 2) * dA * np.sum(Q ** 2) * dA)))
        IH = SxH + SpH - SH
        rows.append(dict(lam=lam, g=g, normW=float(W.sum() * dA), normQ=float(Q.sum() * dA),
                         Wmax=float(W.max()), Wmin=float(W.min()), Qmax=float(Q.max()),
                         marg_err=float(np.abs(Wx - rx).max()), Sx=Sx, Sp=Sp, SxH=SxH, SpH=SpH, SH=SH,
                         Fx=float(Fx), Fp=float(Fp), R2W=float(R2W), R4W=float(R4W), R2H=float(R2H),
                         R4H=float(R4H), KL=KL, DCS=DCS, IH=IH))
        if lam in (1e-4,) or (n < 5 and lam == 0.0):
            np.savez('ps_n%d.npz' % n, x=x, W=W, Q=Q, rx=rx, Qx=Qx, Wx=Wx)
    ch3[n] = rows
    for r in rows:
        print('Ch3 n=%d lam=%.3f normW=%.6f normQ=%.6f W in [%.4f, %.4f] Qmax=%.4f marg_err=%.1e '
              'Sx=%.4f Sp=%.4f SxH=%.4f SpH=%.4f SH=%.4f Fx=%.4f Fp=%.4f R2W=%.5f R4W=%.4f R2H=%.4f '
              'R4H=%.4f KL=%.4f DCS=%.4f IH=%.4f' % (n, r['lam'], r['normW'], r['normQ'], r['Wmin'],
              r['Wmax'], r['Qmax'], r['marg_err'], r['Sx'], r['Sp'], r['SxH'], r['SpH'], r['SH'],
              r['Fx'], r['Fp'], r['R2W'], r['R4W'], r['R2H'], r['R4H'], r['KL'], r['DCS'], r['IH']))
res['ch3'] = ch3

# ------------------------------------------------------------------ exact diagonalisation
def sinc_dvr(V, L, N):
    x = np.linspace(-L, L, N); d = x[1] - x[0]
    i = np.arange(N)
    dif = i[:, None] - i[None, :]
    with np.errstate(divide='ignore'):
        T = np.where(dif == 0, np.pi ** 2 / 3, 2.0 / np.where(dif == 0, 1, dif) ** 2) * (-1.0) ** dif / (2 * d * d)
    e, v = eigh(T + np.diag(V(x)))
    return x, d, e, v

def box_dvr(V, a, b, N):
    """Colbert-Miller DVR on (a,b) with Dirichlet walls; N-1 interior points."""
    i = np.arange(1, N)
    x = a + (b - a) * i / N
    I, J = np.meshgrid(i, i, indexing='ij')
    with np.errstate(divide='ignore', invalid='ignore'):
        off = (1 / np.sin(np.pi * (I - J) / (2 * N)) ** 2 - 1 / np.sin(np.pi * (I + J) / (2 * N)) ** 2)
        diag = (2 * N ** 2 + 1) / 3 - 1 / np.sin(np.pi * i / N) ** 2
    T = np.where(I == J, np.diag(np.diag(np.ones((N - 1, N - 1))) * 0) if False else 0, off)
    T[np.arange(N - 1), np.arange(N - 1)] = diag
    T = T * (-1.0) ** (I - J) * (np.pi ** 2 / (b - a) ** 2) / 4   # hbar^2/2m * (pi/(b-a))^2 /2
    e, v = eigh(T + np.diag(V(x)))
    return x, (b - a) / N, e, v

def measures(x, d, vec, kmax=40.0, nk=8001, E=None, Vx=None):
    """S_x, S_k, F_x = 4<k^2>, F_k = 4<x^2>.  If E and Vx are given, <k^2> = 2(E - <V>) (exact on the DVR grid,
    free of the 1/k^4 momentum-tail truncation); otherwise <k^2> from the k grid."""
    psi = vec / np.sqrt(d)
    psi = psi * np.sign(psi[np.argmax(np.abs(psi))])
    rx = psi ** 2
    k = np.linspace(-kmax, kmax, nk); dk = k[1] - k[0]
    phi = (psi[None, :] * np.exp(-1j * np.outer(k, x))).sum(1) * d / np.sqrt(2 * np.pi)
    rk = np.abs(phi) ** 2
    Sx, Sk = shannon_1d(rx, d), shannon_1d(rk, dk)
    sx2 = float(np.sum(x ** 2 * rx) * d); sk2 = float(np.sum(k ** 2 * rk) * dk)
    if E is not None:
        sk2 = float(2 * (E - np.sum(Vx * vec ** 2)))
    return dict(Sx=Sx, Sk=Sk, Sxk=Sx + Sk, Fx=4 * sk2, Fk=4 * sx2, sx2=sx2, sk2=sk2, FxFk=16 * sx2 * sk2,
                normk=float(rk.sum() * dk))

def V_well(lam):
    def V(x):
        z = lam * x
        out = np.empty_like(z)
        s = np.abs(z) < 1e-4
        out[s] = 1 / 3 + z[s] ** 2 / 45
        out[~s] = (1 - z[~s] / np.tan(z[~s])) / z[~s] ** 2
        return out
    return V

quart = {}
for b in (0.01, 0.1, 1.0, 1.3, 1.5):
    Vq = lambda xx: xx ** 2 / 2 + b * xx ** 4
    x_, d_, e_, v_ = sinc_dvr(Vq, 7.0, 401)
    quart[b] = []
    for n in range(4):
        m = measures(x_, d_, v_[:, n], E=e_[n], Vx=Vq(x_)); m['E'] = float(e_[n]); quart[b].append(m)
        print('quartic exact b=%.2f n=%d E=%.6f Sx=%.5f Sk=%.5f Sx+Sk=%.5f Fx=%.4f Fk=%.4f sx2=%.4f sk2=%.4f FxFk=%.3f normk=%.6f'
              % (b, n, e_[n], m['Sx'], m['Sk'], m['Sxk'], m['Fx'], m['Fk'], m['sx2'], m['sk2'], m['FxFk'], m['normk']))
res['quartic_exact'] = {str(k): v for k, v in quart.items()}

# first-order PT values for the quartic at small b
qpt = {}
xg = np.linspace(-10, 10, 4001); dg = xg[1] - xg[0]
for b in (0.01, 0.05, 0.1):
    qpt[b] = []
    for n in range(4):
        c = aho_coeffs(n, b)
        psi = psi_from(c, xg); phi = phi_from(c, xg, n)
        sx2 = float(np.sum(xg ** 2 * psi ** 2) * dg); sk2 = float(np.sum(xg ** 2 * phi ** 2) * dg)
        r = dict(Sx=shannon_1d(psi ** 2, dg), Sk=shannon_1d(phi ** 2, dg), Fx=4 * sk2, Fk=4 * sx2,
                 E1=n + 0.5 + 0.75 * b * (2 * n * n + 2 * n + 1))
        r['Sxk'] = r['Sx'] + r['Sk']
        qpt[b].append(r)
        print('quartic PT b=%.2f n=%d Sx=%.5f Sk=%.5f Sx+Sk=%.5f Fx=%.4f Fk=%.4f E1=%.4f'
              % (b, n, r['Sx'], r['Sk'], r['Sxk'], r['Fx'], r['Fk'], r['E1']))
res['quartic_pt'] = {str(k): v for k, v in qpt.items()}

well = {}
for lam in (1.0, 1.3, 1.5, 2.0):
    L = np.pi / lam
    x_, d_, e_, v_ = box_dvr(V_well(lam), -L, L, 600)
    well[lam] = []
    for n in range(4):
        m = measures(x_, d_, v_[:, n], kmax=60.0, nk=12001, E=e_[n], Vx=V_well(lam)(x_)); m['E'] = float(e_[n]); well[lam].append(m)
        print('well exact lam=%.2f n=%d E=%.6f Sx=%.5f Sk=%.5f Sx+Sk=%.5f Fx=%.4f Fk=%.4f sx2=%.4f sk2=%.4f FxFk=%.3f normk=%.5f'
              % (lam, n, e_[n], m['Sx'], m['Sk'], m['Sxk'], m['Fx'], m['Fk'], m['sx2'], m['sk2'], m['FxFk'], m['normk']))
res['well_exact'] = {str(k): v for k, v in well.items()}
# exact ground state for V0 = lam^2: sin(lam x)/x, E = lam^2/2  (check with V0 = 1 at lam = 1)
x_, d_, e_, v_ = box_dvr(V_well(1.0), -np.pi, np.pi, 600)
ex = np.sin(x_) / np.where(np.abs(x_) < 1e-12, 1, x_); ex[np.abs(x_) < 1e-12] = 1
ex /= np.sqrt(np.sum(ex ** 2) * d_)
res['well_exact_check'] = dict(E0=float(e_[0]), overlap=float(abs(np.sum(ex * v_[:, 0] / np.sqrt(d_)) * d_)))
print('well lam=1: E0 = %.8f (exact 0.5); overlap with sin(x)/x = %.10f' % (e_[0], res['well_exact_check']['overlap']))

# many states (examiner's suggestion): S_x + S_k up to n = 100
many = {}
x_, d_, e_, v_ = sinc_dvr(lambda xx: xx ** 2 / 2 + xx ** 4, 9.0, 1201)
many['quartic_b1'] = [measures(x_, d_, v_[:, n], kmax=60, nk=6001)['Sxk'] for n in range(0, 101, 2)]
x_, d_, e_, v_ = box_dvr(V_well(1.0), -np.pi, np.pi, 1200)
many['well_lam1'] = [measures(x_, d_, v_[:, n], kmax=160, nk=16001)['Sxk'] for n in range(0, 101, 2)]
res['many'] = many
print('many: quartic', np.round(many['quartic_b1'][:6], 4), '...', round(many['quartic_b1'][-1], 4))
print('many: well   ', np.round(many['well_lam1'][:6], 4), '...', round(many['well_lam1'][-1], 4))
with open('results.json', 'w') as f:
    json.dump(res, f, indent=1)
