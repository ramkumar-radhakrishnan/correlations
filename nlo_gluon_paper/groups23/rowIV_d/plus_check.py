"""Small-xi behaviour of the d_perp term of Group II Row IV (third line of Eq. (fin) of RowIV_FF_userform.tex):
    T(xi0, xi1) = int_{xi0}^{xi1} dxi (1-xi)/xi F(k/xi),
    F(q) = int_{y,y'} e^{-iq.(y'-y)} f^{c,i}(y) f^{c,j}(y') [delta^{ij} - 2 Dh^i Dh^j],   Dh = (y'-y)/|y'-y|,
    f^{c,i}(y) = int_x X^i(x-y) [U^{ce}(y) - U^{ce}(x)] rho^e(x),   X(r) = r/r^2.
Claim: F(q) = (8 pi/q^2) C^{ij}(0) (qh^i qh^j - delta^{ij}/2) + O(1/q^3), C^{ij}(D) = int_y f^{c,i}(y) f^{c,j}(y+D),
so the integrand is O(xi) at xi -> 0: no log(vee/k+), the lower limit can be set to 0 (error O(xi0^2)).
For contrast the pole/log term of the same row, weight dxi (1-xi)/xi^3 and F_E(q) = int e^{-iq.D} [-log D^2] C^{ii}(D),
behaves as 4 pi C^{ii}(0)/q^2 and gives log(1/xi0) = log(vee/k+) (the BFKL piece).
Configuration: SU(2) in the adjoint (SO(3) rotations), smooth U(y) -> 1 at large y, two Gaussian colour charges.
Results in plus_check.json."""
import json
import os

import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import jv

HERE = os.path.dirname(os.path.abspath(__file__))
R = {}
N, L = 3072, 40.0
h = L / N
n = np.arange(N)
c1 = np.where(n < N // 2, n, n - N) * h                      # wrapped coordinates (also used for D)
Y1, Y2 = np.meshgrid(c1, c1, indexing='ij')
p1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
P1, P2 = np.meshgrid(p1, p1, indexing='ij')
PP = P1**2 + P2**2
PP[0, 0] = 1.0


def conv_X(g):
    """int_x X(x-y) g(x) for both components, by FFT: FT[-r/r^2] = 2 pi i p/p^2."""
    G = np.fft.fft2(g)
    out = []
    for P in (P1, P2):
        K = 2j * np.pi * P / PP
        K[0, 0] = 0.0
        out.append(np.fft.ifft2(G * K).real)
    return out


# ---------------------------------------------------------------- 0. calibration: (delta-2DhDh) M e^{-D^2} on the grid
D2 = Y1**2 + Y2**2
with np.errstate(invalid='ignore', divide='ignore'):
    Dh1, Dh2 = Y1 / np.sqrt(D2), Y2 / np.sqrt(D2)
Dh1[0, 0] = Dh2[0, 0] = 0.0
Tm = {(0, 0): 1 - 2 * Dh1 * Dh1, (1, 1): 1 - 2 * Dh2 * Dh2, (0, 1): -2 * Dh1 * Dh2, (1, 0): -2 * Dh1 * Dh2}
Tm[(0, 0)][0, 0] = Tm[(1, 1)][0, 0] = Tm[(0, 1)][0, 0] = Tm[(1, 0)][0, 0] = 0.0   # angular average at D = 0


def ft_line(A, q, ang):
    """int d^2D e^{-i q.D} A(D) for q = q (cos ang, sin ang), direct sum."""
    qx, qy = q * np.cos(ang), q * np.sin(ang)
    return np.exp(-1j * qx * c1) @ A @ np.exp(-1j * qy * c1) * h * h


M = np.array([[0.7, 0.3], [-0.2, -0.4]])
Gcal = sum(M[i, j] * Tm[(i, j)] for i in range(2) for j in range(2)) * np.exp(-D2)
cal = []
for q in (1.0, 5.0, 20.0, 40.0):
    for ang in (0.0, 0.6):
        qh = np.array([np.cos(ang), np.sin(ang)])
        hank = quad(lambda r: r * jv(2, q * r) * np.exp(-r * r), 0, 12, limit=400)[0]
        exact = 4 * np.pi * np.einsum('ij,ij', M, np.outer(qh, qh) - np.eye(2) / 2) * hank
        num = ft_line(Gcal, q, ang)
        cal.append(dict(q=q, ang=ang, exact=exact, grid=num.real, imag=num.imag, q2_exact=q * q * exact,
                        q2_asym=8 * np.pi * np.einsum('ij,ij', M, np.outer(qh, qh) - np.eye(2) / 2)))
R['calibration'] = cal
R['calibration_max_abs_err_times_q2'] = max(abs(c['grid'] - c['exact']) * c['q']**2 for c in cal)

# ---------------------------------------------------------------- 1. configuration
rng = np.random.default_rng(7)


def rot(th):
    """exp(theta^a L^a), (L^a)_{bc} = -eps_{abc}: adjoint of SU(2), Rodrigues; th shape (3, N, N)."""
    t = np.sqrt((th**2).sum(0)) + 1e-300
    k = th / t
    Kx = np.zeros((3, 3) + th.shape[1:])
    Kx[0, 1], Kx[0, 2], Kx[1, 2] = -k[2], k[1], -k[0]
    Kx[1, 0], Kx[2, 0], Kx[2, 1] = k[2], -k[1], k[0]
    K2 = np.einsum('ab...,bc...->ac...', Kx, Kx)
    I = np.eye(3)[:, :, None, None]
    return I + np.sin(t) * Kx + (1 - np.cos(t)) * K2


theta = np.zeros((3, N, N))
for _ in range(6):
    cen = rng.uniform(-2.5, 2.5, 2)
    amp = rng.normal(0, 1.4, 3)
    w = rng.uniform(0.6, 1.4)
    theta += amp[:, None, None] * np.exp(-((Y1 - cen[0])**2 + (Y2 - cen[1])**2) / w**2)
U = rot(theta)                                                # U^{ce}(y), -> 1 at large y
# two Gaussian charges; the colour of the second makes Q^c = int_x [1 - U(x)]^{ce} rho^e(x) vanish, so that
# f ~ 1/y^2 at large y (otherwise f ~ Q/|y| and C(D) ~ log|D| does not decay in the box)
gs = [np.exp(-((Y1 - cx)**2 + (Y2 - cy)**2) / (2 * s * s)) / (2 * np.pi * s * s) for cx, cy, s in ((-0.9, 0.4, 0.25), (1.1, -0.6, 0.35))]
A = [np.einsum('ceab,ab->ce', np.eye(3)[:, :, None, None] - U, g) * h * h for g in gs]
col1 = rng.normal(0, 1, 3)
col2 = -np.linalg.solve(A[1], A[0] @ col1)
R['colours'] = dict(col1=col1.tolist(), col2=col2.tolist(), Q=(A[0] @ col1 + A[1] @ col2).tolist())
rho = col1[:, None, None] * gs[0] + col2[:, None, None] * gs[1]

E = [conv_X(rho[e]) for e in range(3)]                       # int_x X^i(x-y) rho^e(x)
Urho = np.einsum('ce...,e...->c...', U, rho)                 # U^{ce}(x) rho^e(x)
Et = [conv_X(Urho[c]) for c in range(3)]
f = np.zeros((3, 2, N, N))
for c in range(3):
    for i in range(2):
        f[c, i] = sum(U[c, e] * E[e][i] for e in range(3)) - Et[c][i]
del U, Urho, E, Et, theta

# C^{ij}(D) = sum_c int_y f^{c,i}(y) f^{c,j}(y+D)
Fh = np.fft.fft2(f, axes=(-2, -1))
C = np.zeros((2, 2, N, N))
for i in range(2):
    for j in range(2):
        C[i, j] = np.fft.ifft2((np.conj(Fh[:, i]) * Fh[:, j]).sum(0)).real * h * h
del Fh, f
C /= np.trace(C[:, :, 0, 0])                                  # normalization: tr C(0) = 1 (everything is bilinear)
C0 = C[:, :, 0, 0]
R['C0'] = C0.tolist()

# smooth window in D: it changes F(q) only at small q (q <~ 1/Rw); the large-q tail comes from D -> 0
Rw = 7.0
win = np.exp(-D2 / (2 * Rw * Rw))
R['C_at_box_edge_over_C0'] = float(np.abs(C[:, :, N // 2, 0]).max() / np.abs(C0).max())
CT = sum(C[i, j] * Tm[(i, j)] for i in range(2) for j in range(2)) * win
# -log D^2, with its cell average at D = 0: the mean of log(u^2+v^2) over [-1/2,1/2]^2 is log(1/2) - 3 + pi/2
a_sq = np.log(0.5) - 3 + np.pi / 2
R['a_sq_check'] = [a_sq, 4 * dblquad(lambda r, t: r * np.log(r * r), 0, np.pi / 4, 0, lambda t: 0.5 / np.cos(t))[0] * 2]
mlog = -np.log(np.where(D2 > 0, D2, 1.0))
mlog[0, 0] = -(np.log(h * h) + a_sq)
CE = (C[0, 0] + C[1, 1]) * mlog * win
del C

# ---------------------------------------------------------------- 2. q^2 F(q) vs the asymptotic forms
asym = []
for ang in (0.0, 0.9, 2.0):
    qh = np.array([np.cos(ang), np.sin(ang)])
    predT = 8 * np.pi * np.einsum('ij,ij', C0, np.outer(qh, qh) - np.eye(2) / 2)
    predE = 4 * np.pi * np.trace(C0)
    for q in (2.0, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0):
        FT, FE = ft_line(CT, q, ang), ft_line(CE, q, ang)
        asym.append(dict(ang=ang, q=q, q2F_T=[(q * q * FT).real, (q * q * FT).imag], pred_T=predT,
                         q2F_E=[(q * q * FE).real, (q * q * FE).imag], pred_E=predE))
R['asymptotics'] = asym

# ---------------------------------------------------------------- 3. the xi integrals, |k| = kk along ang = 0
kk = 1.5
gx, gw = np.polynomial.legendre.leggauss(48)


def xi_int(A, weight, a, b):
    """int_a^b dxi weight(xi) F(kk/xi), Gauss-Legendre in log xi, panels of one e-fold."""
    tot = 0.0
    edges = np.unique(np.concatenate([np.log(np.geomspace(a, b, max(2, int(np.ceil(np.log(b / a))) + 1)))]))
    for lo, hi in zip(edges[:-1], edges[1:]):
        t = 0.5 * (hi - lo) * gx + 0.5 * (hi + lo)
        for tt, ww in zip(t, gw):
            x = np.exp(tt)
            tot += 0.5 * (hi - lo) * ww * x * weight(x) * ft_line(A, kk / x, 0.0)
    return tot


wT = lambda x: (1 - x) / x
wE = lambda x: (1 - x) / x**3
xi0s = [0.2, 0.1, 0.05, 0.025]
qh = np.array([1.0, 0.0])
cT = 8 * np.pi * np.einsum('ij,ij', C0, np.outer(qh, qh) - np.eye(2) / 2) / kk**2   # F_T ~ cT xi^2
cE = 4 * np.pi * np.trace(C0) / kk**2                                                # F_E ~ cE xi^2
segT, segE = [], []
for a, b in zip(xi0s[1:], xi0s[:-1]):
    segT.append(xi_int(CT, wT, a, b))
    segE.append(xi_int(CE, wE, a, b))
base_T = xi_int(CT, wT, xi0s[0], 1.0)
base_E = xi_int(CE, wE, xi0s[0], 1.0)
rows = []
accT, accE = base_T, base_E
rows.append(dict(xi0=xi0s[0], T=[accT.real, accT.imag], TE=[accE.real, accE.imag]))
for i, x0 in enumerate(xi0s[1:]):
    accT += segT[i]
    accE += segE[i]
    rows.append(dict(xi0=x0, T=[accT.real, accT.imag], TE=[accE.real, accE.imag],
                     dT=[segT[i].real, segT[i].imag],
                     dT_pred=cT * ((xi0s[i]**2 - x0**2) / 2 - (xi0s[i]**3 - x0**3) / 3),
                     dTE=[segE[i].real, segE[i].imag], dTE_pred=cE * (np.log(xi0s[i] / x0) - (xi0s[i] - x0))))
R['xi_integrals'] = dict(k=kk, cT=cT, cE=cE, rows=rows,
                         T_extrapolated_to_0=(accT + cT * (xi0s[-1]**2 / 2 - xi0s[-1]**3 / 3)).real)

# ---------------------------------------------------------------- 4. the xi -> 1 end: upper limit k+/(k++Lambda) vs 1
F1 = ft_line(CT, kk, 0.0)
up = []
for lk in (0.3, 0.1, 0.03, 0.01):
    x1 = 1 / (1 + lk)
    t = 0.5 * (1 - x1) * gx + 0.5 * (1 + x1)
    strip = sum(0.5 * (1 - x1) * ww * (1 - x) / x * ft_line(CT, kk / x, 0.0) for x, ww in zip(t, gw))
    up.append(dict(Lambda_over_kp=lk, strip=[strip.real, strip.imag], pred=(F1 * (1 - x1)**2 / 2).real))
R['upper_end'] = dict(F_at_xi1=[F1.real, F1.imag], rows=up)

json.dump(R, open(os.path.join(HERE, 'plus_check.json'), 'w'), indent=1, default=float)
print('calibration max |err| q^2:', R['calibration_max_abs_err_times_q2'])
print('C0 =', C0)
for r in asym:
    print('ang %.1f q %5.1f  q2F_T %+.5f%+.5fi pred %+.5f | q2F_E %+.5f%+.5fi pred %+.5f' % (
        r['ang'], r['q'], *r['q2F_T'], r['pred_T'], *r['q2F_E'], r['pred_E']))
print('cT', cT, 'cE', cE)
for r in rows:
    print(r)
print('T extrapolated', R['xi_integrals']['T_extrapolated_to_0'])
for r in up:
    print(r)
print('box edge', R['C_at_box_edge_over_C0'], 'a_sq', R['a_sq_check'], 'Q', R['colours']['Q'])
