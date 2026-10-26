"""Small identities used in the note.
 1. G1 III of bfkl.tex: N1 = -p8 K_b O_III;  N2 + (N3a terms 2, 3) = p8 K_h O_III;  N3a term 1 = p8 K_f O_III.
 2. X'.F(w;x,y) = (1/pi) int_z F_z  (the z integral of Group I Row II 2a), numerically.
 3. The symmetric three-rho part of the y' -> z terms vanishes by itself (x <-> x' antisymmetry).
 4. The Fourier lemma behind the y' -> z terms: (1+zeta) int d^2r e^{-i(1+zeta)k.r} K(r).u g(z+r) -> -2 pi i (k.u)/k^2 g(z).
 5. Group I Row V: the third term of 2b finite is the pointwise form of -(2b UV):  (y-x).(y-z)/((y-x)^2 (x-z)^2) = K(z,z,x) - K(y,z,x)."""
import json
import numpy as np
import env  # noqa: F401
from fp import G, PI, sq
import split as SP
import ft
from ft import PI3
import three_rho_test as TT
import three_rho_engine as T3
from core import symmetrize
import rows as R

out = {}
K = lambda v: v / sq(v)
p8 = 1j * G**4 / (8 * PI**5)

# 1
rel = []
for s in range(1, 6):
    c = SP.cfg(s)
    c.kv = np.zeros(2)
    x, y, z, w, xp, wp = (c.pt[n] for n in ('x', 'y', 'z', 'w', 'xp', 'wp'))
    xpwp, yz, xw, xz, zw = xp - wp, y - z, x - w, x - z, z - w
    A, B = sq(xz), sq(xw)
    D = lambda a, b: a @ b
    s6 = (D(xpwp, zw) * D(yz, xw) / (sq(xpwp) * sq(yz) * sq(zw) * (B - A))
          - D(xpwp, xw) * D(yz, zw) * A / (sq(xpwp) * sq(yz) * sq(zw) * (B - A) * B)
          + D(xpwp, yz) * D(zw, xw) * A / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2)
          - D(xpwp, zw) * D(yz, xz) * B / (sq(xpwp) * sq(yz) * sq(zw) * (B - A) * A)
          + D(xpwp, xz) * D(yz, zw) / (sq(xpwp) * sq(yz) * sq(zw) * (B - A))
          - D(xpwp, yz) * D(zw, xz) * B / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2))
    n2 = -(-D(xpwp, xw) * D(yz, zw) / (sq(xpwp) * sq(yz) * sq(zw) * B) + s6)          # in units of p8 L_v
    n3_23 = -(-D(xpwp, yz) * D(zw, xw) * A / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2)
              + D(xpwp, yz) * D(zw, xz) * B / (sq(xpwp) * sq(yz) * sq(zw) * (B - A)**2))
    Kh = D(K(xpwp), K(zw)) * D(K(yz), K(xz))
    rel.append(abs(n2 + n3_23 - Kh) / abs(Kh))
out['G1 III: N2 + N3a(2,3) = p8 K_h'] = max(rel)

# 2
P = {k: np.random.default_rng(7).normal(size=2) for k in ()}
rng = np.random.default_rng(11)
res2 = []
for trial in range(2):
    P = {k: rng.normal(size=2) for k in ('xp', 'wp', 'x', 'y', 'w')}
    import cw_test as CW

    def fz(P_, z):
        Xp = K(P_['xp'] - P_['wp'])
        x, y, w = P_['x'], P_['y'], P_['w']
        return ((Xp @ K(x - w)) * (K(y - z) @ (K(x - z) - K(w - z))) - (Xp @ K(y - w)) * (K(x - z) @ (K(y - z) - K(w - z))))
    I = CW.integrate(fz, P, ['w', 'x', 'y'])
    res2.append(abs(I / np.pi - CW.F_dot(P)) / abs(CW.F_dot(P)))
out["X'.F = (1/pi) int_z F_z"] = max(res2)

# 3
S = {k: [s for t in f() for s in symmetrize(t, keep=(3,))] for k, f in R.THREE_YZ.items()}
r3 = []
for seed in (1, 2, 3):
    cfg = T3.Cfg3(seed)
    tot, sc = 0.0, 0.0
    for ts in S.values():
        v, s_ = T3.value3(ts, cfg)
        tot += v
        sc += s_
    r3.append(abs(tot) / sc)
out["y'->z terms: symmetric three-rho part / sum|terms|"] = max(r3)

# 4  exact for a Gaussian g: int d^2r e^{-iq.r} (r.u/r^2) e^{-r^2/(2 s^2)} = -2 pi i (q.u)/q^2 (1 - e^{-q^2 s^2/2})
#    numerically for g = e^{-r^2/2} (1 + r.c) with polar quadrature
def ft_num(q, u, cvec, nr=4000, nt=512):
    r = np.linspace(0, 12, nr + 1)[1:] - 6 / nr
    th = (np.arange(nt) + 0.5) * 2 * np.pi / nt
    R_, T_ = np.meshgrid(r, th, indexing='ij')
    X, Y = R_ * np.cos(T_), R_ * np.sin(T_)
    g = np.exp(-(X**2 + Y**2) / 2) * (1 + X * cvec[0] + Y * cvec[1])
    f = np.exp(-1j * (q[0] * X + q[1] * Y)) * (X * u[0] + Y * u[1]) / (X**2 + Y**2) * g
    return np.sum(f * R_) * (12 / nr) * (2 * np.pi / nt)


k = np.array([0.7, -0.4])
u = np.array([0.3, 1.1])
cvec = np.array([0.5, -0.2])
lem = []
for zeta in (5, 10, 20, 40):
    q = (1 + zeta) * k
    val = ft_num(q, u, cvec) * (1 + zeta)
    lim = -2j * np.pi * (k @ u) / (k @ k)
    lem.append(dict(zeta=zeta, ratio=[val.real / lim.real if lim.real else 0, val.imag / lim.imag]))
out['Fourier lemma (1+zeta) FT / limit'] = lem

# 5
rel5 = []
for s in range(5):
    x, y, z = (np.random.default_rng(s).normal(size=2) for _ in range(3))
    rng5 = np.random.default_rng(100 + s)
    x, y, z = rng5.normal(size=2), rng5.normal(size=2), rng5.normal(size=2)
    lhs = ((y - x) @ (y - z)) / (sq(y - x) * sq(x - z))
    rhs = K(z - x) @ K(z - x) - K(y - x) @ K(z - x)
    rel5.append(abs(lhs - rhs) / abs(lhs))
out['2b finite term 3 = K(z,z,x) - K(y,z,x)'] = max(rel5)
for k_, v in out.items():
    print(k_, v)
json.dump(out, open('checks.json', 'w'), indent=1, default=float)
