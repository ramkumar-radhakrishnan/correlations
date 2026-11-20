import numpy as np, texeval as te
from lattice import Lattice
lat = Lattice(L=5, N=3, seed=4); U, f, N = lat.U, lat.f, lat.N
rng = np.random.default_rng(1)
def colour(x, xp, z, w, wp, n):
    U = lat.U
    # Group II (7.3)-(7.8) colour factor with <rho^a(u) rho^b(v)> = delta^{ab} n(u,v)/(N^2-1)
    d = np.eye(lat.NA) / lat.NA
    P = n[xp, x] * d          # <rho^{d'}(x') rho^{d}(x)>
    # primes -> capitals: c'->C, d'->D, e'->E ; P[D,d] = <rho^{D}(x') rho^{d}(x)>
    t1 = np.einsum("CbD,cbd,aC,ac,Dd->", f, f, U[wp], U[w], P)
    t2 = np.einsum("CbD,acd,aC,cb,de,De->", f, f, U[wp], U[z], U[x], P)
    t3 = np.einsum("aCD,cbd,ac,Cb,DE,Ed->", f, f, U[w], U[z], U[xp], P)
    t4 = N * np.einsum("dE,de,Ee->", U[xp], U[x], P)
    return t1 - t2 - t3 + t4
n = rng.normal(size=(lat.S, lat.S)); n = n + n.T
mx = 0
for _ in range(200):
    x, xp, z, w, wp = rng.integers(0, lat.S, 5)
    a = colour(x, xp, z, w, wp, n); b = colour(xp, x, z, wp, w, n)
    mx = max(mx, abs(a - b) / (abs(a) + abs(b) + 1e-30))
print('Group II colour factor: max relative asymmetry under (x,w)<->(x\',w\'):', mx)
U0 = lat.U.copy(); lat.set_U_identity(); c1 = colour(1, 2, 3, 4, 5, n)
lat.set_U_constant(); c2 = colour(1, 2, 3, 4, 5, n); lat.U = U0
print('Group II colour factor at U=1:', abs(c1), '  at U=const:', abs(c2))

# U -> 1 and U -> constant for the parsed BFKL / JIMWLK blocks
B = te.blocks('../Final_terms_original.tex')
INT = ['x', "x'", 'z', 'w', "w'"]; PH = [('w', 1), ("w'", -1)]
P = rng.normal(size=(lat.S, lat.NA, lat.S, lat.NA))
for mode in ('identity', 'constant'):
    getattr(lat, 'set_U_' + mode)()
    ev = te.Evaluator(lat, P2=P)
    vals = [abs(ev.total(te.parse(te.body(B[i])), INT, PH)) for i in (0, 1, 3, 4, 5)]
    print(f'U={mode}: |(1.1)|,|(1.2)|,|(2.1)|,|(2.2)|,|(2.3)| =', ['%.1e' % v for v in vals])
