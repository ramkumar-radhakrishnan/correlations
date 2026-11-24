# c2_colour.py -- colour checks for Explain_Row_I_II_3rho.tex (Group I Row II, symmetrised three-rho term)
# Run: python3 -I c2_colour.py -> c2_colour.log
# C1  d^{abc}d^{abc} = (N^2-4)(N^2-1)/N,  f^{hbc} d^{abc} = 0,  Tr T^a T^b = N delta^{ab}
# C2  f^{gde} U^{da} U^{eb} = U^{gh} f^{hab} (adjoint Wilson line),  U real orthogonal
# C3  trace rule  f^{ijk} d^{lmn} O^{il} B^{jm} C^{kn} = i Tr{T^h O^T C D^h B^T O}
# C4  printed coefficient x bracket  ==  2 C^{abc} d^{abc} x normalisation   (direct contraction), SU(3), SU(4), SU(5)
# C5  properties: U(w) drops out; real; antisymmetric under x <-> y; zero at x = y, at x' = w', for constant U; V -> V^* flips sign; eps^3 for a weak target; SU(2)
# C6  operator identity  rho^a(1){rho^b(2),rho^c(3)} = 2{...}_sym + contact + one-rho  (explicit charges on sites)
# C7  exact projectile states: single quark (o = d.d/(4N)), and the average formula
import numpy as np, itertools
from scipy.linalg import expm
out = []
def pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

def gens(N):
    t = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex); m[i, j] = m[j, i] = 0.5; t.append(m)
            m = np.zeros((N, N), complex); m[i, j] = -0.5j; m[j, i] = 0.5j; t.append(m)
    for l in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(l): m[i, i] = 1
        m[l, l] = -l; t.append(m / np.sqrt(2 * l * (l + 1)))
    return np.array(t)
def fd(t):
    comm = np.einsum('aij,bjk->abik', t, t) - np.einsum('bij,ajk->abik', t, t)
    anti = np.einsum('aij,bjk->abik', t, t) + np.einsum('bij,ajk->abik', t, t)
    f = np.real(-2j * np.einsum('abik,cki->abc', comm, t)); d = np.real(2 * np.einsum('abik,cki->abc', anti, t))
    return f, d
def randV(N, rng, scale=1.0, t=None):
    c = rng.normal(size=N * N - 1) * scale
    return expm(1j * np.einsum('a,aij->ij', c, t)), c
def adj(V, t):
    return np.real(2 * np.einsum('aij,jk,bkl,li->ab', t, V, t, V.conj().T))

rng = np.random.default_rng(2026)
for N in (3, 4, 5):
    t = gens(N); f, d = fd(t); n = N * N - 1
    T = -1j * f; D = d.astype(complex)                       # (T^h)_{bc} = -i f^{hbc}, (D^h)_{bc} = d^{hbc}
    norm = N / ((N * N - 4) * (N * N - 1))
    pr('--- SU(%d)' % N)
    pr('C1  d.d - (N^2-4)(N^2-1)/N = %.1e ; max|f^{hbc} d^{abc}| = %.1e ; max|Tr T^a T^b - N delta| = %.1e'
       % (abs(np.einsum('abc,abc', d, d) - (N * N - 4) * (N * N - 1) / N), np.abs(np.einsum('hbc,abc->ha', f, d)).max(),
          np.abs(np.einsum('aij,bji->ab', T, T) - N * np.eye(n)).max()))
    devs = {k: 0.0 for k in ['C2', 'C3', 'C4', 'w', 'im', 'anti', 'xy', 'xw', 'const', 'cc']}
    vals = []
    for conf in range(5):
        V = {k: randV(N, rng, 1.0, t)[0] for k in ['x', 'y', 'xp', 'w', 'wp', 'w2']}
        Ux, Uy, Uxp, Uw, Uwp, Uw2 = [adj(V[k], t) for k in ['x', 'y', 'xp', 'w', 'wp', 'w2']]
        devs['C2'] = max(devs['C2'], np.abs(np.einsum('gde,da,eb->gab', f, Ux, Ux) - np.einsum('gh,hab->gab', Ux, f)).max(),
                         np.abs(Ux @ Ux.T - np.eye(n)).max())
        # C3 with arbitrary B, C
        B = rng.normal(size=(n, n)); C = rng.normal(size=(n, n))
        lhs = np.einsum('ijk,lmn,il,jm,kn', f, d, Ux, B, C)
        rhs = 1j * np.einsum('hab,bc,hcd,da', T, Ux.T @ C, D, B.T @ Ux)
        devs['C3'] = max(devs['C3'], abs(lhs - rhs) / abs(lhs))
        def Ctens(Uxp, Uwp, Uw, Ux, Uy):                     # C^{abc} of the row (Sec. 6 of the source)
            return np.einsum('ga,gbc->abc', Uxp - Uwp, np.einsum('hbc,gh->gbc', f, Uw) - np.einsum('gde,db,ec->gbc', f, Ux, Uy))
        def tr(Ua, Ub, Uc):                                   # Tr{T^h U^T(a) U(b) D^h U^T(c) U(a)}
            return np.einsum('hab,bc,hcd,da', T, Ua.T @ Ub, D, Uc.T @ Ua)
        def printed(Uxp, Uwp, Ux, Uy):
            return norm * 2j * (tr(Ux, Uy, Uxp) - tr(Ux, Uy, Uwp))
        direct = 2 * norm * np.einsum('abc,abc', Ctens(Uxp, Uwp, Uw, Ux, Uy), d)
        pv = printed(Uxp, Uwp, Ux, Uy); vals.append(direct)
        devs['C4'] = max(devs['C4'], abs(pv - direct) / abs(direct))
        devs['w'] = max(devs['w'], abs(np.einsum('abc,abc', Ctens(Uxp, Uwp, Uw2, Ux, Uy) - Ctens(Uxp, Uwp, Uw, Ux, Uy), d)))
        devs['im'] = max(devs['im'], abs(pv.imag))
        devs['anti'] = max(devs['anti'], abs(printed(Uxp, Uwp, Uy, Ux) + pv) / abs(pv))
        devs['xy'] = max(devs['xy'], abs(printed(Uxp, Uwp, Ux, Ux)))
        devs['xw'] = max(devs['xw'], abs(printed(Uwp, Uwp, Ux, Uy)))
        devs['const'] = max(devs['const'], abs(printed(Ux, Ux, Ux, Ux)), abs(np.einsum('abc,abc', Ctens(Ux, Ux, Ux, Ux, Ux), d)))
        Us = [adj(V[k].conj(), t) for k in ['xp', 'wp', 'x', 'y']]
        devs['cc'] = max(devs['cc'], abs(printed(*Us) + pv) / abs(pv))
    pr('C2  f U U = U f and U U^T = 1: %.1e' % devs['C2'])
    pr('C3  trace rule f d O B C = i Tr{T^h O^T C D^h B^T O} (random B, C; O adjoint Wilson line): rel dev %.1e' % devs['C3'])
    pr('C4  printed [N/((N^2-4)(N^2-1))] (2i)[Tr{T U^T(x)U(y) D U^T(x\')U(x)} - Tr{.. U^T(w\') ..}]  vs  2 C^{abc}d^{abc} N/((N^2-4)(N^2-1)):'
       ' max rel dev %.1e ; values (5 configurations): %s' % (devs['C4'], ' '.join('%+.5f' % v for v in vals)))
    pr('C5  U(w) -> other U(w): change of C.d %.1e ; |Im printed| %.1e ; x<->y antisymmetry %.1e ; at x=y %.1e ; at x\'=w\' %.1e ;'
       ' constant U: %.1e ; V -> V^*: printed -> -printed, rel dev %.1e'
       % (devs['w'], devs['im'], devs['anti'], devs['xy'], devs['xw'], devs['const'], devs['cc']))
    # weak target
    cs = {k: rng.normal(size=n) for k in ['x', 'y', 'xp', 'wp']}
    def val(eps):
        Us = {k: adj(expm(1j * eps * np.einsum('a,aij->ij', cs[k], t)), t) for k in cs}
        return (norm * 2j * (np.einsum('hab,bc,hcd,da', T, Us['x'].T @ Us['y'], D, Us['xp'].T @ Us['x'])
                             - np.einsum('hab,bc,hcd,da', T, Us['x'].T @ Us['y'], D, Us['wp'].T @ Us['x']))).real
    e = [0.2, 0.1, 0.05, 0.025]
    pr('C5  weak target V = exp(i eps a.t): ratios v(eps)/v(eps/2) =', ' '.join('%.3f' % (val(a) / val(a / 2)) for a in e), '(8 = eps^3)')

# SU(2): d = 0
t2 = gens(2); f2, d2 = fd(t2)
pr('--- SU(2): max|d^{abc}| = %.1e' % np.abs(d2).max())

# ---------------- C6: operator identity with charges on sites
def site_ops(t, nsites):
    N = t.shape[1]; ops = []
    for s in range(nsites):
        mats = []
        for a in range(t.shape[0]):
            m = np.array([[1.0 + 0j]])
            for s2 in range(nsites): m = np.kron(m, t[a] if s2 == s else np.eye(N))
            mats.append(m)
        ops.append(np.array(mats))
    return ops
def sym3(A, B, Cc):
    return sum(X @ Y @ Z for X, Y, Z in itertools.permutations([A, B, Cc])) / 6
def sym2(A, B): return (A @ B + B @ A) / 2
for N, nsites in ((2, 3), (3, 3)):
    t = gens(N); f, d = fd(t); n = N * N - 1; rho = site_ops(t, nsites); dim = N ** nsites
    worst = 0
    for (s1, s2, s3) in itertools.product(range(nsites), repeat=3):
        Cabc = rng.normal(size=(n, n, n))
        d12, d13 = float(s1 == s2), float(s1 == s3)
        lhs = np.zeros((dim, dim), complex); rhs = np.zeros((dim, dim), complex)
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    cc = Cabc[a, b, c]
                    lhs += cc * rho[s1][a] @ (rho[s2][b] @ rho[s3][c] + rho[s3][c] @ rho[s2][b])
                    rhs += 2 * cc * sym3(rho[s1][a], rho[s2][b], rho[s3][c])
        # contact terms: i C [f^{abr} d12 {rho^r(1) rho^c(3)}sym + f^{acr} d13 {rho^r(1) rho^b(2)}sym]
        X1 = np.einsum('abc,abr->rc', Cabc, f); X2 = np.einsum('abc,acr->rb', Cabc, f)
        for r in range(n):
            for c in range(n):
                rhs += 1j * d12 * X1[r, c] * sym2(rho[s1][r], rho[s3][c]) + 1j * d13 * X2[r, c] * sym2(rho[s1][r], rho[s2][c])
        # one-rho term: -(1/6) C [f^{abr} f^{rcs} + f^{acr} f^{rbs}] d12 d13 rho^s(1)
        Y = np.einsum('abc,abr,rcs->s', Cabc, f, f) + np.einsum('abc,acr,rbs->s', Cabc, f, f)
        rhs += -(1 / 6) * d12 * d13 * np.einsum('s,sij->ij', Y, rho[s1])
        worst = max(worst, np.abs(lhs - rhs).max() / np.abs(lhs).max())
    pr('C6  SU(%d), %d sites, all %d assignments of (x\',x,y) to sites, random C^{abc}: rho^a{rho^b,rho^c} = 2{..}_sym + contact + one-rho,'
       ' max rel dev %.1e' % (N, nsites, nsites ** 3, worst))

# ---------------- C7: single quark, colour-averaged state rho = 1/N
for N in (3, 4):
    t = gens(N); f, d = fd(t); n = N * N - 1
    S = np.zeros((n, n, n), complex)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                S[a, b, c] = np.trace(sym3(t[a], t[b], t[c])) / N
    o = np.einsum('abc,abc', d, S)
    pr('C7  SU(%d) single quark: <{t^a t^b t^c}_sym> - d^{abc} N/((N^2-4)(N^2-1)) o : %.1e ; o = %.6f = (N^2-4)(N^2-1)/(4N^2) = %.6f ; f-part f^{abc}S^{abc} = %.1e'
       % (N, np.abs(S - d * N / ((N * N - 4) * (N * N - 1)) * o).max(), o.real, (N * N - 4) * (N * N - 1) / (4 * N * N), abs(np.einsum('abc,abc', f, S))))
open(__file__.replace('.py', '.log'), 'w').write('\n'.join(out) + '\n')
