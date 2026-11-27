"""Does the omitted term  ΔA^c(w) = [ N̄2 , Ā1^{†(1)c}(w) ]  contribute to the cross section
   δσ = Σ_{c,i} < A_LO^{c}_i(w')^† ΔA^c_i(w) > + c.c.   for colour-singlet projectiles?
Exact colour algebra with explicit charges (quarks: t^a, antiquarks: -t^{a*}, gluons: adjoint T^a)."""
import numpy as np, itertools
from functools import reduce
from lattice import su_n_generators, structure_constants, adjoint, random_su_n
rng = np.random.default_rng(3)
t = su_n_generators(3); f = structure_constants(3 and t); TA = -1j*f    # adjoint generators
def charges(kinds):
    dims = [3 if k in 'qa' else 8 for k in kinds]
    Q = []
    for s, k in enumerate(kinds):
        gen = t if k == 'q' else ([-g.conj() for g in t] if k == 'a' else list(TA))
        Q.append([reduce(np.kron, [gen[a] if r == s else np.eye(dims[r]) for r in range(len(kinds))]) for a in range(8)])
    return Q, int(np.prod(dims))
def singlet(Q, dim):
    C = sum(sum(Q[s][a] for s in range(len(Q))) @ sum(Q[s][a] for s in range(len(Q))) for a in range(8))
    w, v = np.linalg.eigh(C); return v[:, np.argmin(np.abs(w))]
def K(x, w): d = x - w; return d / (d @ d)
def run(kinds, trials=3):
    Q, dim = charges(kinds); S = singlet(Q, dim); ns = len(kinds)
    out = []
    for _ in range(trials):
        X = rng.normal(size=(ns, 2)); w, wp = rng.normal(size=2), rng.normal(size=2)
        U = {s: adjoint(random_su_n(3, rng), t) for s in range(ns)}; Uw, Uwp = (adjoint(random_su_n(3, rng), t) for _ in range(2))
        Qbar = [[sum(U[s][a, b] * Q[s][b] for b in range(8)) for a in range(8)] for s in range(ns)]
        G = np.array([[0.0 if r == s else 0.37 * np.log(1/np.linalg.norm(X[r]-X[s])) for s in range(ns)] for r in range(ns)])
        Nbar = sum(G[r, s] * Qbar[r][a] @ Qbar[s][a] for r in range(ns) for s in range(ns) for a in range(8))
        N    = sum(G[r, s] * Q[r][a] @ Q[s][a]       for r in range(ns) for s in range(ns) for a in range(8))
        def Abar1(c, i, ww): return sum(K(X[s], ww)[i] * Qbar[s][c] for s in range(ns))          # Ā1^†(1)
        def A1(c, i, ww):    return sum(K(X[s], ww)[i] * Q[s][c] for s in range(ns))             # A1^†(1)
        def ALO(c, i, ww, Uww): return Abar1(c, i, ww) - sum(Uww[c, d] * A1(d, i, ww) for d in range(8))
        dsig = ref = 0
        for c, i in itertools.product(range(8), range(2)):
            a_lo = ALO(c, i, wp, Uwp)
            dA = Nbar @ Abar1(c, i, w) - Abar1(c, i, w) @ Nbar
            ref_term = Abar1(c, i, w) @ N                          # a term you did keep, (IV.18): Ā1^† N2
            dsig += S.conj() @ a_lo.conj().T @ dA @ S
            ref += S.conj() @ a_lo.conj().T @ ref_term @ S
        out.append((dsig + np.conj(dsig), ref + np.conj(ref)))
    return out
for name, kinds in [('q qbar dipole', 'qa'), ('q qbar g', 'qag'), ('baryon q q q', 'qqq')]:
    for d, r in run(kinds):
        print(f'{name:14s}: δσ from omitted commutator = {d.real:+.3e}   (for scale: a kept N2 term = {r.real:+.3e})')
