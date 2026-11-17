# sanity: symmetrisation formulas used in nc.py, tested with su(2) spin-2 matrices (Lie algebra = 3-dim subspace of 25-dim)
import numpy as np, itertools
rng = np.random.default_rng(1)
j = 2; d = 2 * j + 1
m = np.arange(j, -j - 1, -1)
Jz = np.diag(m).astype(complex)
Jp = np.zeros((d, d), complex)
for k in range(1, d): Jp[k - 1, k] = np.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
Jx = (Jp + Jp.conj().T) / 2; Jy = (Jp - Jp.conj().T) / 2j
G = [Jx, Jy, Jz]
def rnd(): return sum((rng.normal() + 1j * rng.normal()) * g for g in G)
def com(a, b): return a @ b - b @ a
def S(*ops):
    n = len(ops); tot = 0
    for p in itertools.permutations(range(n)):
        t = np.eye(d, dtype=complex)
        for i in p: t = t @ ops[i]
        tot = tot + t
    import math
    return tot / math.factorial(n)
def Z3(a, b, c): return 0.25 * com(com(a, b), c) + com(a, com(b, c)) / 12 + com(b, com(a, c)) / 12
def in_alg(R):
    # distance of R from span(G)
    B = np.array([g.flatten() for g in G]).T
    c, *_ = np.linalg.lstsq(B, R.flatten(), rcond=None)
    return np.linalg.norm(B @ c - R.flatten())
A, B, C, D = rnd(), rnd(), rnd(), rnd()
# n=3
R3 = A @ B @ C - (S(A, B, C) + 0.5 * S(com(A, B), C) + 0.5 * S(com(A, C), B) + 0.5 * S(com(B, C), A))
print("n=3: remainder - Z3 =", np.linalg.norm(R3 - Z3(A, B, C)), " |remainder| =", np.linalg.norm(R3))
# n=4
X = [A, B, C, D]
tot = S(*X)
for i, k in itertools.combinations(range(4), 2):
    r = [q for q in range(4) if q not in (i, k)]
    tot = tot + 0.5 * S(com(X[i], X[k]), X[r[0]], X[r[1]])
two = 0
for (i, k), (p, q) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
    two = two + S(0.5 * com(X[i], X[k]), 0.5 * com(X[p], X[q]))
for l in range(4):
    i, k, p = [q for q in range(4) if q != l]
    two = two + S(Z3(X[i], X[k], X[p]), X[l])
R4 = A @ B @ C @ D - tot - two
print("n=4: remainder after 4,3,2-operator terms: |R| =", np.linalg.norm(R4), " distance from Lie algebra =", in_alg(R4),
      " (without the two-operator terms:", in_alg(A @ B @ C @ D - tot), ")")
