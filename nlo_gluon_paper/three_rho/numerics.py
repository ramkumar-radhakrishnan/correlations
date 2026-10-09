"""Numerical colour algebra for SU(N): structure constants, adjoint Wilson lines,
evaluation of fully contracted colour factors."""
import numpy as np
import itertools
import string


def gellmann_like(N):
    """Hermitian generators t^a of SU(N), Tr t^a t^b = delta/2."""
    ts = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex); m[i, j] = m[j, i] = 0.5; ts.append(m)
            m = np.zeros((N, N), complex); m[i, j] = -0.5j; m[j, i] = 0.5j; ts.append(m)
    for k in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(k):
            m[i, i] = 1
        m[k, k] = -k
        m *= 1 / np.sqrt(2 * k * (k + 1))
        ts.append(m)
    return ts


def structure_constants(N):
    ts = gellmann_like(N)
    n = len(ts)
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            com = ts[a] @ ts[b] - ts[b] @ ts[a]
            for c in range(n):
                f[a, b, c] = np.real(-2j * np.trace(com @ ts[c]))
    return ts, f


def random_su(N, rng):
    z = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diag(r) / np.abs(np.diag(r))
    q = q * d
    q = q / np.linalg.det(q) ** (1 / N)
    return q


def adjoint(V, ts):
    n = len(ts)
    U = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            U[a, b] = np.real(2 * np.trace(ts[a] @ V @ ts[b] @ V.conj().T))
    return U


def random_adjoint(N, ts, rng, symmetric):
    if symmetric:
        # V = W diag(1,-1,-1,...) W^dagger with det = 1 needs an even number of -1's
        W = random_su(N, rng)
        dg = np.ones(N)
        if N == 2:
            # rotation by pi in SO(3): V = i n.sigma, V^2 = -1 (centre)
            return adjoint(W @ np.diag([1j, -1j]) @ W.conj().T, ts)
        dg[1] = dg[2] = -1
        return adjoint(W @ np.diag(dg) @ W.conj().T, ts)
    return adjoint(random_su(N, rng), ts)


def eval_colour(facs, f, Umats, rhovecs, ncol):
    """Evaluate a fully contracted colour factor.  rho's are classical vectors."""
    letters = {}
    pool = iter(string.ascii_letters)
    ops, arrs = [], []
    for fac in facs:
        idx = [fac[1], fac[2]] if fac[0] in ('U', 'd') else ([fac[1], fac[2], fac[3]] if fac[0] == 'f' else [fac[1]])
        sub = ''
        for x in idx:
            if x not in letters:
                letters[x] = next(pool)
            sub += letters[x]
        ops.append(sub)
        if fac[0] == 'U':
            arrs.append(Umats[fac[3]])
        elif fac[0] == 'f':
            arrs.append(f)
        elif fac[0] == 'd':
            arrs.append(np.eye(ncol))
        elif fac[0] == 'rho':
            arrs.append(rhovecs[fac[2]])
    expr = ','.join(ops) + '->'
    return np.einsum(expr, *arrs, optimize='greedy') if arrs else 1.0
