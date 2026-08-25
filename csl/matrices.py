"""Explicit symbolic matrices (small cases) and spectra / CSL observables."""
from fractions import Fraction
import numpy as np
from csl_lattice import clock, pi_squared, embed, hamiltonian
from cases import pieces, COUPL


def sym(x, name):
    f = Fraction(x).limit_denominator(64)
    if f == 0:
        return ""
    if f == 1:
        return f"+{name}"
    if f == -1:
        return f"-{name}"
    if f.denominator == 1:
        return f"{f.numerator:+d}{name}"
    return f"{f.numerator:+d}{name}/{f.denominator}"


def symbolic_matrix(n_q, N, bc="open"):
    P = pieces(n_q, N, bc)
    Q, dim = 2 ** n_q, (2 ** n_q) ** N
    ent = [["" for _ in range(dim)] for _ in range(dim)]
    for c in COUPL:
        M = P[c].real
        for i in range(dim):
            for j in range(dim):
                ent[i][j] += sym(M[i, j], c)
    labels = []
    for idx in range(dim):
        ks, r = [], idx
        for _ in range(N):
            ks.append(r // Q ** (N - 1 - len(ks)))
            r = r % Q ** (N - 1 - len(ks) + 1)
        ks = [(idx // Q ** (N - 1 - s)) % Q for s in range(N)]
        labels.append("".join(str(k) for k in ks))
    return labels, [[e if e else "." for e in row] for row in ent]


def show(n_q, N, bc="open", width=None):
    labels, M = symbolic_matrix(n_q, N, bc)
    w = width or max(len(e) for row in M for e in row) + 1
    hdr = " " * (len(labels[0]) + 2) + "".join(f"{l:>{w}s}" for l in labels)
    print(f"\n  n_q={n_q}, N={N}, {bc} BC   basis |k_1...k_N>, "
          f"phi = 2*pi*k/{2**n_q}")
    print(hdr)
    for l, row in zip(labels, M):
        print(f"  {l}  " + "".join(f"{e:>{w}s}" for e in row))
