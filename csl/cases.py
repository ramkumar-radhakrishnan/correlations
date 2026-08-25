"""Explicit Pauli-string Hamiltonians for the requested (n_q, N) cases."""
from fractions import Fraction
import numpy as np
from csl_lattice import (kron, clock, pi_squared, embed, pauli_decompose)

COUPL = ["t", "J", "h", "k"]     # t, J, h, kappa


def pieces(n_q, N, bc="open"):
    """Return the four coupling-independent operators {t,J,h,k} -> matrix."""
    Q = 2 ** n_q
    dim = Q ** N
    Z = clock(Q); Zd = Z.conj().T; P2 = pi_squared(Q)
    out = {c: np.zeros((dim, dim), dtype=complex) for c in COUPL}
    for n in range(1, N + 1):
        out["t"] += embed({n: P2}, N, Q)
        out["h"] += embed({n: np.eye(Q) - 0.5 * (Z + Zd)}, N, Q)
    bonds = [(n, n + 1) for n in range(1, N)]
    if bc == "periodic":
        bonds.append((N, 1))
    for (n, m) in bonds:
        ZZ = embed({m: Z, n: Zd}, N, Q)
        out["J"] += np.eye(dim) - 0.5 * (ZZ + ZZ.conj().T)
        out["k"] += -(ZZ - ZZ.conj().T) / (2j)
    return out


def label(s, n_q):
    if set(s) == {"I"}:
        return "I"
    parts = []
    for q, ch in enumerate(s):
        if ch == "I":
            continue
        n, p = divmod(q, n_q)
        parts.append(f"{ch}{{{n+1},{n_q-1-p}}}" if n_q > 1 else f"{ch}{{{n+1}}}")
    return " ".join(parts)


def frac(x):
    f = Fraction(x).limit_denominator(64)
    assert abs(float(f) - x) < 1e-9, x
    if f.denominator == 1:
        return f"{f.numerator:+d}"
    return f"{f.numerator:+d}/{f.denominator}"


def report(n_q, N, bc="open"):
    Q, nq_tot = 2 ** n_q, n_q * N
    print("=" * 74)
    print(f"  n_q = {n_q} qubit(s)/site,  Q = {Q},  N = {N} sites,  "
          f"{nq_tot} qubits,  dim = {Q**N},  BC = {bc}")
    print("=" * 74)
    P = pieces(n_q, N, bc)
    total = {}
    for c in COUPL:
        d = pauli_decompose(P[c], nq_tot)
        if not d:
            continue
        rows = sorted(d.items(), key=lambda kv: (kv[0].count("I") * -1, kv[0]))
        print(f"\n  -- {c} sector ({len(rows)} Pauli terms) --")
        for s, v in rows:
            print(f"     {frac(v.real):>8s} {c}   {label(s, n_q)}")
            total[s] = total.get(s, 0) + 1
    print(f"\n  distinct Pauli strings in H: {len(total)}"
          f"   (identity + {len(total)-1} non-trivial)")
