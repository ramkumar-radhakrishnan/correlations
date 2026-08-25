"""
Digitized chiral-soliton-lattice (CSL) Hamiltonian on a 1D spatial lattice.

Generalizes the N = 2 site, Q = 2 (one qubit/site) derivation of
`CSL_quantumcomp.pdf` to arbitrary  n_q  qubits per site  (Q = 2**n_q  field
values) and arbitrary  N  sites.

Continuum starting point (one spatial dimension):

    E = f^2 \int dx [ Pi^2/(2 f^4) + (1/2)(d_x phi)^2
                      + m^2 (1 - cos phi) - H d_x phi ]

Compactified lattice form (spacing a, sites n = 1..N):

    H = t Sum_n Pi_n^2
      + J Sum_bonds (1 - cos(phi_{n+1} - phi_n))
      + h Sum_n (1 - cos phi_n)
      - kappa Sum_bonds sin(phi_{n+1} - phi_n)

    t     = 1/(2 f^2 a)      kinetic
    J     = f^2/a            gradient / bond
    h     = f^2 m^2 a        pion mass
    kappa = f^2 H = mu B/(4 pi^2)   magnetic (Wess-Zumino) term

Digitization: phi in [0, 2pi) -> phi_k = 2 pi k / Q, k = 0..Q-1, and the
site register stores k in binary, k = sum_j 2^j b_j.

Basis / operator ordering conventions
-------------------------------------
  * full ket = |site 1> (x) |site 2> (x) ... (x) |site N>   (site 1 leftmost)
  * within a site, |b_{n_q-1} ... b_1 b_0>  (MSB leftmost), so the site index
    k is the plain binary reading of the register.
  * qubit label (n, j) = site n (1-based), bit j (0 = LSB).
"""

import itertools
import numpy as np

# ---------------------------------------------------------------------------
# single-qubit Paulis
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PY = np.array([[0, -1j], [1j, 0]], dtype=complex)
PZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": PX, "Y": PY, "Z": PZ}


def kron(mats):
    out = np.array([[1.0 + 0j]])
    for m in mats:
        out = np.kron(out, m)
    return out


# ---------------------------------------------------------------------------
# single-site (Q-dimensional) clock / shift operators and Pi^2
# ---------------------------------------------------------------------------
def clock(Q):
    """Z_Q |k> = omega^k |k>,  omega = exp(2 pi i / Q).   (Z_Q = exp(i phi))"""
    w = np.exp(2j * np.pi / Q)
    return np.diag(w ** np.arange(Q))


def shift(Q):
    """X_Q |k> = |k+1 mod Q>  (the increment / ladder operator)."""
    X = np.zeros((Q, Q), dtype=complex)
    for k in range(Q):
        X[(k + 1) % Q, k] = 1.0
    return X


def brillouin_zone(Q):
    """Integer rotor momenta in DFT-index order j = l mod Q.

    Symmetric window l in {-Q/2, ..., Q/2 - 1} for even Q (matches eq. (1.11)
    of the note, which uses l in {0, -1} for Q = 2).
    """
    return np.array([j if j < Q / 2 else j - Q for j in range(Q)])


def pi_squared(Q):
    """Pi^2 = sum_l l^2 |l><l| in the *position* (computational) basis."""
    l = brillouin_zone(Q)
    F = np.exp(2j * np.pi * np.outer(np.arange(Q), l) / Q) / np.sqrt(Q)  # <k|l>
    return F @ np.diag(l.astype(float) ** 2) @ F.conj().T


def pi_squared_coeffs(Q):
    """Closed-form expansion  Pi^2 = sum_{m=0}^{Q-1} c_m X_Q^m .

        c_0 = (Q^2 + 2)/12
        c_m = (-1)^m / (2 sin^2(pi m / Q))      (m != 0)
    """
    c = np.zeros(Q)
    c[0] = (Q ** 2 + 2) / 12.0
    for m in range(1, Q):
        c[m] = (-1) ** m / (2.0 * np.sin(np.pi * m / Q) ** 2)
    return c


# ---------------------------------------------------------------------------
# many-site embedding
# ---------------------------------------------------------------------------
def embed(op_by_site, N, Q):
    """Place per-site Q x Q operators (dict site->matrix, sites 1-based)."""
    return kron([op_by_site.get(n, np.eye(Q, dtype=complex)) for n in range(1, N + 1)])


def hamiltonian(n_q, N, t=1.0, J=1.0, h=1.0, kappa=0.0, bc="open"):
    """Full digitized CSL Hamiltonian, dimension Q**N = 2**(N n_q)."""
    Q = 2 ** n_q
    dim = Q ** N
    Z, P2 = clock(Q), pi_squared(Q)
    Zd = Z.conj().T
    H = np.zeros((dim, dim), dtype=complex)

    # kinetic + mass (single site)
    for n in range(1, N + 1):
        H += t * embed({n: P2}, N, Q)
        H += h * embed({n: np.eye(Q) - 0.5 * (Z + Zd)}, N, Q)

    # bonds
    bonds = [(n, n + 1) for n in range(1, N)]
    if bc == "periodic":
        bonds.append((N, 1))
    for (n, m) in bonds:
        ZZ = embed({m: Z, n: Zd}, N, Q)          # exp(i(phi_m - phi_n))
        H += J * (np.eye(dim) - 0.5 * (ZZ + ZZ.conj().T))
        H += -kappa * (ZZ - ZZ.conj().T) / (2j)
    return H


# ---------------------------------------------------------------------------
# Pauli decomposition
# ---------------------------------------------------------------------------
def pauli_decompose(M, nqubits, tol=1e-10):
    """Fast Pauli decomposition. Returns {pauli-string: coeff}, qubit 0 leftmost.

    c_P = Tr(P M)/2^n is computed as a sequence of local 4x4 basis changes on
    the (row, col) index pair of each qubit -- O(4^n) instead of O(8^n).
    """
    n = nqubits
    T = np.zeros((4, 2, 2), dtype=complex)          # T[p, j, i] = P_p[j, i]/2
    for p, ch in enumerate("IXYZ"):
        T[p] = PAULI[ch].T / 2.0                    # trace(P M) = sum_ij P[j,i] M[i,j]
    A = M.reshape([2] * (2 * n))                    # i_0..i_{n-1}, j_0..j_{n-1}
    A = np.transpose(A, [k for q in range(n) for k in (q, n + q)])
    A = A.reshape([4] * n)                          # per-qubit (i,j) composite index
    for q in range(n):
        A = np.tensordot(T.reshape(4, 4), A, axes=([1], [q]))
        A = np.moveaxis(A, 0, q)
    out = {}
    it = np.nditer(A, flags=["multi_index"])
    for v in it:
        if abs(v) > tol:
            out["".join("IXYZ"[k] for k in it.multi_index)] = complex(v)
    return out


def pretty_pauli(terms, n_q, N, tol=1e-10):
    """Render a Pauli dict with (site, bit) labels, e.g. 'X_{2,0} Z_{3,1}'."""
    lines = []
    for s, c in sorted(terms.items(), key=lambda kv: (kv[0].count("I") * -1, kv[0])):
        if abs(c) < tol:
            continue
        if set(s) == {"I"}:
            name = "I"
        else:
            parts = []
            for q, ch in enumerate(s):
                if ch == "I":
                    continue
                n, j = divmod(q, n_q)
                parts.append(f"{ch}_{{{n+1},{n_q-1-j}}}")
            name = " ".join(parts)
        coef = c.real if abs(c.imag) < tol else c
        lines.append((name, coef))
    return lines
