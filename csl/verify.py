"""Checks of the generalized formulas against the N=2, Q=2 derivation."""
import numpy as np
from csl_lattice import (I2, PX, PY, PZ, kron, clock, shift, pi_squared,
                         pi_squared_coeffs, hamiltonian, pauli_decompose)

ok = lambda name, c: print(f"[{'PASS' if c else 'FAIL':4}] {name}")

# ---- 1. Pi^2 closed form  c_0=(Q^2+2)/12, c_m=(-1)^m/(2 sin^2(pi m/Q)) -----
for Q in (2, 4, 8, 16):
    X = shift(Q); c = pi_squared_coeffs(Q)
    rec = sum(c[m] * np.linalg.matrix_power(X, m) for m in range(Q))
    ok(f"Q={Q:2d}: Pi^2 = sum_m c_m X^m", np.allclose(rec, pi_squared(Q)))
    ok(f"Q={Q:2d}: spectrum of Pi^2 = l^2", np.allclose(
        np.sort(np.linalg.eigvalsh(pi_squared(Q))),
        np.sort(np.array([j if j < Q/2 else j-Q for j in range(Q)])**2.0)))

# ---- 2. Q=2 reproduces eqs (1.10), (1.15), (1.17) --------------------------
Z2 = clock(2)
ok("Q=2: 1-cos phi = I - Z",            np.allclose(np.eye(2)-0.5*(Z2+Z2.conj().T), I2-PZ))
ok("Q=2: Pi^2 = (1/2)(I - X)  [eq 1.15]", np.allclose(pi_squared(2), 0.5*(I2-PX)))
ok("Q=2: sin(dphi) = 0        [eq 1.18]", np.allclose(
    (kron([Z2.conj().T, Z2]) - kron([Z2, Z2.conj().T]))/(2j), np.zeros((4,4))))
ok("Q=2: cos(dphi) = Z x Z    [eq 1.17]", np.allclose(
    0.5*(kron([Z2.conj().T,Z2])+kron([Z2,Z2.conj().T])), kron([PZ,PZ])))

# ---- 3. full N=2, Q=2 Hamiltonian vs eq (1.21) -----------------------------
f, a, m = 1.3, 0.7, 0.9
t, J, h = 1/(2*f**2*a), f**2/a, f**2*m**2*a
q = 1/(4*f**2*a)
paper = np.array([
    [2*q,        -q,             -q,            0     ],
    [-q,         2*q+2*h+2*J,     0,           -q     ],
    [-q,          0,             2*q+2*h+2*J,  -q     ],
    [0,          -q,             -q,            2*q+4*h]], dtype=complex)
Hgen = hamiltonian(n_q=1, N=2, t=t, J=J, h=h, kappa=0.0, bc="open")
ok("N=2, Q=2 Hamiltonian == eq (1.21)", np.allclose(Hgen, paper))

# ---- 4. Q=4 operator identities -------------------------------------------
Z4, X4 = clock(4), shift(4)
Z1, Z0 = kron([PZ, I2]), kron([I2, PZ])   # (MSB, LSB) = (bit1, bit0)
X1, X0 = kron([PX, I2]), kron([I2, PX])
ok("Q=4: cos phi = (Z1 + Z1 Z0)/2", np.allclose(0.5*(Z4+Z4.conj().T), 0.5*(Z1+Z1@Z0)))
ok("Q=4: sin phi = (Z1 - Z1 Z0)/2", np.allclose((Z4-Z4.conj().T)/(2j), 0.5*(Z1-Z1@Z0)))
ok("Q=4: X + X^dag = X0 + X1 X0",   np.allclose(X4+X4.conj().T, X0 + X1@X0))
ok("Q=4: X^2 = X1",                 np.allclose(X4@X4, X1))
ok("Q=4: Pi^2 = 3/2 I - X0 - X1 X0 + X1/2", np.allclose(
    pi_squared(4), 1.5*np.eye(4) - X0 - X1@X0 + 0.5*X1))

# Q=4 bond identities:  cos d = (A a + A B a b)/2 ,  sin d = A a (b - B)/2
A, B = kron([PZ,I2,I2,I2]), kron([I2,PZ,I2,I2])   # site n+1: bit1, bit0
aa, bb = kron([I2,I2,PZ,I2]), kron([I2,I2,I2,PZ]) # site n  : bit1, bit0
ZZ = kron([Z4, np.eye(4)]) @ kron([np.eye(4), Z4.conj().T])   # exp(i(phi_{n+1}-phi_n))
ok("Q=4: cos dphi = (Z11 Z21 + Z11 Z10 Z21 Z20)/2",
   np.allclose(0.5*(ZZ+ZZ.conj().T), 0.5*(A@aa + A@B@aa@bb)))
ok("Q=4: sin dphi = Z11 Z21 (Z20 - Z10)/2",
   np.allclose((ZZ-ZZ.conj().T)/(2j), 0.5*(A@aa@(bb-B))))

# ---- 5. Q=4 reduces to Q=2 physics on the {0, pi} subspace -----------------
ok("Q=4 mass term diag = m^2 a (0,1,2,1)", np.allclose(
    np.diag(np.eye(4)-0.5*(Z4+Z4.conj().T)).real, [0,1,2,1]))

# ---- 6. Hermiticity / reality of all requested cases -----------------------
for nq, N in [(1,4),(2,2),(2,4)]:
    H = hamiltonian(nq, N, t=t, J=J, h=h, kappa=0.31, bc="open")
    ok(f"n_q={nq}, N={N}: Hermitian, dim={H.shape[0]}", np.allclose(H, H.conj().T))
    terms = pauli_decompose(H, nq*N)
    im = max(abs(c.imag) for c in terms.values())
    ok(f"n_q={nq}, N={N}: all Pauli coefficients real (max Im={im:.1e})", im < 1e-10)

# ---- 7. general product formula for Z_Q = exp(i phi) ----------------------
#   Z_Q = e^{i pi (Q-1)/Q} prod_j [ cos(pi 2^j/Q) I - i sin(pi 2^j/Q) Z_j ]
for n_q in (1, 2, 3, 4):
    Q = 2 ** n_q
    prod = np.eye(Q, dtype=complex)
    for j in range(n_q):
        th = np.pi * 2 ** j / Q
        Zj = kron([PZ if p == n_q - 1 - j else I2 for p in range(n_q)])
        prod = prod @ (np.cos(th) * np.eye(Q) - 1j * np.sin(th) * Zj)
    prod *= np.exp(1j * np.pi * (Q - 1) / Q)
    ok(f"n_q={n_q}: Z_Q = e^{{i pi (Q-1)/Q}} prod_j [c I - i s Z_j]",
       np.allclose(prod, clock(Q)))
    # tensor-product (phase-gate) form
    tp = kron([np.diag([1, np.exp(2j*np.pi*2**(n_q-1-p)/Q)]) for p in range(n_q)])
    ok(f"n_q={n_q}: Z_Q = (x)_j diag(1, omega^{{2^j}})", np.allclose(tp, clock(Q)))

# ---- 8. classical / continuum limits ---------------------------------------
#   Pi^2 diagonal element -> Q^2/12 + 1/6  (Fornberg spectral 2nd-derivative)
for Q in (4, 8, 16, 32):
    ok(f"Q={Q:2d}: <k|Pi^2|k> = Q^2/12 + 1/6",
       abs(pi_squared(Q)[0, 0].real - (Q**2/12 + 1/6)) < 1e-9)

# ---- 9. closed-form Pauli coefficients of cos phi and sin phi --------------
#   cos phi = sum_{S subset {0..n_q-2}} sin(pi/Q + |S| pi/2) A_S Z_{n_q-1} prod_{j in S} Z_j
#   sin phi = sum_S           cos(pi/Q + |S| pi/2) A_S  (same string)
#   A_S = prod_{j in S} sin(pi 2^j/Q) * prod_{j not in S} cos(pi 2^j/Q)
import itertools as _it
for n_q in (1, 2, 3, 4, 5):
    Q = 2 ** n_q
    C = np.zeros((Q, Q), dtype=complex); S_ = np.zeros((Q, Q), dtype=complex)
    rest = list(range(n_q - 1))
    for r in range(len(rest) + 1):
        for Sset in _it.combinations(rest, r):
            A = 1.0
            for j in rest:
                A *= np.sin(np.pi*2**j/Q) if j in Sset else np.cos(np.pi*2**j/Q)
            ops = [PZ if (p == 0 or (n_q-1-p) in Sset) else I2 for p in range(n_q)]
            P = kron(ops)
            C += np.sin(np.pi/Q + r*np.pi/2) * A * P
            S_ += np.cos(np.pi/Q + r*np.pi/2) * A * P
    Zc = clock(Q)
    ok(f"n_q={n_q}: closed-form Pauli series for cos phi",
       np.allclose(C, 0.5*(Zc + Zc.conj().T)))
    ok(f"n_q={n_q}: closed-form Pauli series for sin phi",
       np.allclose(S_, (Zc - Zc.conj().T)/(2j)))

# ---- 10. Pauli-string counts ----------------------------------------------
for n_q in (1, 2, 3, 4):
    Q = 2 ** n_q
    ok(f"n_q={n_q}: Pi^2 has 3^(n_q-1)+1 = {3**(n_q-1)+1} Pauli strings",
       len(pauli_decompose(pi_squared(Q), n_q)) == 3**(n_q-1) + 1)
    ok(f"n_q={n_q}: cos phi has Q/2 = {Q//2} Pauli-Z strings",
       len(pauli_decompose(0.5*(clock(Q)+clock(Q).conj().T), n_q)) == Q // 2)
