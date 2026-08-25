"""Low-lying spectra and the CSL order parameter for the three cases."""
import numpy as np
from csl_lattice import clock, embed, hamiltonian


def winding_op(n_q, N, bc="open"):
    """(1/N_b) sum_bonds sin(phi_{n+1} - phi_n)  -- the CSL / pion-current density."""
    Q = 2 ** n_q
    Z = clock(Q); dim = Q ** N
    bonds = [(n, n + 1) for n in range(1, N)] + ([(N, 1)] if bc == "periodic" else [])
    W = np.zeros((dim, dim), dtype=complex)
    for (n, m) in bonds:
        ZZ = embed({m: Z, n: Z.conj().T}, N, Q)
        W += (ZZ - ZZ.conj().T) / (2j)
    return W / len(bonds)


def scan(n_q, N, kappas, t=0.5, J=1.0, h=1.0, bc="open", nlev=4):
    print(f"\n  n_q={n_q}  N={N}  ({n_q*N} qubits, dim {2**(n_q*N)})  "
          f"t={t} J={J} h={h}  {bc} BC")
    print("    kappa |      E0        E1        E2        E3   |  <sin dphi>   gap")
    print("    ------+-----------------------------------------+-------------------")
    W = winding_op(n_q, N, bc)
    for k in kappas:
        H = hamiltonian(n_q, N, t, J, h, k, bc)
        ev, V = np.linalg.eigh(H)
        g = V[:, 0]
        w = (g.conj() @ W @ g).real
        lev = "  ".join(f"{e:8.4f}" for e in ev[:nlev])
        print(f"    {k:5.2f} | {lev} | {w:10.5f}  {ev[1]-ev[0]:8.4f}")


if __name__ == "__main__":
    ks = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
    scan(1, 2, ks)          # the note's original case
    scan(1, 4, ks)
    scan(2, 2, ks)
    scan(2, 4, ks)
