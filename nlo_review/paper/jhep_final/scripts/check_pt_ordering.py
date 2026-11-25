#!/usr/bin/env python3
"""Check of Eqs. (pt2) and (pt3) of Section 2.5: perturbation theory with operator-valued
(non-commuting) matrix elements V_nm, including the order of the factors in the normalisation terms.

Model: L soft levels with energies E_n (numbers), valence space of dimension d; V_nm are d x d blocks.
The exact unitary Omega with (i) Omega^dag H Omega block diagonal and (ii) P_n Omega P_n Hermitian
positive is constructed from the exact spectral projectors and compared with the series.
The error of the series truncated at order N must scale as g^(N+1).
In Eq. (pt3) the components along excited levels n>=1 are of third-order accuracy; the vacuum component
is given to second order (in the model a term V_0n V_nm V_m0 of third order exists; in the CGC it starts at g^4).
Run:  python3 -I check_pt_ordering.py
"""
import numpy as np
from scipy.linalg import sqrtm, eigh
rng = np.random.default_rng(7)
L, d = 5, 3
E = np.array([0.0, 1.0, 2.3, 3.1, 4.7])
def herm(n):
    a = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)); return a + a.conj().T
V = herm(L*d)
def blk(M, n, m): return M[n*d:(n+1)*d, m*d:(m+1)*d]
V[0:d, 0:d] = 0                      # V_00 = 0
V2 = V.copy()
for n in range(L): V2[n*d:(n+1)*d, n*d:(n+1)*d] = 0   # all diagonal blocks zero (for pt2 with generic i)

def exact_omega(Vm, g):
    H = np.kron(np.diag(E), np.eye(d)) + g*Vm
    w, U = eigh(H)
    Om = np.zeros((L*d, L*d), complex)
    for n in range(L):
        idx = np.argsort(np.abs(w - E[n]))[:d]
        Q = U[:, idx] @ U[:, idx].conj().T
        P = np.zeros((L*d, L*d)); P[n*d:(n+1)*d, n*d:(n+1)*d] = np.eye(d)
        S = P @ Q @ P
        Sb = S[n*d:(n+1)*d, n*d:(n+1)*d]
        inv = np.zeros((L*d, L*d), complex); inv[n*d:(n+1)*d, n*d:(n+1)*d] = np.linalg.inv(sqrtm(Sb))
        Om += Q @ P @ inv
    return Om

def series2(Vm, g, i):
    """<n|Psi_i> blocks from Eq. (pt2)."""
    out = {}
    W = g*Vm
    others = [n for n in range(L) if n != i]
    norm = np.eye(d, dtype=complex)
    for m in others: norm -= 0.5*blk(W, i, m) @ blk(W, m, i)/(E[i]-E[m])**2
    out[i] = norm
    for n in others:
        t = blk(W, n, i)/(E[i]-E[n])
        for m in others: t = t + blk(W, n, m) @ blk(W, m, i)/((E[i]-E[n])*(E[i]-E[m]))
        out[n] = t
    return out

def series3(Vm, g):
    """<n|Psi_0> blocks from Eq. (pt3); swap=True puts the norm factor to the LEFT (wrong order)."""
    res = {}
    for swap in (False, True):
        W = g*Vm; o = range(1, L); out = {}
        nrm = sum(blk(W, 0, m) @ blk(W, m, 0)/E[m]**2 for m in o)
        out[0] = np.eye(d) - 0.5*nrm
        for n in o:
            t = -blk(W, n, 0)/E[n]
            for m in o:
                t = t + blk(W, n, m) @ blk(W, m, 0)/(E[n]*E[m])
                for l in o: t = t - blk(W, n, m) @ blk(W, m, l) @ blk(W, l, 0)/(E[n]*E[m]*E[l])
                a = blk(W, 0, m) @ blk(W, m, 0)
                x = blk(W, n, 0)
                p1 = (a @ x) if swap else (x @ a)
                t = t + p1/(E[n]**2*E[m]) + 0.5*p1/(E[n]*E[m]**2)
            out[n] = t
        res[swap] = out
    return res

print("unitarity of exact Omega:", np.abs(exact_omega(V, 0.01).conj().T @ exact_omega(V, 0.01) - np.eye(L*d)).max())
print("\nEq.(pt2), initial level i=2 (E_i != 0):  g, max error, error/g^3")
for g in (0.02, 0.01, 0.005):
    Om = exact_omega(V2, g); s = series2(V2, g, 2)
    err = max(np.abs(blk(Om, n, 2) - s[n]).max() for n in range(L))
    print(f"  {g:7.4f}  {err:.3e}  {err/g**3:.3f}")
print("\nEq.(pt3), vacuum (i=0):  g, max error, error/g^4 ; and with the norm factor on the LEFT: error/g^3")
for g in (0.02, 0.01, 0.005):
    Om = exact_omega(V, g); s = series3(V, g)
    err = max(np.abs(blk(Om, n, 0) - s[False][n]).max() for n in range(1, L))
    errw = max(np.abs(blk(Om, n, 0) - s[True][n]).max() for n in range(1, L))
    err0 = np.abs(blk(Om, 0, 0) - s[False][0]).max()
    print(f"  {g:7.4f}  n>=1: {err:.3e}  {err/g**4:.3f}   wrong order: {errw:.3e}  {errw/g**3:.3f}   n=0 (2nd order): {err0:.2e} {err0/g**3:.2f}")
