# Referee 3, error item (Section 7.3): normalisation of the dipole emission probability and consistency with Eq. (7.2).
# (1) colour charges of the q qbar g colour singlet, exactly in charge space (explicit SU(N) matrices);
# (2) real - virtual with symbolic F_ab (sympy), finite N_c;
# (3) prefactors: emission probability from the paper's own Eq. (N2)/(PsiLO), LO prefactor from Eq. (LOsinglet).
# usage: python3 -I emission.py
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, sympy as sp
from colour import SUN


def qqg_charges(G):
    N, n = G.N, G.n; t = G.t; tb = -t.transpose(0, 2, 1); T = G.T
    I, Ia = np.eye(N), np.eye(n)
    Q1 = np.array([np.kron(np.kron(m, I), Ia) for m in t])
    Q2 = np.array([np.kron(np.kron(I, m), Ia) for m in tb])
    Qg = np.array([np.kron(np.kron(I, I), m) for m in T])
    # q qbar g singlet: psi_{i j a} = (t^a)_{i j} (quark index i, antiquark index j, gluon a), normalised
    psi = np.einsum('aij->ija', t).reshape(N * N * n); psi = psi / np.linalg.norm(psi)
    tot = Q1 + Q2 + Qg
    sing = np.abs(np.einsum('aij,j->ai', tot, psi)).max()
    if sing > 1e-12:                       # other index convention for the gluon
        psi = np.einsum('aij->jia', t).reshape(N * N * n); psi /= np.linalg.norm(psi)
        sing = np.abs(np.einsum('aij,j->ai', tot, psi)).max()
    Qs = {'1': Q1, '2': Q2, 'g': Qg}
    dot = {(a, b): np.real(np.einsum('i,eij,ejk,k->', psi.conj(), Qs[a], Qs[b], psi)) for a in Qs for b in Qs}
    # q qbar singlet
    phi = np.eye(N).reshape(N * N) / np.sqrt(N)
    q1 = np.array([np.kron(m, I) for m in t]); q2 = np.array([np.kron(I, m) for m in tb])
    dot0 = {('1', '1'): np.real(np.einsum('i,eij,ejk,k->', phi, q1, q1, phi)),
            ('2', '2'): np.real(np.einsum('i,eij,ejk,k->', phi, q2, q2, phi)),
            ('1', '2'): np.real(np.einsum('i,eij,ejk,k->', phi, q1, q2, phi))}
    return sing, dot, dot0


def main():
    Nc = sp.Symbol('N_c', positive=True)
    for N in (2, 3, 4):
        G = SUN(N); sing, dot, dot0 = qqg_charges(G)
        CF = (N * N - 1) / (2 * N)
        pred = {('1', '1'): CF, ('2', '2'): CF, ('g', 'g'): N, ('1', '2'): 1 / (2 * N), ('1', 'g'): -N / 2, ('2', 'g'): -N / 2}
        dev = max(abs(dot[k] - v) for k, v in pred.items())
        print('SU(%d): q qbar g singlet |Q_tot psi| = %.1e;  Q1.Q1 = %.6f  Q2.Q2 = %.6f  Qg.Qg = %.6f  Q1.Q2 = %.6f  Q1.Qg = %.6f  Q2.Qg = %.6f'
              '  (dev from C_F, C_F, N_c, 1/(2N_c), -N_c/2, -N_c/2: %.1e);  q qbar singlet: Q1.Q2 = %.6f (-C_F = %.6f)' % (
                  N, sing, dot[('1', '1')], dot[('2', '2')], dot[('g', 'g')], dot[('1', '2')], dot[('1', 'g')], dot[('2', 'g')],
                  dev, dot0[('1', '2')], -CF))
    # ---- (2) real - virtual, symbolic ----
    F = {(a, b): sp.Symbol('F_%s%s' % (a, b)) for a in '12z' for b in '12z'}
    CFs = (Nc**2 - 1) / (2 * Nc)
    Qd = {('1', '1'): CFs, ('2', '2'): CFs, ('z', 'z'): Nc, ('1', '2'): 1 / (2 * Nc), ('2', '1'): 1 / (2 * Nc),
          ('1', 'z'): -Nc / 2, ('z', '1'): -Nc / 2, ('2', 'z'): -Nc / 2, ('z', '2'): -Nc / 2}
    real = sum(Qd[k] * F[k] for k in Qd)
    virt = CFs * (F[('1', '1')] + F[('2', '2')] - F[('1', '2')] - F[('2', '1')])
    D = lambda a, b: (F[(a, a)] + F[(b, b)] - F[(a, b)] - F[(b, a)]) / (2 * Nc)
    bk = D('1', 'z') + D('z', '2') - D('1', '2')
    print('real - virtual - N_c^2 [D_1z + D_z2 - D_12] =', sp.simplify(sp.expand(real - virt - Nc**2 * bk)))
    print('LO colour sum for the dipole: C_F[F11+F22-F12-F21] - 2 N_c C_F D_12 =', sp.simplify(virt - 2 * Nc * CFs * D('1', '2')),
          ';  2 N_c C_F = N_c^2 - 1:', sp.simplify(2 * Nc * CFs - (Nc**2 - 1)))
    # ---- (3) prefactors ----
    g, als, pi, kp, L, KD = sp.symbols('g alpha_s pi k_plus L K_D', positive=True)
    gsq = 4 * sp.pi * als
    # emission probability per d^2z: from Eq. (N2) [N_2 = -(1/2) x probability] or |Psi_LO|^2 with <g|g> = delta:
    #   P = g^2/(4 pi^3) ln(vee/Lambda) int_z K(x,y,z) <rho^a(x) rho^a(y)>, and for the dipole int K n = C_F K_D
    P = gsq / (4 * sp.pi**3) * L * CFs * KD
    print('emission probability g^2 C_F/(4 pi^3) K_D L =', sp.simplify(P), ' (= alpha_s C_F/pi^2 K_D L)')
    # LO spectrum (LOsinglet): alpha_s/(4 pi^4 (N^2-1) k+) sum <Q_s Q_s'> F(x_s, x_s')  -> dipole: alpha_s/(4 pi^4 k+) D_12
    PLO = als / (4 * sp.pi**4 * (Nc**2 - 1) * kp)
    lo_dip = sp.simplify(PLO * 2 * Nc * CFs)
    print('LO for a dipole = (', lo_dip, ') x D_12   -> "LO prefactor" of (7.2) = alpha_s/(4 pi^4 k+)')
    delta = sp.simplify(P / KD * PLO * Nc**2)                     # coefficient of K_D [D1z+Dz2-D12]
    target = als * Nc / (2 * sp.pi**2) * L * lo_dip
    print('delta coefficient / [alpha_s N_c/(2 pi^2) L x LO prefactor] =', sp.simplify(delta / target))
    Pref = gsq / (4 * sp.pi**4) * L * CFs * KD                    # the form quoted by the referee
    print('with the referee-quoted g^2 C_F/(4 pi^4): ratio =', sp.simplify(Pref / KD * PLO * Nc**2 / target))
    # ---- Parseval: WW density in coordinate space from the Gunion-Bertsch spectrum ----
    # FT convention of the paper: f(k) = int d^2x e^{ikx} f(x); int d^2k/(2pi)^2 e^{-ikz} k^i/k^2 = -i z^i/(2 pi z^2)
    # => int d^2k |k^i/k^2|^2 = (2pi)^2 int d^2z |z^i/(2pi z^2)|^2 = int d^2z/z^2 : dN/(dY d^2k) = alpha_s C_F/(pi^2 k^2)
    #    <-> dN/(dY d^2z) = alpha_s C_F/(pi^2 z^2), same normalisation alpha_s C_F/pi^2.
    print('int_0^inf dk J_1(k z) = 1/z (standard)  -> int d^2k/(2pi)^2 e^{-ikz} k^i/k^2 = -i z^i/(2 pi z^2): Parseval gives the same alpha_s C_F/pi^2')

if __name__ == '__main__':
    main()
