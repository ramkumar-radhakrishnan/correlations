"""Is [N^(2), A^(1)] zero?  N^(2) ∝ ∫_{x,y} G(x-y) rho^a(x) rho^a(y),  A^(1)(p) ∝ rho^b(-p) (linear).
Test with an explicit projectile: two quarks at transverse positions x1 != x2,
rho^a(x) = t^a_(1) δ(x-x1) + t^a_(2) δ(x-x2)  acting on C^3 ⊗ C^3."""
import numpy as np
from lattice import su_n_generators
t = su_n_generators(3); I3 = np.eye(3)
T1 = [np.kron(ta, I3) for ta in t]          # charge of quark 1
T2 = [np.kron(I3, ta) for ta in t]          # charge of quark 2
G11, G22, G12 = 0.7, 0.7, 0.3               # any kernel values G(x1-x1), G(x2-x2), G(x1-x2), with G12 != 0
N2 = sum(G11*T1[a]@T1[a] + G22*T2[a]@T2[a] + 2*G12*T1[a]@T2[a] for a in range(8))
# A^(1)(p) ∝ rho^b(-p) = e^{-ip x1} t_1^b + e^{-ip x2} t_2^b ; take one Fourier mode
p_x1, p_x2 = 0.4, -1.1
for b in (0, 3, 7):
    A = np.exp(-1j*p_x1)*T1[b] + np.exp(-1j*p_x2)*T2[b]
    C = N2 @ A - A @ N2
    print(f'b={b}:  ||[N2, A]|| = {np.linalg.norm(C):.4f}   (||N2 A|| = {np.linalg.norm(N2@A):.4f})')
# where it does vanish: total charge, and a single point source
Qb = T1[0] + T2[0]
print('with the total charge Q^b = rho^b(p=0):  ||[N2, Q^b]|| =', f'{np.linalg.norm(N2@Qb - Qb@N2):.1e}')
N_single = 0.7*sum(T1[a]@T1[a] for a in range(8))
print('single point source (N ∝ Casimir):  ||[N, t_1^b]|| =', f'{np.linalg.norm(N_single@T1[0]-T1[0]@N_single):.1e}')
