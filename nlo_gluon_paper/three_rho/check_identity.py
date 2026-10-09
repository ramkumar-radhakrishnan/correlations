"""Check Eqs. (plain), (Aanti), (antiA) with random matrices, including the one-rho remainder."""
import itertools
import numpy as np
rng = np.random.default_rng(5)
n = 6
A, B, C = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(3)]
cm = lambda X, Y: X @ Y - Y @ X
ac = lambda X, Y: X @ Y + Y @ X
ops = [A, B, C]
S = sum(ops[p[0]] @ ops[p[1]] @ ops[p[2]] for p in itertools.permutations(range(3))) / 6
two = 0.25 * (ac(cm(A, B), C) + ac(cm(A, C), B) + ac(cm(B, C), A))
R1 = (3 * cm(cm(A, B), C) + 2 * cm(A, cm(B, C)) + 2 * cm(B, cm(A, C)) + cm(cm(A, C), B) + cm(cm(B, C), A)) / 12
e1 = np.abs(A @ B @ C - (S + two + R1)).max()
# anticommutator forms: the remainders are double commutators; check that the difference is
# a combination of double commutators by verifying it vanishes for commuting-commutator algebras
lhs2 = A @ ac(B, C) - (2 * S + 0.5 * ac(cm(A, B), C) + 0.5 * ac(cm(A, C), B))
R2 = (cm(cm(A, B), C) + cm(cm(A, C), B)) / 6 + 0.5 * (cm(A, cm(B, C)) * 0)  # derived below
# exact one-rho remainder of A{B,C}: R1(A,B,C)+R1(A,C,B) + (1/4)([[B,C],A]+[[C,B],A]) -> compute directly
R1b = (3 * cm(cm(A, C), B) + 2 * cm(A, cm(C, B)) + 2 * cm(C, cm(A, B)) + cm(cm(A, B), C) + cm(cm(C, B), A)) / 12
e2 = np.abs(A @ ac(B, C) - (2 * S + 0.5 * ac(cm(A, B), C) + 0.5 * ac(cm(A, C), B) + R1 + R1b
                              + 0.25 * (ac(cm(B, C), A) + ac(cm(C, B), A)) - 0.25 * 0)).max()
# {B,C}A = A{B,C} - {[A,B],C} - {B,[A,C]}
e3 = np.abs(ac(B, C) @ A - (A @ ac(B, C) - ac(cm(A, B), C) - ac(B, cm(A, C)))).max()
print('ABC identity with explicit R1:', e1)
print('A{B,C} identity with explicit remainders:', e2)
print('{B,C}A = A{B,C} - [A,{B,C}]:', e3)
