"""Check: X1X2X3X4 = S(X1..X4) + 1/2 sum_{i<j} S([Xi,Xj], Xk, Xl) + (terms with two or more commutators)."""
import itertools
import numpy as np
rng = np.random.default_rng(1)
n = 9
X = [rng.normal(size=(n, n)) for _ in range(4)]
cm = lambda A, B: A @ B - B @ A
def S(ops):
    ps = list(itertools.permutations(range(len(ops))))
    out = 0
    for p in ps:
        m = np.eye(n)
        for k in p:
            m = m @ ops[k]
        out = out + m
    return out / len(ps)
lhs = X[0] @ X[1] @ X[2] @ X[3]
first = S(X)
second = 0
for i, j in itertools.combinations(range(4), 2):
    rest = [X[k] for k in range(4) if k not in (i, j)]
    second = second + 0.5 * S([cm(X[i], X[j])] + rest)
D = lhs - first - second
# span of products with two Lie factors or one: S([[Xi,Xj],Xk], Xl), S([Xi,Xj],[Xk,Xl]), nested triple commutators
span = []
idx = range(4)
for i, j, k in itertools.permutations(idx, 3):
    l = [m for m in idx if m not in (i, j, k)][0]
    span.append(S([cm(cm(X[i], X[j]), X[k]), X[l]]))
    span.append(cm(cm(cm(X[i], X[j]), X[k]), X[l]))
for i, j in itertools.combinations(idx, 2):
    k, l = [m for m in idx if m not in (i, j)]
    span.append(S([cm(X[i], X[j]), cm(X[k], X[l])]))
A = np.array([s.ravel() for s in span]).T
coef, *_ = np.linalg.lstsq(A, D.ravel(), rcond=None)
res = np.abs(A @ coef - D.ravel()).max()
# control: a wrong coefficient (1/3 instead of 1/2) must fail
Dw = lhs - first - (2 / 3) * second
coefw, *_ = np.linalg.lstsq(A, Dw.ravel(), rcond=None)
print('residual with 1/2:', res, '   with 1/3 (control):', np.abs(A @ coefw - Dw.ravel()).max(), '  |D| =', np.abs(D).max())
