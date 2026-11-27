"""Symbolic first-order perturbation theory for Chapter 4 (quartic potential and symmetric well).

Units hbar = m = 1.  Harmonic basis |m> with y = (a + a^dagger)/sqrt(2).
Quartic:        H = p^2/2 + x^2/2 + b x^4                 (omega = 1, y = x)
Symmetric well: H = p^2/2 + 1/3 + lam^2 x^2/45 + 2 lam^4 x^4/945 + lam^6 x^6/4725
                omega = lam*sqrt(2/45), y = sqrt(omega) x,
                H'/omega = beta4 y^4 + beta6 y^6, beta4 = sqrt(10) lam/14, beta6 = 3 lam^2/28.
Position:  psi_n(x) = N [ |n> + sum_m c_m |m> ],  <y|m> = H_m(y) e^{-y^2/2}/sqrt(2^m m! sqrt(pi))
Momentum:  <q|m> = (-i)^m H_m(q) e^{-q^2/2}/sqrt(2^m m! sqrt(pi)),  q = k/sqrt(omega)
"""
import sympy as sp

D = 16
a = sp.zeros(D, D)
for m in range(1, D):
    a[m - 1, m] = sp.sqrt(m)
Y = (a + a.T) / sp.sqrt(2)
Y4 = sp.simplify(Y ** 4)
Y6 = sp.simplify(Y ** 6)
b, lam, y, q = sp.symbols('b lambda y q', positive=True)

def herm_fn(m, var):
    return sp.hermite(m, var) / sp.sqrt(2 ** m * sp.factorial(m))

def first_order(n, Hp):
    """coefficients c_m (m != n) of the first-order correction, Hp in units of omega."""
    c = {}
    for m in range(D - 6):
        if m == n:
            continue
        el = sp.simplify(Hp[m, n])
        if el != 0:
            c[m] = sp.simplify(el / (n - m))
    return c

def report(name, n, c, var_sym):
    print('=' * 70)
    print(name, 'n =', n)
    print('  kets:  |%d>' % n, ' + '.join('(%s)|%d>' % (sp.nsimplify(v), m) for m, v in sorted(c.items())))
    norm2 = sp.simplify(1 + sum(v ** 2 for v in c.values()))
    print('  1/N^2 =', sp.expand(norm2))
    P = sp.expand(sp.simplify(herm_fn(n, var_sym) + sum(v * herm_fn(m, var_sym) for m, v in c.items())))
    print('  position polynomial P_n(y) [psi = N pi^-1/4 e^{-y^2/2} P_n]:')
    print('    ', sp.collect(P, var_sym))
    # momentum: (-i)^m relative to (-i)^n: factor (-1)^((m-n)/2)
    Q = sp.expand(sp.simplify(herm_fn(n, var_sym) + sum(v * sp.Integer(-1) ** ((m - n) // 2) * herm_fn(m, var_sym)
                                                      for m, v in c.items())))
    print('  momentum polynomial Q_n(q) [phi = (-i)^n N pi^-1/4 e^{-q^2/2} Q_n]:')
    print('    ', sp.collect(Q, var_sym))
    return P, Q, norm2

results = {}
for n in range(4):
    c = first_order(n, b * Y4)
    P, Q, N2 = report('QUARTIC', n, c, y)
    E1 = sp.simplify(b * Y4[n, n])
    print('  E_n^(1) =', E1)
    results[('q', n)] = (c, P, Q, N2, E1)

beta4 = sp.sqrt(10) * lam / 14
beta6 = sp.Rational(3, 28) * lam ** 2
for n in range(4):
    c = first_order(n, beta4 * Y4 + beta6 * Y6)
    P, Q, N2 = report('SYMMETRIC WELL', n, c, y)
    E1 = sp.simplify(beta4 * Y4[n, n] + beta6 * Y6[n, n])
    print('  E_n^(1)/omega =', sp.expand(E1))
    results[('s', n)] = (c, P, Q, N2, E1)

import pickle
with open('derive_results.pkl', 'wb') as f:
    pickle.dump({k: tuple(sp.srepr(x) if not isinstance(x, dict) else {kk: sp.srepr(vv) for kk, vv in x.items()}
                          for x in v) for k, v in results.items()}, f)
