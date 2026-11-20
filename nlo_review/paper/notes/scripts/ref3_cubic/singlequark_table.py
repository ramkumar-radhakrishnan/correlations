# Single quark at a (all three charge points at a): closed forms of the cubic colour weights
#   <C rho rho rho> = [Delta_loc + i Phi_loc * theta_ord] / (4 N_c),
# Delta_loc = sum_orderings C d (at the triple coincidence), Phi_loc = C f (ordered structures only; zero contact for
# the anticommutator types since tr(t^a {t^b, t^c}) = d^{abc}/2).  Checked against the exact trace (1/N) tr(...) .
# usage: python3 -I singlequark_table.py
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colour import SUN, structures, ORDER, POINTS, tdelta, Dinv, Finv, Ginv, Qdinv
from quark_dipole import State, exact


def formulas(G, U):
    N = G.N; a = 'x'
    D = lambda p, q: Dinv(U, p, q); F = lambda p, q, r: Finv(G, U, p, q, r)
    g = lambda p, q, r: Ginv(G, U, p, q, r); qd = lambda p, q, r, s: Qdinv(G, U, p, q, r, s)
    C3 = N * (N * N - 1)
    Dl = {'I.II.N2b': 0.0, 'I.III(z)': g(a, 'z', 'wp'), 'I.III(w)': g(a, 'w', 'wp'),
          'I.V.N1': -2 * (g('xx', a, 'z') + g('xx', 'z', 'wp')),
          'II.I.N1a': 2 * g(a, 'z', 'wp') - 2 * qd('z', a, a, 'w'), 'II.I.N1b': 2 * g(a, 'z', 'w') - 2 * qd('z', a, a, 'wp'),
          'III.II': qd('z', a, a, 'w'), "III.III.1(y')": 2 * g('yp', 'z', 'w') - 2 * g(a, 'z', 'w'), 'III.III.1(z)': -2 * g(a, 'z', 'w'),
          "III.V(y')": g(a, 'yp', 'w'), 'III.V(z)': g(a, 'z', 'w'), "III.VI(y')": -g(a, 'yp', 'z'), 'III.VI(z)': 0.0}
    Pl = {'I.III(z)': C3 - N * D(a, 'wp') - N * D(a, 'z') + F(a, 'z', 'wp'),
          'I.III(w)': C3 - N * D(a, 'wp') - N * D(a, 'w') + F(a, 'w', 'wp'),
          'III.II': C3 + N * D(a, 'w') - N * D(a, 'z') - F(a, 'z', 'w'),
          "III.V(y')": -C3 + N * D(a, 'yp') + N * D(a, 'w') - F(a, 'yp', 'w'),
          'III.V(z)': -C3 + N * D(a, 'z') + N * D(a, 'w') - F(a, 'z', 'w'),
          "III.VI(y')": C3 - N * D(a, 'yp') - N * D(a, 'z') + F(a, 'yp', 'z'),
          'III.VI(z)': C3 - N * D(a, 'z')}
    return Dl, Pl


def rel(a, b, c):
    m = max(abs(a), abs(b))
    return abs(a - b) / (m if m > 1e-10 * c else c)   # identically vanishing values: deviation relative to max|C|


def main():
    for N in (3, 4):
        G = SUN(N); st = State(G, 'quark'); rng = np.random.default_rng(31 + N); w = 0.0
        for _ in range(4):
            Vs = {p: G.adj(G.randV(rng)) for p in POINTS + ['xx']}
            # the 'x' of Row I.V is not a charge point: it is renamed 'xx'; its charge points x', y, y' sit at a = 'x'
            for name in ORDER:
                C, typ, pts = structures(G, {p: np.eye(G.n) for p in POINTS})[name]
                U = dict(Vs)
                for p in pts: U[p] = Vs['x']
                if name == 'I.V.N1': U['x'] = Vs['xx']
                C = structures(G, U)[name][0]
                ex = exact(st, C, typ, (0, 0, 0))
                Dl, Pl = formulas(G, {**U, 'x': Vs['x'], 'xx': Vs['xx']})
                pred = (Dl[name] + 1j * Pl.get(name, 0.0) * (typ == 'ord')) / (4 * N)
                dl_dir = tdelta(typ) * np.einsum('abc,abc->', C, G.d); pl_dir = np.einsum('abc,abc->', C, G.f)
                c = np.abs(C).max()
                w = max(w, rel(ex, pred, c), rel(dl_dir, Dl[name], c))
                if typ == 'ord': w = max(w, rel(pl_dir, Pl[name], c))
                else: w = max(w, rel(ex.imag, 0.0, c) if abs(ex) < 1e-10 * c else abs(ex.imag) / abs(ex))
        print('SU(%d): single quark, 13 structures x 4 targets: exact (1/N) tr = [Delta_loc + i Phi_loc theta_ord]/(4N) with the'
              ' closed forms: max rel. dev %.1e' % (N, w))


if __name__ == '__main__':
    main()
