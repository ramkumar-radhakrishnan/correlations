# Referee 3, Section 4.3 and Major 7: charge conjugation.
#  C: rho^a t^a -> -(rho^a t^a)^T, i.e. rho^a -> -eta_a rho^a ;  V -> V^*, i.e. U^{ab} -> eta_a eta_b U^{ab}.
# (1) C is an automorphism of the charge algebra; n is C-even, o is C-odd; f-type tensors pick up -1, d-type +1.
# (2) the invariant three-index colour tensors are d and f only (null space of the adjoint action), and the only
#     C-even coefficient tensor of rho rho rho is f: a target-independent (c-number) coefficient cannot produce o.
# (3) target invariants: D, F, Q C-even; G, Q_d C-odd.
# (4) exact states: dipole and its charge conjugate (q <-> qbar): symmetrised part flips, contact part does not;
#     (dipole -> conjugate dipole, V -> V^*) together leave every cubic term invariant (C-invariance of the result);
#     one-step BFKL: the q qbar g charge products are symmetric under q <-> qbar, so delta n is C-even.
# usage: python3 -I cparity.py
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colour import SUN, structures, ORDER, POINTS, Dinv, Finv, Ginv, Qinv, Qdinv
from quark_dipole import State, exact, split, Ufor


def main():
    for N in (3, 4):
        G = SUN(N); eta = G.eta; rng = np.random.default_rng(7 + N)
        e3 = np.einsum('a,b,c->abc', eta, eta, eta)
        print('SU(%d): max|eta eta eta f + f| = %.1e   max|eta eta eta d - d| = %.1e' % (
            N, np.abs(e3 * G.f + G.f).max(), np.abs(e3 * G.d - G.d).max()))
        # (1) automorphism: rho' = -eta rho satisfies [rho'^a, rho'^b] = i f^{abc} rho'^c  (on the fundamental: t' = -t^T)
        tp = -G.t.transpose(0, 2, 1)
        comm = np.einsum('aij,bjk->abik', tp, tp) - np.einsum('bij,ajk->abik', tp, tp)
        print('   [t^C_a, t^C_b] - i f t^C_c : %.1e ;  t^C_a = -eta_a t_a : %.1e' % (
            np.abs(comm - 1j * np.einsum('abc,cij->abij', G.f, tp)).max(), np.abs(tp + eta[:, None, None] * G.t).max()))
        # (2) invariant 3-tensors (null space computed for SU(3) only: the SVD for SU(4) is 3375 x 50625)
        n = G.n; T = G.T; I = np.eye(n)
        M = []
        if N > 3:
            M = None
        for e in (range(n) if M is not None else []):
            Ae = (np.einsum('ij,kl,mn->ikmjln', T[e], I, I) + np.einsum('ij,kl,mn->ikmjln', I, T[e], I)
                  + np.einsum('ij,kl,mn->ikmjln', I, I, T[e])).reshape(n**3, n**3)
            M.append(Ae)
        if M is not None:
            M = np.vstack(M)
            sv = np.linalg.svd(M, compute_uv=False)
            null = int((sv < 1e-10).sum())
            print('   dimension of the space of invariant tensors X^{abc}: %d (d and f)' % null)
        # C acts on a coefficient tensor c_{abc} of rho^a rho^b rho^c as c -> -eta_a eta_b eta_c c
        print('   C on coefficient tensors: f -> %+.0f f, d -> %+.0f d   => only f is C-even; f^{abc} d^{abc} = %.1e (no o)' % (
            np.sum(-e3 * G.f * G.f) / np.sum(G.f * G.f), np.sum(-e3 * G.d * G.d) / np.sum(G.d * G.d), np.abs(np.sum(G.f * G.d))))
        # (3) target invariants
        Vs = {p: G.randV(rng) for p in POINTS}
        U = {p: G.adj(v) for p, v in Vs.items()}; Uc = {p: G.adj(v.conj()) for p, v in Vs.items()}
        print('   U(V*) = eta U eta : %.1e' % max(np.abs(Uc[p] - eta[:, None] * U[p] * eta[None, :]).max() for p in POINTS))
        vals = [('D(x,y)', Dinv(U, 'x', 'y'), Dinv(Uc, 'x', 'y')), ('F(x,y,z)', Finv(G, U, 'x', 'y', 'z'), Finv(G, Uc, 'x', 'y', 'z')),
                ('Q(x,y;z,w)', Qinv(G, U, 'x', 'y', 'z', 'w'), Qinv(G, Uc, 'x', 'y', 'z', 'w')),
                ('G(x,y,z)', Ginv(G, U, 'x', 'y', 'z'), Ginv(G, Uc, 'x', 'y', 'z')),
                ('Q_d(x,y;z,w)', Qdinv(G, U, 'x', 'y', 'z', 'w'), Qdinv(G, Uc, 'x', 'y', 'z', 'w'))]
        print('   ' + '  '.join('%s: %+.5f -> %+.5f' % (a, b.real, c.real) for a, b, c in vals))
        # (4) exact dipole states
        sd, sc = State(G, 'dipole'), State(G, 'dipole_swapped')
        print('   dipole o(q,q,q) = %+.6f, conjugate dipole o(q,q,q) = %+.6f ; n identical: %.1e' % (
            sd.o[0, 0, 0], sc.o[0, 0, 0], np.abs(sd.n2 - sc.n2).max()))
        wsymP = wsymT = wcon = winv = 0.0
        Vsite = [G.randV(rng) for _ in range(2)]
        for name in ORDER:
            rho_pts = structures(G, {p: np.eye(n) for p in POINTS})[name][2]
            for sites in itertools.product(range(2), repeat=3):
                Vp = {p: G.randV(rng) for p in POINTS}
                for p, s in zip(rho_pts, sites): Vp[p] = Vsite[s]
                U = {p: G.adj(v) for p, v in Vp.items()}; Uc = {p: G.adj(v.conj()) for p, v in Vp.items()}
                C, typ, _ = structures(G, U)[name]; Cc = structures(G, Uc)[name][0]
                s1, c1 = split(sd, C, typ, sites); s2, c2 = split(sc, C, typ, sites); s3, c3 = split(sd, Cc, typ, sites)
                e_ = exact(sd, C, typ, sites); e_C = exact(sc, Cc, typ, sites)
                sc_ = abs(s1) + abs(c1) + 1e-300
                wsymP = max(wsymP, abs(s1 + s2) / sc_); wsymT = max(wsymT, abs(s1 + s3) / sc_)
                wcon = max(wcon, abs(c1 - c2) / sc_, abs(c1 - c3) / sc_)
                winv = max(winv, abs(e_ - e_C) / (abs(e_) + 1e-300))
        print('   all 13 structures, 8 assignments: sym part odd under projectile C (%.1e) and under target C (%.1e);'
              ' contact part even under each (%.1e); exact <C rho rho rho> invariant under joint C: %.1e' % (wsymP, wsymT, wcon, winv))


if __name__ == '__main__':
    main()
