# Referee 3, Major 7 (a), (b): the cubic finite terms for a single quark and for a colour-singlet dipole,
# computed EXACTLY in colour (charge) space, and compared with the split
#   <C rho rho rho> = [symmetrised part  N_c/((N_c^2-4)(N_c^2-1)) Delta o]  +  [contact terms of eq. cubicavg].
# usage: python3 -I quark_dipole.py
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colour import SUN, structures, ORDER, POINTS, tdelta, delta_table, Ginv

np.set_printoptions(precision=6, suppress=True)


class State:
    """Charges Q_s^a on sites s (matrices on the colour Hilbert space) and a density matrix."""
    def __init__(s, G, kind):
        N = G.N; I = np.eye(N); t = G.t; tb = -t.transpose(0, 2, 1)      # antiquark generators -t^T
        if kind == 'quark':
            s.Q = [t]; s.rho = I / N; s.sign = [1]
        elif kind == 'antiquark':
            s.Q = [tb]; s.rho = I / N; s.sign = [-1]
        elif kind in ('dipole', 'dipole_swapped'):
            q = np.array([np.kron(m, I) for m in t]); qb = np.array([np.kron(I, m) for m in tb])
            v = np.eye(N).reshape(N * N) / np.sqrt(N)                     # singlet sum_i |i i>/sqrt(N)
            s.rho = np.outer(v, v.conj())
            s.Q = [q, qb] if kind == 'dipole' else [qb, q]                # swapped: antiquark at site 0
            s.sign = [1, -1] if kind == 'dipole' else [-1, 1]
            Qtot = q + qb
            assert np.abs(np.einsum('aij,j->ai', Qtot, v)).max() < 1e-14  # colour singlet
        s.P = len(s.Q); s.G = G
        P, n = s.P, G.n
        s.E3 = np.einsum('paij,qbjk,rckl,li->paqbrc', s.Q, s.Q, s.Q, s.rho, optimize=True)   # <Q_p^a Q_q^b Q_r^c>
        s.E2 = np.einsum('paij,qbjk,ki->paqb', s.Q, s.Q, s.rho, optimize=True)
        s.n2 = np.einsum('paqa->pq', s.E2).real
        # explicit symmetrisation: {A B C}_sym with A=(p,a), B=(q,b), C=(r,c)
        E = s.E3
        S3 = (E + E.transpose(0, 1, 4, 5, 2, 3) + E.transpose(2, 3, 0, 1, 4, 5) + E.transpose(2, 3, 4, 5, 0, 1)
              + E.transpose(4, 5, 0, 1, 2, 3) + E.transpose(4, 5, 2, 3, 0, 1)) / 6
        s.S3 = S3
        s.o = np.einsum('paqbrc,abc->pqr', S3, G.d).real
        s.gf = np.einsum('paqbrc,abc->pqr', S3, G.f)


def exact(st, C, typ, sites):
    s1, s2, s3 = sites; E = st.E3
    if typ == 'ord':
        return np.einsum('abc,abc->', C, E[s1, :, s2, :, s3, :])
    if typ == 'A23':
        return np.einsum('abc,abc->', C, E[s1, :, s2, :, s3, :] + E[s1, :, s3, :, s2, :].transpose(0, 2, 1))
    if typ == 'A12':
        return np.einsum('abc,abc->', C, E[s1, :, s2, :, s3, :] + E[s2, :, s1, :, s3, :].transpose(1, 0, 2))


def split(st, C, typ, sites):
    """(symmetrised part, contact part) of eq. cubicavg, with S^{abc} = d-form of eq. A2d and n of the state."""
    G = st.G; N = G.N; s1, s2, s3 = sites
    sym = tdelta(typ) * N / ((N * N - 4) * (N * N - 1)) * np.einsum('abc,abc->', C, G.d) * st.o[s1, s2, s3] if N > 2 else 0.0
    Phi = np.einsum('abc,abc->', C, G.f); n = st.n2; dl = lambda i, j: 1.0 if i == j else 0.0
    if typ == 'ord':
        con = 1j * Phi / (2 * (N * N - 1)) * (dl(s1, s2) * n[s1, s3] + dl(s2, s3) * n[s1, s2] - dl(s1, s3) * n[s1, s2])
    elif typ == 'A23':
        con = 1j * Phi / (N * N - 1) * (dl(s1, s2) * n[s1, s3] - dl(s1, s3) * n[s1, s2])
    else:
        con = 1j * Phi / (N * N - 1) * (dl(s2, s3) * n[s1, s2] - dl(s1, s3) * n[s1, s2])
    return sym, con


def scale(vals, C):
    """relative scale; values that vanish identically (e.g. I.II N2b at a triple coincidence) are measured against max|C|"""
    m = max(vals); c = np.abs(C).max()
    return m if m > 1e-12 * c else c


def Ufor(G, rng, rho_pts, sites, site_U, eps=None):
    U = {p: G.adj(G.randV(rng, eps)) for p in POINTS}
    for p, s in zip(rho_pts, sites): U[p] = site_U[s]
    return U


def main():
    worst = {}
    for N in (2, 3, 4):
        G = SUN(N); rng = np.random.default_rng(100 + N)
        print('=' * 100); print('SU(%d):  C_F = %.6f' % (N, G.CF))
        # ---------- the states: n and o ----------
        for kind in ('quark', 'antiquark', 'dipole', 'dipole_swapped'):
            st = State(G, kind)
            oq = (N * N - 4) * (N * N - 1) / (4 * N * N)
            sg = np.array(st.sign)
            o_pred = oq * np.einsum('p,q,r->pqr', sg, sg, sg)
            n_pred = G.CF * np.outer(sg, sg)
            print('  %-15s n = %s (pred C_F s_i s_j: dev %.1e)   o = %s (pred o_q s_i s_j s_k, o_q=(N^2-4)(N^2-1)/(4N^2)=%.6f: dev %.1e)'
                  '   max|g_f-part| %.1e   max|S3 - d-form| %.1e' % (
                      kind, st.n2.ravel(), np.abs(st.n2 - n_pred).max(), st.o.ravel(), oq, np.abs(st.o - o_pred).max(),
                      np.abs(st.gf).max(),
                      np.abs(st.S3 - (np.einsum('abc,pqr->paqbrc', G.d, st.o) * N / ((N * N - 4) * (N * N - 1)) if N > 2 else 0)).max()))
        # single-quark identities
        t = G.t
        dttt = np.einsum('abc,aij,bjk,cki->', G.d, t, t, t) / N
        fttt = np.einsum('abc,aij,bjk,cki->', G.f, t, t, t) / N
        print('  quark: <d^{abc} t^a t^b t^c> = %.6f   C_F (N^2-4)/(2N) = %.6f ;  <f^{abc} t^a t^b t^c> = %.6fi   i N C_F/2 = %.6fi' % (
            dttt.real, G.CF * (N * N - 4) / (2 * N), fttt.imag, N * G.CF / 2))
        tr3 = np.einsum('aij,bjk,cki->abc', t, t, t)
        print('  quark: max| tr(t^a t^b t^c) - (d^{abc} + i f^{abc})/4 | = %.1e' % np.abs(tr3 - (G.d + 1j * G.f) / 4).max())
        # ---------- all 13 structures: exact vs split ----------
        for kind in ('quark', 'antiquark', 'dipole', 'dipole_swapped'):
            st = State(G, kind); w = 0.0
            site_U = [G.adj(G.randV(rng)) for _ in range(st.P)]
            print('  -- %s: exact <C rho rho rho> vs symmetrised + contact, all structures, all %d assignments' % (kind, st.P ** 3))
            for name in ORDER:
                rows = []
                for sites in itertools.product(range(st.P), repeat=3):
                    for rep in range(3):
                        # rho points get the Wilson line of their site; other points independent
                        rho_pts = structures(G, {p: np.eye(G.n) for p in POINTS})[name][2]
                        U = Ufor(G, rng, rho_pts, sites, site_U)
                        C, typ, pts = structures(G, U)[name]
                        ex = exact(st, C, typ, sites); sy, co = split(st, C, typ, sites)
                        dev = abs(ex - sy - co) / scale([abs(ex), abs(sy) + abs(co)], C)
                        w = max(w, dev)
                        if rep == 0: rows.append((sites, ex, sy, co))
                if kind == 'dipole' and N == 3:
                    print('     %-14s' % name + '  '.join('%s:[ex %+.4f%+.4fi | sym %+.4f | con %+.4fi]' % (
                        ''.join('qQ'[s] for s in r[0]), r[1].real, r[1].imag, r[2], r[3].imag) for r in rows[:4]))
            print('     largest relative deviation exact vs (sym + contact): %.1e' % w)
            worst[(N, kind)] = w
        # ---------- quark vs antiquark: C-odd = sym, C-even = contact ----------
        stq, sta = State(G, 'quark'), State(G, 'antiquark'); wc = 0.0
        std, sts = State(G, 'dipole'), State(G, 'dipole_swapped')
        site_U = [G.adj(G.randV(rng)) for _ in range(2)]
        for name in ORDER:
            for sites in itertools.product(range(2), repeat=3):
                rho_pts = structures(G, {p: np.eye(G.n) for p in POINTS})[name][2]
                U = Ufor(G, rng, rho_pts, sites, site_U)
                C, typ, pts = structures(G, U)[name]
                # dipole and its charge conjugate (q <-> qbar at the same positions)
                ed, es = exact(std, C, typ, sites), exact(sts, C, typ, sites)
                sy, co = split(std, C, typ, sites)
                sc = scale([abs(ed) + abs(es), abs(sy) + abs(co)], C)
                wc = max(wc, abs((ed - es) / 2 - sy) / sc, abs((ed + es) / 2 - co) / sc)
                if sites == (0, 0, 0):
                    Us = {p: (site_U[0] if p in rho_pts else U[p]) for p in POINTS}
                    C1 = structures(G, Us)[name][0]
                    eq_, ea = exact(stq, C1, typ, (0, 0, 0)), exact(sta, C1, typ, (0, 0, 0))
                    sq, cq = split(stq, C1, typ, (0, 0, 0))
                    sc = scale([abs(eq_) + abs(ea), abs(sq) + abs(cq)], C1)
                    wc = max(wc, abs((eq_ - ea) / 2 - sq) / sc, abs((eq_ + ea) / 2 - cq) / sc)
        print('  C-parity: (P - P^C)/2 = symmetrised part, (P + P^C)/2 = contact part (quark/antiquark and dipole/swapped dipole): largest rel. dev %.1e' % wc)
        worst[(N, 'Cparity')] = wc
        # ---------- ordering: six orderings of the same operators ----------
        if N == 3:
            st = State(G, 'dipole'); site_U = [G.adj(G.randV(rng)) for _ in range(2)]
            name = 'I.III(z)'; rho_pts = ('xp', 'x', 'y')
            print('  -- ordering test (I.III(z), dipole, assignment qqQ): the six orderings of the same three operators C^{abc} rho^a(p1) rho^b(p2) rho^c(p3)')
            sites = (0, 0, 1)
            U = Ufor(G, rng, rho_pts, sites, site_U); C = structures(G, U)[name][0]
            wo = 0.0; syms = []
            for perm in itertools.permutations(range(3)):
                # operator product rho^{a_perm0}(p_perm0) rho^{..} rho^{..}: an 'ord' product with permuted tensor/sites
                Cp = np.transpose(C, perm); sp = tuple(sites[i] for i in perm)
                ex = exact(st, Cp, 'ord', sp); sy, co = split(st, Cp, 'ord', sp)
                syms.append(sy); wo = max(wo, abs(ex - sy - co) / abs(ex))
                print('     ordering %s: exact %+.6f%+.6fi = sym %+.6f + contact %+.6fi' % (perm, ex.real, ex.imag, sy, co.imag))
            print('     symmetrised part identical for all orderings: spread %.1e; exact = sym + contact for each ordering: %.1e' % (
                np.ptp(syms), wo))
            worst[(N, 'ordering')] = wo
    print('=' * 100)
    print('SUMMARY (largest relative deviations):', {k: float('%.2e' % v) for k, v in worst.items()})


if __name__ == '__main__':
    main()
