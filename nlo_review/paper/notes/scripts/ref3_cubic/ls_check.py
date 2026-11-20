# Referee 3, Major 7 (c): comparison with Li-Skokov (arXiv:2104.01879, JHEP 06 (2021) 141) and McLerran-Skokov
# (arXiv:1611.09870).  Read from the arXiv sources (export.arxiv.org); the formulas below are transcribed from
# LS eqs. (NLO_single_gluon_production), (coefficientfunc_E1), (effective_vertex_g3_second) [labels of the source].
#  (1) transcription check: E_2 (printed, five lines) = Gamma_2^j L^j(p3,-k) (printed factorised form);
#  (2) k-parity of LS's n^(4): E_1 real and even under reflection of all momenta => its contribution is EVEN in k;
#      E_2 = Re + i sgn(k x p)(...): Re part even, Im part odd (LS's own statement);
#  (3) d-projection (A2'): 2-Wilson-line term f^{ade} rho^d rho^e rho^{e2} (U U^dag)^{a e2} -> f^{ade} d^{d e e2} = 0;
#      3-Wilson-line term -> G-type invariant (totally antisymmetric under simultaneous exchange of (q_i, p_i));
#      the only exchange that is a symmetry of the integration domain is (p -> k - p, p1 <-> p2), so the surviving part
#      of E_2 is E_2^A = [E_2(p,p1,p2,p3,k) - E_2(k-p,p2,p1,p3,k)]/2.  Is Re E_2^A non-zero?  (=> k-even o G term in LS)
#  (4) our side: in the sum of Rows II.I and III.II the four-Wilson-line invariants Q_d (odd part) and Q (contact part)
#      cancel; checked with the tensors of app:fin_sym.
# usage: python3 -I ls_check.py
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colour import SUN, structures, POINTS, Ginv, Dinv, Finv

cr = lambda a, b: a[0] * b[1] - a[1] * b[0]
dt = lambda a, b: a @ b
n2 = lambda a: a @ a
epsk = lambda k: np.array([-k[1], k[0]])          # (eps k)^j = eps^{lj} k^l, eps^{12} = 1


def Gamma2(p, p1, p2, k):
    kp = k - p; skp = cr(k, p); a = abs(skp)
    A = (dt(p - p1, p1) / n2(p1) * dt(kp, p2) / (n2(kp) * n2(p2))
         + cr(p, p1) / n2(p1) * dt(kp - p2, p2) / (n2(kp) * n2(p2)) * skp / n2(p) * (1j * dt(k, kp) / a - 1))
    B = ((cr(p, p1) / n2(p1) * dt(k, p) / n2(p) + 0.5 * dt(p, p1) / n2(p1) * skp / n2(p)) * dt(kp, p2) / (n2(kp) * n2(p2))
         + 0.5 * dt(p - p1, p1) / n2(p1) * dt(kp - p2, p2) / (n2(p2) * n2(kp)) * skp / n2(p) * (1j * dt(p, p - k) / a + 1)
         + 0.5j * cr(p, p1) / n2(p1) * cr(kp, p2) / (n2(p2) * n2(kp)) * skp / n2(p) * (n2(k) + n2(p) - dt(p, k)) / a)
    return -A * k / n2(k) - B * epsk(k) / n2(k)


def Lv(p3, k): return p3 / n2(p3) + k / n2(k)        # L^j(p3, -k)


def E2_printed(p, p1, p2, p3, k):
    kp = k - p; skp = cr(k, p); a = abs(skp)
    t1 = -dt(p - p1, p1) / n2(p1) * dt(kp, p2) / (n2(kp) * n2(p2)) * dt(k + p3, p3) / (n2(k) * n2(p3))
    t2 = -cr(p, p1) / n2(p1) * dt(kp - p2, p2) / n2(p2) * skp / (n2(p) * n2(kp)) * (1j * dt(k, kp) / a - 1) * dt(k + p3, p3) / (n2(k) * n2(p3))
    t3 = -(dt(k, p) * cr(p, p1) / (n2(p) * n2(p1)) + 0.5 * skp * dt(p, p1) / (n2(p) * n2(p1))) * dt(kp, p2) / (n2(kp) * n2(p2)) * cr(k, p3) / (n2(k) * n2(p3))
    t4 = -0.5 * dt(p - p1, p1) / n2(p1) * dt(kp - p2, p2) / n2(p2) * skp / (n2(p) * n2(kp)) * (1j * dt(p, p - k) / a + 1) * cr(k, p3) / (n2(k) * n2(p3))
    t5 = 0.5 * cr(p, p1) / n2(p1) * cr(kp, p2) / n2(p2) * (dt(p, kp) - n2(k)) / (n2(p) * n2(kp)) * 1j * skp / a * cr(k, p3) / (n2(k) * n2(p3))
    return t1 + t2 + t3 + t4 + t5


def E1(p, p1, p2, k):
    return 0.5 * cr(p1, p) / (n2(p - p1) * n2(p1)) * (cr(p, p2) / (n2(p) * n2(p2)) + cr(p, k) / (n2(p) * n2(k)))


def main():
    rng = np.random.default_rng(5)
    R = lambda: rng.normal(size=2)
    # (1) transcription
    dev = 0; mx = 0
    for _ in range(2000):
        p, p1, p2, p3, k = R(), R(), R(), R(), R()
        a = E2_printed(p, p1, p2, p3, k); b = Gamma2(p, p1, p2, k) @ Lv(p3, k)
        dev = max(dev, abs(a - b) / max(abs(a), 1e-300)); mx = max(mx, abs(a))
    print('(1) LS E_2 printed vs Gamma_2.L: max rel. dev %.1e (2000 random points)' % dev)
    # (2) parity properties
    d1 = max(abs(E1(-p, -p1, -p2, -k) - E1(p, p1, p2, k)) for p, p1, p2, k in [(R(), R(), R(), R()) for _ in range(500)])
    print('(2) E_1 real; E_1(-all) - E_1(all): %.1e  => E_1 term + c.c. is even under k -> -k' % d1)
    c = [(R(), R(), R(), R(), R()) for _ in range(500)]
    dre = max(abs(E2_printed(-p, -p1, -p2, -p3, -k).real - E2_printed(p, p1, p2, p3, k).real) for p, p1, p2, p3, k in c)
    dim = max(abs(E2_printed(-p, -p1, -p2, -p3, -k).imag - E2_printed(p, p1, p2, p3, k).imag) for p, p1, p2, p3, k in c)
    print('    E_2 under reflection of all momenta: Re part changes by %.1e, Im part by %.1e (both reflection-even)' % (dre, dim))
    print('    => Re E_2 gives a k-even, Im E_2 (the sgn(k x p) terms) a k-odd contribution, as stated by LS.')
    # (3) d-projections
    for N in (3, 4):
        G = SUN(N)
        V = [G.adj(G.randV(rng)) for _ in range(3)]
        two = np.einsum('ade,dex,ac,xc->', G.f, G.d, V[0], V[1])
        three = np.einsum('cde,xyz,xd,ye,zc->', G.f, G.d, V[0], V[1], V[2])
        three_sw = np.einsum('cde,xyz,xd,ye,zc->', G.f, G.d, V[1], V[0], V[2])
        print('(3) SU(%d): 2-WL colour f^{ade} d^{dee2} (U U^T)^{a e2} = %.1e ;  3-WL colour = %+.5f, with (q1,p1)<->(q2,p2): %+.5f' % (
            N, abs(two), three, three_sw))
    reA = []; imA = []; ref = []
    for _ in range(4000):
        p, p1, p2, p3, k = R(), R(), R(), R(), R()
        e = E2_printed(p, p1, p2, p3, k); es = E2_printed(k - p, p2, p1, p3, k)
        reA.append(abs((e - es).real / 2)); imA.append(abs((e - es).imag / 2)); ref.append(abs(e))
    reA, imA, ref = map(np.array, (reA, imA, ref))
    print('    antisymmetrised E_2^A at 4000 random points: median |Re E_2^A|/|E_2| = %.3f, median |Im E_2^A|/|E_2| = %.3f;'
          ' fraction with |Re E_2^A| < 1e-12 |E_2|: %.3f' % (np.median(reA / ref), np.median(imA / ref), np.mean(reA < 1e-12 * ref)))
    # (4) Q_d and Q cancellation between Rows II.I and III.II
    for N in (3, 4):
        G = SUN(N); w4 = 0.0
        for _ in range(5):
            U = {p: G.adj(G.randV(rng)) for p in POINTS}
            S = structures(G, U)
            D1a = 2 * np.einsum('abc,abc->', S['II.I.N1a'][0], G.d); D32 = np.einsum('abc,abc->', S['III.II'][0], G.d)
            pred = Ginv(G, U, 'x', 'xp', 'z') + Ginv(G, U, 'x', 'z', 'wp')
            w4 = max(w4, abs(D1a / 2 + D32 - pred) / (abs(D1a) + abs(D32)))
            # contacts: y -> x and x' -> x
            for (a, b), predf in ((('y', 'x'), lambda U: N * Dinv(U, 'w', 'wp') - N * Dinv(U, 'xp', 'w') - Finv(G, U, 'x', 'z', 'wp') + Finv(G, U, 'x', 'xp', 'z')),
                                  (('xp', 'x'), lambda U: -N * Dinv(U, 'w', 'wp') + Finv(G, U, 'x', 'z', 'wp') - N * Dinv(U, 'x', 'z') + N * Dinv(U, 'x', 'w'))):
                Uc = dict(U); Uc[a] = U[b]; Sc = structures(G, Uc)
                P1a = np.einsum('abc,abc->', Sc['II.I.N1a'][0], G.f); P32 = np.einsum('abc,abc->', Sc['III.II'][0], G.f)
                # (A c i Phi): II.I N1a: A=-i, c=+1 (y=x) / -1 (x'=x);  III.II: A=-2i, c=-1/2 (y=x) / +1/2 (x'=x)
                comb = (P1a - P32) if a == 'y' else (P32 - P1a)
                w4 = max(w4, abs(comb - predf(Uc)) / (abs(P1a) + abs(P32)))
        print('(4) SU(%d): II.I + III.II, odd part: Delta_1a/2 + Delta_III.II = G(x,x\',z) + G(x,z,w\'); contacts y->x, x\'->x free of Q:'
              ' max rel. dev %.1e' % (N, w4))


if __name__ == '__main__':
    main()
