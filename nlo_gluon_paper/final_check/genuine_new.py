"""The genuine two-rho terms of Sec. 1 of Evolution_terms.tex, encoded directly from the new text.
Units N = (1/(2pi)^3) g^4/(16 pi^5) (1/k+) log(vee/Lambda).
Dimensionally regulated brackets are written back as z-integrands (exact in d_perp = 2 - eps, MS-bar):
   2/eps + 2gamma + log((a-b)^2 mu^2/4)               = -(1/pi) int_z [K(a,b,z) - K(b,b,z)]
   -d/eps - d gamma - (d/2) log((a-b)^2 mu^2/4) - 1    =  (1/pi) int_z [K(a,b,z) - K(a,a,z)/2 - K(b,b,z)/2] + O(eps)
   and, in Row V 2b, the x-integral of the counterterm   2/eps + 2gamma + log((y-z)^2mu^2/4) = (1/pi) int_x [K(z,z,x) - K(y,z,x)].
Group I and Group III: + c.c. (conj_swap); Group II: no c.c."""
import numpy as np
import env  # noqa
from core import term, symmetrize, rename, conj_swap
import sets as S2

NCP = 1


def KK(a, b, c, d):
    return ('KK', a, b, c, d)


def with_cc(t):
    return symmetrize(t) + symmetrize(conj_swap(t))


def G1_II():
    """+ pi Nc [2/eps+...](x-w) KK(x'-w').KK(x-w) [U(x')rho(x') - U(w')rho(x')]^a [U(w)rho(x) - U(x)rho(x)]^a  + c.c."""
    out = []
    lo = KK("x'", "w'", 'x', 'w')
    zk = [(-1.0, KK('x', 'z', 'w', 'z')), (+1.0, KK('w', 'z', 'w', 'z'))]   # pi*(-1/pi)[K(x,w,z) - K(w,w,z)]
    for s1, p1 in ((1, "x'"), (-1, "w'")):
        for s2, p2 in ((1, 'w'), (-1, 'x')):
            for kc, k in zk:
                t = term(kc * s1 * s2, [('a', "b'", p1), ('a', 'b', p2)], [], [("b'", "x'"), ('b', 'x')], [lo, k], 'n:G1.II', ncp=NCP)
                out += with_cc(t)
    return out


def G1_V(include_uv_ct=True):
    out = []
    # 2a and 2b endpoint kernels: (coefficient, [kernel]) with the phase label
    k_z = [(1.0, [KK("x'", "w'", 'y', 'z'), KK('z', 'x', 'x', 'z')]), (0.5, [KK("x'", "w'", 'y', 'z'), KK('z', 'x', 'y', 'x')])]
    k_x = [(1.0, [KK("x'", "w'", 'y', 'x'), KK('z', 'x', 'x', 'z')]), (-0.5, [KK("x'", "w'", 'y', 'x'), KK('z', 'x', 'y', 'z')])]
    col2a = [(-1.0, [('d', 'b', 'x'), ('c', 'e', 'z')], [('a', 'd', 'c'), ('f', 'b', 'e')], 0), (1.0, [('a', 'f', 'x')], [], 1)]
    col2b = [(1.0, [('a', 'f', 'y')], [], 1), (-1.0, [('a', 'f', 'x')], [], 1)]
    canon = {'z': {'z': 'w', 'x': 'z'}, 'x': {'x': 'w'}}
    for cols, lab in ((col2a, '2a'), (col2b, '2b.fin')):
        for ph, ks in (('z', k_z), ('x', k_x)):
            for kc, kern in ks:
                for cc, Uc, fc, ncp in cols:
                    for s1, p1 in ((1, "x'"), (-1, "w'")):
                        t = term(-1.0 * kc * cc * s1, [('a', "b'", p1)] + Uc, fc, [("b'", "x'"), ('f', 'y')], kern,
                                 'n:G1.V.%s(%s)' % (lab, ph), ncp=ncp)
                        out += with_cc(rename(t, canon[ph]))
    if include_uv_ct:
        S = [(1.0, KK('z', 'x', 'z', 'x')), (-1.0, KK('y', 'x', 'z', 'x'))]   # K(z,z,x) - K(y,z,x)
        lo = KK("x'", "w'", 'y', 'z')
        # counterterm (finite part): -(g^4/8pi^5) -> -2 ;  colour Nc [U(y) - U(z)]
        # UV part: + g^4 Nc/(8 pi^4) E(y-z) -> +2 pi Nc (1/pi) int_x S = +2 Nc S
        for kc_tot, lab in ((-2.0, '2b.CT'), (+2.0, '2b.UV')):
            for kc, k in S:
                for s1, p1 in ((1, "x'"), (-1, "w'")):
                    for s2, p2 in ((1, 'y'), (-1, 'z')):
                        t = term(kc_tot * kc * s1 * s2, [('a', "b'", p1), ('a', 'f', p2)], [], [("b'", "x'"), ('f', 'y')], [lo, k],
                                 'n:G1.V.%s' % lab, ncp=1)
                        out += with_cc(rename(t, canon['z']))
    return out


def G2_IV():
    """(g^4/8pi^4) Nc [ -d/eps ... - 1 ](y'-y) KK(x-y).KK(x'-y') [U(y')-U(x')]^{e'c} [U(y)-U(x)]^{ec} rho^{e'}(x') rho^e(x); y->w, y'->w'"""
    out = []
    lo = KK('x', 'w', "x'", "w'")
    zk = [(1.0, KK('w', 'z', "w'", 'z')), (-0.5, KK('w', 'z', 'w', 'z')), (-0.5, KK("w'", 'z', "w'", 'z'))]
    for s1, p1 in ((1, "w'"), (-1, "x'")):
        for s2, p2 in ((1, 'w'), (-1, 'x')):
            for kc, k in zk:
                t = term(2.0 * kc * s1 * s2, [("e'", 'c', p1), ('e', 'c', p2)], [], [("e'", "x'"), ('e', 'x')], [lo, k], 'n:G2.IV', ncp=1)
                out += symmetrize(t)
    return out


def G3_III():
    """-(g^4/8pi^5) [KK(x-w).KK(x'-y') KK(z-w).KK(y'-z) + 1/2 KK(x-z).KK(y'-z) KK(x-w).KK(x'-y')] x colour + c.c.; y' -> w'"""
    out = []
    ks = [(1.0, [KK('x', 'w', "x'", "w'"), KK('z', 'w', "w'", 'z')]), (0.5, [KK('x', 'z', "w'", 'z'), KK('x', 'w', "x'", "w'")])]
    left = [(+1, [("e'", "c'", "w'"), ("d'", 'b', 'z')]), (-1, [('b', "d'", 'z'), ("c'", "e'", "x'")])]
    right = [(+1, [('a', 'c', 'w')], [('c', 'b', 'd')], ('d', 'x')), (-1, [('b', 'c', 'z'), ('d', 'e', 'x')], [('a', 'c', 'd')], ('e', 'x'))]
    for kc, kern in ks:
        for s1, Ul in left:
            for s2, Ur, fr, rr in right:
                t = term(-2.0 * kc * s1 * s2, Ul + Ur, [("c'", "d'", 'a')] + fr, [("e'", "x'"), rr], kern, 'n:G3.III2')
                out += with_cc(t)
    return out


def genuine_new():
    A = S2.A_sets()
    return {'g:I.II': G1_II(), 'g:I.V': G1_V(), 'g:II.I': A['g:II.I'], 'g:II.IV': G2_IV(), 'g:III.III': G3_III()}


if __name__ == '__main__':
    from core import Config, pointwise_total
    A = S2.A_sets()
    N = genuine_new()
    for k in N:
        r = []
        for seed in (3, 7, 11):
            cfg = Config(seed, True)
            a = pointwise_total(N[k], cfg)[0].real
            b = pointwise_total(A[k], cfg)[0].real
            r.append(abs(a - b) / (1 + abs(b)))
        print('%-10s new text vs earlier corrected encoding: %.2e' % (k, max(r)))
    uvct = [t for t in G1_V() if t['tag'].startswith(('n:G1.V.2b.UV', 'n:G1.V.2b.CT'))]
    cfg = Config(3, True)
    print('Row V 2b: UV + counterterm (pointwise in x): %.2e (scale %.2e)' % pointwise_total(uvct, cfg))
