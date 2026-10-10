"""All two-rho terms proportional to log(vee/Lambda) of the notes.

R3 : two-rho parts of the three-rho terms (four-rho remainders A1..B2, Groups I-III), by exact
     Weyl symmetrization of the ordered products (the content of ThreeRho_to_TwoRho.pdf);
R4 : two-rho parts of the four-rho terms of Eq. (4rho) (incl. c.c.);
G  : the genuine two-rho terms of Evolution_two_rho.tex (Group I Rows II, V; Group II Rows I, IV;
     Group III Row III part 2);
corrections of the three-rho note (Row II alpha part, Row III part 2 -> D).
All coefficients in units of N."""
import sympy as sp
import numpy as np
from core import term, symmetrize, rename, conj_swap
from sources import SOURCES, Source, KK
from cancel3 import CANON
import classes3 as C3
import fourrho

I = sp.I


def cnum(x):
    return complex(sp.N(x))


def native_kernel(src):
    out = []
    for kc, fs in C3.monomials(src):
        kern = []
        for f in fs:
            if f[0] == 'KK':
                kern.append(f)
            elif f[0] == 'raw' and 'Phi' in f[1]:
                kern.append(('PHI', f[2][0], 'z', 'w'))
            else:
                raise ValueError(f)
        out.append((kc, kern))
    return out


def ordered_words(form, word):
    if form == 'plain':
        return [list(word)]
    if form == 'A{BC}':
        A, B, C = word
        return [[A, B, C], [A, C, B]]
    if form == '{BC}A':
        B, C, A = word
        return [[B, C, A], [C, B, A]]
    raise ValueError(form)


def source_terms(src, keep=(2,), monos=None, subs=None, cc_list=None):
    """ordered terms of one source (each kernel monomial separately), symmetrized."""
    out = []
    base = dict(CANON.get(src.key, {}))
    kern_list = native_kernel(src)
    for n, (sign, colour, word) in enumerate(src.subterms):
        if subs is not None and n + 1 not in subs:
            continue
        Us = [(u[1], u[2], u[3]) for u in colour if u[0] == 'U']
        fs = [tuple(u[1:]) for u in colour if u[0] == 'f']
        for m, (kc, kern) in enumerate(kern_list):
            if monos is not None and m not in monos:
                continue
            c = cnum(src.pref * sign * kc)
            tag = '%s.%d%s' % (src.key, n + 1, 'abcd'[m] if len(kern_list) > 1 else '')
            for w in ordered_words(src.form, word):
                t = term(c, Us, fs, w, kern, tag)
                t = rename(t, base)
                variants = [t]
                if src.cc:
                    variants.append(conj_swap(t))
                for v in variants:
                    if cc_list is not None and (v['tag'].endswith('*') not in cc_list):
                        continue
                    out += symmetrize(v, keep)
    return out


def R3(keep=(2,)):
    out = []
    for src in SOURCES:
        out += source_terms(src, keep)
    return out


def R4(keep=(2,)):
    out = []
    for c, row, lab, Us, word, kern, cc in fourrho.words():
        kk = [('KK', a, b, c_, d) for (a, b), (c_, d) in kern]
        t = term(c, [(u[1], u[2], u[3]) for u in Us], [], word, kk, '4rho:%s:%s' % (row, lab))
        out += symmetrize(t, keep)
    return out


# ---------------------------------------------------------------- genuine two-rho terms
def U(i, j, p):
    return (i, j, p)


def LOcol(xp, wp, x, w, s1=1, s2=1):
    """[U(xp)-U(wp)]^{a b'} rho^{b'}(xp) [U(x)-U(w)]^{a b} rho^b(x) as a list of (sign, Us)"""
    return [(+1, [U('a', "b'", xp), U('a', 'b', x)]), (-1, [U('a', "b'", wp), U('a', 'b', x)]),
            (-1, [U('a', "b'", xp), U('a', 'b', w)]), (+1, [U('a', "b'", wp), U('a', 'b', w)])]


def G1_II():
    """-pi Nc KK(x'-w').KK(x-w) BG1(x,w) [U(x')-U(w')]rho(x') [U(w)-U(x)]rho(x) + c.c."""
    out = []
    kern = [('KK', "x'", "w'", 'x', 'w'), ('BG1', 'x', 'w')]
    for s, Us in LOcol("x'", "w'", 'x', 'w'):
        part = 'dir' if Us[1][2] == 'w' else 'bar'
        t = term(-np.pi * (-s), Us, [], [("b'", "x'"), ('b', 'x')], kern, 'G.I.II.' + part, ncp=1)
        out += symmetrize(t) + symmetrize(conj_swap(t))
    return out


NC = 3


def G1_V():
    """Row V, parts 2a + log part of 2b:
    -N [e^{-ik(w'-z)} K_low + e^{-ik(w'-x)} K_up] [U(x')-U(w')]^{ab'} rho^{b'}(x') C^{af} rho^f(y) + c.c.
    C^{af} = -f^{adc} f^{fbe} U^{db}(x) U^{ce}(z) + Nc U^{af}(x)  (2a)  +  Nc U^{af}(y) - Nc U^{af}(x)  (2b)"""
    out = []
    K_low = [(1.0, [('KK', "x'", "w'", 'y', 'z'), ('KK', 'z', 'x', 'x', 'z')]),
             (0.5, [('KK', "x'", "w'", 'y', 'z'), ('KK', 'z', 'x', 'y', 'x')])]
    K_up = [(1.0, [('KK', "x'", "w'", 'y', 'x'), ('KK', 'z', 'x', 'x', 'z')]),
            (-0.5, [('KK', "x'", "w'", 'y', 'x'), ('KK', 'z', 'x', 'y', 'z')])]
    cols = [('2a.dir', -1, [U('d', 'b', 'x'), U('c', 'e', 'z')], [('a', 'd', 'c'), ('f', 'b', 'e')], 0),
            ('2a.x', 1, [U('a', 'f', 'x')], [], 1),
            ('2b.bar', 1, [U('a', 'f', 'y')], [], 1),
            ('2b.x', -1, [U('a', 'f', 'x')], [], 1)]
    for phase, Kset, canon in (('z', K_low, {'z': 'w', 'x': 'z'}), ('x', K_up, {'x': 'w'})):
        for kc, kern in Kset:
            for part, cc, Uc, fc, ncp in cols:
                for s1, p1 in ((1, "x'"), (-1, "w'")):
                    t = term(-1.0 * kc * cc * s1, [U('a', "b'", p1)] + Uc, fc, [("b'", "x'"), ('f', 'y')], kern,
                             'G.I.V%s(%s)' % (part, phase), ncp=ncp)
                    t = rename(t, canon)
                    out += symmetrize(t) + symmetrize(conj_swap(t))
    return out


def G2_I():
    """2N KK(x'-w').KK(x-w)[K(w,w',z) - K(x',w,z)/2 - K(x,w',z)/2 + K(x,x',z)/4] x colour (no c.c.)"""
    out = []
    lo = ('KK', "x'", "w'", 'x', 'w')
    kerns = [(1.0, ('KK', 'w', 'z', "w'", 'z')), (-0.5, ('KK', "x'", 'z', 'w', 'z')),
             (-0.5, ('KK', 'x', 'z', "w'", 'z')), (0.25, ('KK', 'x', 'z', "x'", 'z'))]
    cols = [(+1, [U('a', "c'", "w'"), U('a', 'c', 'w')], [("c'", 'b', "d'"), ('c', 'b', 'd')], [("d'", "x'"), ('d', 'x')]),
            (-1, [U('a', "c'", "w'"), U('b', 'c', 'z'), U('d', 'e', 'x')], [("c'", 'b', "d'"), ('a', 'c', 'd')], [("d'", "x'"), ('e', 'x')]),
            (-1, [U('a', 'c', 'w'), U('b', "c'", 'z'), U("e'", "d'", "x'")], [('a', "c'", "d'"), ('c', 'b', 'd')], [("e'", "x'"), ('d', 'x')]),
            ('Nc', [U("e'", 'd', "x'"), U('d', 'e', 'x')], [], [("e'", "x'"), ('e', 'x')])]
    for kc, kf in kerns:
        for cc, Us, fs, rs in cols:
            if cc == 'Nc':
                t = term(2.0 * kc, Us, fs, rs, [lo, kf], 'G.II.I', ncp=1)
            else:
                t = term(2.0 * kc * cc, Us, fs, rs, [lo, kf], 'G.II.I')
            out += symmetrize(t)
    return out


def G2_IV():
    """2 Nc N ELL(w,w') KK(x-w).KK(x'-w') [U(w')-U(x')]^{e'c} [U(w)-U(x)]^{ec} rho^{e'}(x') rho^e(x)  (no c.c.)"""
    out = []
    kern = [('KK', 'x', 'w', "x'", "w'"), ('ELL', 'w', "w'")]
    for s1, p1 in ((1, "w'"), (-1, "x'")):
        for s2, p2 in ((1, 'w'), (-1, 'x')):
            t = term(2.0 * s1 * s2, [U("e'", 'c', p1), U('e', 'c', p2)], [], [("e'", "x'"), ('e', 'x')], kern, 'G.II.IV', ncp=1)
            out += symmetrize(t)
    return out


def G3_III():
    """2N KK(x-w).KK(x'-w')[-K(w,w',z) + K(x,w',z)/2] x colour + c.c.  (y' -> w')"""
    out = []
    lo = ('KK', 'x', 'w', "x'", "w'")
    kerns = [(1.0, ('KK', 'z', 'w', "w'", 'z')), (0.5, ('KK', 'x', 'z', "w'", 'z'))]
    left = [(+1, [U("e'", "c'", "w'"), U("d'", 'b', 'z')]), (-1, [U('b', "d'", 'z'), U("c'", "e'", "x'")])]
    right = [(+1, [U('a', 'c', 'w')], [('c', 'b', 'd')], ('d', 'x')),
             (-1, [U('b', 'c', 'z'), U('d', 'e', 'x')], [('a', 'c', 'd')], ('e', 'x'))]
    for kc, kf in kerns:
        for s1, Ul in left:
            for s2, Ur, fr, rr in right:
                t = term(2.0 * kc * s1 * s2, Ul + Ur, [("c'", "d'", 'a')] + fr, [("e'", "x'"), rr], [lo, kf], 'G.III.III2')
                out += symmetrize(t) + symmetrize(conj_swap(t))
    return out


def genuine():
    return dict(G1_II=G1_II(), G1_V=G1_V(), G2_I=G2_I(), G2_IV=G2_IV(), G3_III=G3_III())


# ---------------------------------------------------------------- corrections of the three-rho note
def D_terms(keep=(2,)):
    import cancel_data as CD
    return source_terms(CD.D_SOURCE, keep)


def RowII_alpha_U(keep=(2,)):
    """the alpha part (monomials b, d) of Row II in the U(w)A3 sub-terms 1, 3 (to be removed)"""
    src = [s for s in SOURCES if s.key == 'GI.II'][0]
    return source_terms(src, keep, monos=(1, 3), subs=(1, 3))


def RowIII_part2(keep=(2,)):
    src = [s for s in SOURCES if s.key == 'GI.III'][0]
    return source_terms(src, keep, monos=(1,))


def neg(ts):
    return [dict(t, c=-t['c']) for t in ts]


def RowII_alpha_all(keep=(2,)):
    src = [s for s in SOURCES if s.key == 'GI.II'][0]
    return source_terms(src, keep, monos=(1, 3))


def commutator_term(keep=(2,)):
    """-i N KK(x'-w').KK(s-w) K(u,s,z) f^{eag} [U^{ab'}(x')-U^{ab'}(w')] rho^{b'}(x') {U^{eh}(u) rho^h(u), U^{gl}(s) rho^l(s)} + c.c."""
    out = []
    kern = [('KK', "x'", "w'", 's', 'w'), ('KK', 'u', 'z', 's', 'z')]
    for sgn, P in ((+1, "x'"), (-1, "w'")):
        Us = [('a', "b'", P), ('e', 'h', 'u'), ('g', 'l', 's')]
        for w in ([("b'", "x'"), ('h', 'u'), ('l', 's')], [("b'", "x'"), ('l', 's'), ('h', 'u')]):
            t = term(-1j * sgn, Us, [('e', 'a', 'g')], w, kern, 'COMM')
            out += symmetrize(t, keep) + symmetrize(conj_swap(t), keep)
    return out


def D_terms_order(order, keep=(2,)):
    """D with the charges ordered as given (permutation of the Row III order x', x, y)."""
    import cancel_data as CD
    src = CD.D_SOURCE
    out = []
    for n, (sign, colour, word) in enumerate(src.subterms):
        Us = [(u[1], u[2], u[3]) for u in colour if u[0] == 'U']
        fs = [tuple(u[1:]) for u in colour if u[0] == 'f']
        kc, kern = native_kernel(src)[0]
        w = [word[k] for k in order]
        t = term(cnum(src.pref * sign * kc), Us, fs, w, kern, 'D.%d' % (n + 1))
        out += symmetrize(t, keep) + symmetrize(conj_swap(t), keep)
    return out
