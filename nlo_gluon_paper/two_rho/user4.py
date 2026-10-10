"""The two-rho lists of 4rho.tex ('Two rho contribution', first part (x Nc) and second part dN(1), dN(2)).
Kernels: LO = KK(x'-w').KK(y-w); LOs = KK(y-w').KK(x'-w); K(a,b) = KK(a-z).KK(b-z).  Units of N."""
from core import term

NC = 3
LO = ('KK', "x'", "w'", 'y', 'w')
LOs = ('KK', 'y', "w'", "x'", 'w')
LOx = ('KK', "x'", "w'", "x'", 'w')     # KK(x'-w').KK(x'-w)


def K(a, b):
    return ('KK', a, 'z', b, 'z')


def U(i, j, p):
    return (i, j, p)


F1 = [  # first part, times Nc
    (2, [LO, K("x'", "x'")], [], [U('c', 'a', "x'"), U('c', 'b', 'w')], [('a', "x'"), ('b', 'y')]),
    (-2, [LO, K("x'", "x'")], [], [U('c', 'a', "x'"), U('c', 'b', 'y')], [('a', "x'"), ('b', 'y')]),
    (2, [LO, K('y', 'y')], [], [U('c', 'a', 'y'), U('c', 'b', "w'")], [('b', "x'"), ('a', 'y')]),
    (-2, [LO, K('y', 'y')], [], [U('c', 'a', 'y'), U('c', 'b', "x'")], [('b', "x'"), ('a', 'y')]),
    (1.5, [LO, K("x'", 'y')], [], [U('b', 'c', "x'"), U('a', 'c', 'y')], [('b', "x'"), ('a', 'y')]),
    (-8, [LOx, ('KK', 'y', 'z', "x'", 'z')], [], [U('b', 'c', "x'"), U('c', 'a', 'y')], [('b', "x'"), ('a', 'y')]),
]
F2 = [  # second part dN(1)
    (-1, [LO, K("x'", 'y')], [('a', 'b', 'c'), ('d', 'e', 'f')], [U('e', 'b', 'w'), U("b'", 'f', "x'"), U('d', 'a', 'y')], [("b'", "x'"), ('c', 'y')]),
    (3, [LO, K('y', 'y')], [('a', 'b', 'c'), ('d', 'e', "b'")], [U('a', "b'", 'y'), U('d', "b'", "w'"), U('e', 'b', 'z')], [("b'", "x'"), ('c', 'y')]),
    (0.5, [LO, K("x'", 'y')], [('a', 'b', 'c'), ('a', 'd', 'e')], [U('d', "b'", "w'"), U("b'", 'b', 'w')], [('e', "x'"), ('c', 'y')]),
    (-1, [LO, K("x'", 'y')], [('a', 'b', 'c'), ('a', 'd', 'e')], [U("b'", 'b', 'w'), U("b'", 'd', "x'")], [('e', "x'"), ('c', 'y')]),
    (2, [LO, K("x'", 'y')], [('a', 'b', 'c'), ('a', 'd', 'e')], [U("b'", 'b', 'y'), U("b'", 'd', "x'")], [('e', "x'"), ('c', 'y')]),
    (1, [LOs, K("x'", 'y')], [('a', 'b', 'c'), ('d', 'a', 'e')], [U("b'", 'b', "w'"), U("b'", 'd', "x'")], [('e', "x'"), ('c', 'y')]),
    (-3, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('a', 'f', 'z'), U('b', "d'", 'w'), U('d', 'f', 'y'), U("b'", "d'", "w'")], [("b'", "x'"), ('e', 'y')]),
    (2, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('a', 'f', 'z'), U('b', "d'", 'y'), U('d', 'f', 'y'), U("b'", "d'", "w'")], [("b'", "x'"), ('e', 'y')]),
    (-3, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('a', 'f', 'w'), U('b', "d'", 'y'), U('d', "d'", 'z'), U("b'", 'f', "w'")], [("b'", "x'"), ('e', 'y')]),
    (3, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('a', "d'", 'z'), U('d', 'f', 'w'), U('f', "b'", "w'"), U("d'", 'b', 'y')], [("b'", "x'"), ('e', 'y')]),
    (1, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'd', 'w'), U('f', "b'", "x'"), U("d'", 'a', 'z'), U("d'", 'b', 'y')], [("b'", "x'"), ('e', 'y')]),
    (0.5, [LO, K("x'", 'y')], [('a', 'b', 'c'), ("d'", 'f', 'l')], [U('a', "d'", 'z'), U('f', 'b', 'w'), U("b'", 'l', "x'")], [("b'", "x'"), ('c', 'y')]),
    (0.5, [LO, K("x'", 'y')], [('a', 'b', 'c'), ("d'", 'f', 'l')], [U('a', "d'", 'z'), U('b', 'f', "w'"), U("b'", 'l', 'y')], [('c', "x'"), ("b'", 'y')]),
    (-1, [LOs, K("x'", 'y')], [('a', 'b', 'c'), ("d'", 'f', 'l')], [U('f', 'b', 'y'), U("b'", 'l', "x'"), U("d'", 'a', 'z')], [("b'", "x'"), ('c', 'y')]),
    (-1, [LOs, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'b', 'y'), U('f', "b'", 'w'), U("d'", 'a', 'z'), U("d'", 'd', 'y')], [("b'", "x'"), ('e', 'y')]),
    (1, [LOs, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'b', 'y'), U('f', "b'", "x'"), U("d'", 'a', 'z'), U("d'", 'd', 'y')], [("b'", "x'"), ('e', 'y')]),
    (2, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('a', "d'", 'z'), U('b', 'f', 'y'), U('d', "d'", 'y'), U("b'", 'f', "x'")], [("b'", "x'"), ('e', 'y')]),
    (1, [LOs, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'd', 'y'), U('f', "b'", 'w'), U("d'", 'a', 'z'), U("d'", 'b', 'y')], [("b'", "x'"), ('e', 'y')]),
    (-2, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'd', 'y'), U('f', "b'", "w'"), U("d'", 'a', 'z'), U("d'", 'b', 'y')], [("b'", "x'"), ('e', 'y')]),
    (-2, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'd', 'y'), U('f', "b'", "x'"), U("d'", 'a', 'z'), U("d'", 'b', 'y')], [("b'", "x'"), ('e', 'y')]),
    (-1, [LOs, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U('f', 'd', 'y'), U('f', "b'", "x'"), U("d'", 'a', 'z'), U("d'", 'b', 'y')], [("b'", "x'"), ('e', 'y')]),
]
F3 = [  # second part dN(2)
    (1, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U("d'", 'a', 'w'), U("d'", "b'", "x'"), U('f', 'b', 'z'), U('f', 'd', 'y')], [("b'", "x'"), ('e', 'y')]),
    (1, [LO, K('y', 'y')], [('a', 'b', 'c'), ('c', 'd', 'e')], [U("d'", 'b', 'w'), U("d'", "b'", "x'"), U('f', 'a', 'y'), U('f', 'd', 'z')], [("b'", "x'"), ('e', 'y')]),
    (-1, [LOs, K("x'", 'y')], [('f', "d'", 'l'), ('d', 'e', "b'")], [U('d', 'f', 'z'), U('c', 'l', 'y'), U("d'", 'e', "x'")], [("b'", "x'"), ('c', 'y')]),
    (1, [LO, K('y', 'y')], [('a', 'b', 'c'), ('d', 'e', "b'")], [U('a', 'd', 'z'), U('b', "b'", 'y'), U('e', 'f', "x'")], [('f', "x'"), ('c', 'y')]),
]


# second part dN(1), term 2 as printed has b' four times; the index assignment that reproduces the exact
# result is  f^{abc} f^{deg} rho^{b'}(x') rho^c(y) U^{ag}(y) U^{db'}(w') U^{eb}(z)
F2_2_FIXED = (3, [LO, K('y', 'y')], [('a', 'b', 'c'), ('d', 'e', 'g')], [U('a', 'g', 'y'), U('d', "b'", "w'"), U('e', 'b', 'z')],
              [("b'", "x'"), ('c', 'y')])


def corrected_terms():
    """the lists of 4rho.tex with the two corrections: F1.6 (-8 Nc) dropped, F2.2 index fixed"""
    out, _ = terms()
    out = [t for t in out if t['tag'] != '4rho.tex F1.6']
    c, kern, fs, Us, rs = F2_2_FIXED
    out.append(term(c, Us, fs, rs, kern, '4rho.tex F2.2 (fixed)'))
    return out


def index_problems(fs, Us, rs):
    from collections import Counter
    cnt = Counter()
    for f in fs:
        cnt.update(f)
    for i, j, p in Us:
        cnt.update([i, j])
    for i, p in rs:
        cnt.update([i])
    return {k: v for k, v in cnt.items() if v != 2}


def terms(which=('F1', 'F2', 'F3')):
    out, bad = [], []
    for name, L, ncp in (('F1', F1, 1), ('F2', F2, 0), ('F3', F3, 0)):
        if name not in which:
            continue
        for n, (c, kern, fs, Us, rs) in enumerate(L):
            pb = index_problems(fs, Us, rs)
            tag = '4rho.tex %s.%d' % (name, n + 1)
            if pb:
                bad.append((tag, pb))
                continue
            out.append(term(c, Us, fs, rs, kern, tag, ncp=ncp))
    return out, bad


if __name__ == '__main__':
    ts, bad = terms()
    print(len(ts), 'terms;', len(bad), 'with index problems:')
    for b in bad:
        print('   ', b)
