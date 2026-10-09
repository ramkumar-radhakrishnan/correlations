"""The four-rho terms of the notes (Eq. (4rho) of 4rho.tex), transcribed line by line.

Common factor N * KK(x'-w').KK(u-w) KK(x-z).KK(y-z); charges at x', u, x, y.
Each line: coefficient * [conjugate bracket] x [amplitude bracket].
A bracket is a list of (sign, U-factors, list of rho words); a list of several
words means their sum (anticommutators and (AB+BA) are written out)."""

def U(i, j, p):
    return ('U', i, j, p)


CONJ_LO = [(+1, [U('a', "b'", "x'")], [[("b'", "x'")]]),
           (-1, [U('a', "b'", "w'")], [[("b'", "x'")]])]
CONJ_X1 = [(+1, [U('a', "b'", "x'"), U("c'", "a'", 'x'), U("c'", 'b', 'z')], [[("b'", "x'"), ("a'", 'x')]]),
           (-1, [U('b', "d'", 'z'), U("d'", "a'", 'x'), U("c'", 'a', "w'")], [[("c'", "x'"), ("a'", 'x')]])]
CONJ_X2 = [(+1, [U("c'", "d'", 'x'), U('a', "e'", "x'"), U('b', "c'", 'z')], [[("d'", 'x'), ("e'", "x'")]]),
           (-1, [U('a', "b'", "x'")], [[('b', 'x'), ("b'", "x'")]])]
CONJ_B = [(+1, [U("c'", 'a', "w'")], [[('b', 'x'), ("c'", "x'")], [("c'", "x'"), ('b', 'x')]]),
          (-1, [U("c'", "d'", 'x'), U('a', 'd', "x'"), U('b', "c'", 'z')], [[("d'", 'x'), ('d', "x'")], [('d', "x'"), ("d'", 'x')]])]

AMP_L1 = [(+1, [U('a', 'b', 'u')], [[('b', 'u'), ('c', 'x'), ('c', 'y')]]),
          (-1, [U('a', 'b', 'w'), U('c', 'd', 'x'), U('c', 'e', 'y')], [[('d', 'x'), ('e', 'y'), ('b', 'u')]])]
AMP_L2 = [(+1, [U('a', 'b', 'w')], [[('c', 'x'), ('c', 'y'), ('b', 'u')]]),
          (-1, [U('c', 'd', 'x'), U('c', 'e', 'y'), U('a', 'b', 'u')], [[('d', 'x'), ('e', 'y'), ('b', 'u')]])]
AMP_L3 = [(+1, [U('a', 'c', 'w'), U('b', 'd', 'z'), U('b', 'e', 'y')], [[('e', 'y'), ('d', 'x'), ('c', 'u')], [('e', 'y'), ('c', 'u'), ('d', 'x')]]),
          (-1, [U('b', 'c', 'y'), U('b', 'd', 'x'), U('a', 'e', 'u')], [[('c', 'y'), ('d', 'x'), ('e', 'u')], [('c', 'y'), ('e', 'u'), ('d', 'x')]])]
AMP_L4 = [(+1, [U('b', 'c', 'x'), U('a', 'd', 'u'), U('b', 'e', 'y')], [[('c', 'x'), ('d', 'u'), ('e', 'y')], [('d', 'u'), ('c', 'x'), ('e', 'y')]]),
          (-1, [U('b', 'c', 'x'), U('a', 'd', 'u'), U('b', 'e', 'z')], [[('c', 'x'), ('d', 'u'), ('e', 'y')], [('d', 'u'), ('c', 'x'), ('e', 'y')]])]
AMP_B = [(+1, [U('a', 'c', 'w')], [[('c', 'u'), ('b', 'y')], [('b', 'y'), ('c', 'u')]]),
         (-1, [U('b', 'c', 'z'), U('c', 'e', 'y'), U('a', "b'", 'u')], [[('e', 'y'), ("b'", 'u')], [("b'", 'u'), ('e', 'y')]])]
AMP_B_e = [(+1, [U('a', 'c', 'w')], [[('c', 'u'), ('b', 'y')], [('b', 'y'), ('c', 'u')]]),
           (-1, [U('b', 'c', 'z'), U('c', 'e', 'y'), U('a', "e'", 'u')], [[('e', 'y'), ("e'", 'u')], [("e'", 'u'), ('e', 'y')]])]
AMP_B_a = [(+1, [U('a', 'c', 'w')], [[('c', 'u'), ('b', 'y')], [('b', 'y'), ('c', 'u')]]),
           (-1, [U('b', 'c', 'z'), U('c', 'e', 'y'), U('a', "a'", 'u')], [[('e', 'y'), ("a'", 'u')], [("a'", 'u'), ('e', 'y')]])]
AMP_X1 = [(+1, [U('b', 'c', 'z'), U('c', 'd', 'y'), U('a', 'e', 'u')], [[('d', 'y'), ('e', 'u')]]),
          (-1, [U('a', 'c', 'w'), U('d', 'e', 'y'), U('d', 'b', 'z')], [[('e', 'y'), ('c', 'u')]])]
AMP_X2 = [(+1, [U('b', 'c', 'z'), U('a', 'd', 'u'), U('c', 'e', 'y')], [[('d', 'u'), ('e', 'y')]]),
          (-1, [U('a', 'c', 'u')], [[('c', 'u'), ('b', 'y')]])]

# (label, row, coefficient, conj bracket, amp bracket, add c.c.)
LINES = [
    ('L1', 'G1', -1, CONJ_LO, AMP_L1, True),
    ('L2', 'G1', +1, CONJ_LO, AMP_L2, True),
    ('L3', 'G1', -1, CONJ_LO, AMP_L3, True),
    ('L4', 'G1', -1, CONJ_LO, AMP_L4, True),
    ('L5', 'G2-I', 0.5, CONJ_B, AMP_B, False),
    ('L6', 'G2-II', 2, CONJ_X1, AMP_X1, False),
    ('L7', 'G2-III', 2, CONJ_X2, AMP_X2, False),
    ('L8', 'G3-I', 1, CONJ_X1, AMP_B_e, True),
    ('L9', 'G3-II', 1, CONJ_X2, AMP_B_a, True),
    ('L10', 'G3-IV', 2, CONJ_X2, AMP_X1, True),
]
KERNEL = [(("x'", "w'"), ('u', 'w')), (('x', 'z'), ('y', 'z'))]


def words():
    """all four-rho words: (coef, row, line, colour U's, rho word, kernel, cc)."""
    from engine import check_contracted
    out = []
    for lab, row, c, conj, amp, cc in LINES:
        for k1, (s1, U1, W1) in enumerate(conj):
            for k2, (s2, U2, W2) in enumerate(amp):
                for w1 in W1:
                    for w2 in W2:
                        facs = U1 + U2
                        word = w1 + w2
                        check_contracted(facs + [('rho', i, p) for i, p in word], lab)
                        out.append((c * s1 * s2, row, '%s.%d.%d' % (lab, k1 + 1, k2 + 1), facs, word, KERNEL, False))
                        if cc:
                            sw = {'w': "w'", "w'": 'w'}
                            g = lambda p: sw.get(p, p)
                            facs_c = [('U', u[1], u[2], g(u[3])) for u in facs]
                            word_c = [(i, p) for i, p in word[::-1]]
                            kern_c = [((g(a), g(b)), (g(cc_), g(d))) for (a, b), (cc_, d) in KERNEL]
                            out.append((c * s1 * s2, row, '%s.%d.%d*' % (lab, k1 + 1, k2 + 1), facs_c, word_c, kern_c, True))
    return out
