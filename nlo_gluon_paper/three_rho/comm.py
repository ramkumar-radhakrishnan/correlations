"""The commutator -[Nbar_2, Abar^(1)] of Row II times the conjugate LO amplitude."""
import numpy as np
import compare_ref as CR


def comm_terms():
    """-[Nbar_2, Abar^(1)] in Row II times the conjugate LO amplitude, three-rho (symmetric) part:
       -i N int KK(x'-w').KK(s-w) K(u,s,z) f^{eag}[U^{ab'}(x')-U^{ab'}(w')] rho^{b'}(x') {Ubar rho(u)^e, Ubar rho(s)^g}
       -> symmetric part: 2 x (-i N) with commuting charges.  Plus c.c."""
    out = []
    kern = [(("x'", "w'"), ('s', 'w')), (('u', 'z'), ('s', 'z'))]
    for sgn, P in ((+1, "x'"), (-1, "w'")):
        col = [('f', 'e', 'a', 'g'), ('U', 'a', "b'", P), ('U', 'e', 'h', 'u'), ('U', 'g', 'l', 's')]
        rh = [("b'", "x'"), ('h', 'u'), ('l', 's')]
        c = -2j * sgn
        out.append(CR.C3(c, col, rh, kern, 'COMM'))
        sw = {'w': "w'", "w'": 'w'}
        g = lambda p: sw.get(p, p)
        out.append(CR.C3(np.conj(c), [('U', x[1], x[2], g(x[3])) if x[0] == 'U' else x for x in col], rh,
                         [((g(a), g(b)), (g(cc), g(d))) for (a, b), (cc, d) in kern], 'COMM*'))
    return out

