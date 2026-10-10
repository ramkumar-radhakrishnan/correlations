"""Region-P limits (coefficient of log(vee/k+) + its T-companion) of the rows of the notes.

For every row,  d3N_row = c_T * dY_T + c_P * dY_P + finite,  dY_T = log(k+/Lambda), dY_P = log(vee/k+).
P_row = lim_{p+ -> infinity} p+ x (integrand)  (units N_P = N with log(vee/Lambda) -> log(vee/k+)).
With the user's rule log(Lambda/k+) = -log(vee/Lambda) + log(vee/k+), the coefficient of log(vee/k+)
of a row is  L_row = P_row - T_row.

Three-rho rows with P_row != T_row (Sources, same colour/word as the T sources of sources.py):
  GI.III.P   i [ -K_b - 2 K_h ]                   (Group I Row III, from the p+ -> vee ends)
  GI.IV.P    2i K_PB                                (Group I Row IV)
  GII.I1.P   -i K_PB  ({BC}A)    GII.I2.P  +i Kbar_PB  (A{BC})   (Group II Row I part 1)
  GIII.I.P, GIII.II.P   -2i K_PB                    (Group III Rows I, II)
  GIII.III.P -i K_PC (A{BC}),  GIII.V.P, GIII.VI.P  -2i K_PC    (Group III Rows III(1), V, VI),
             colour with U(y') -> U(z); the overall sign is that of Evolution_three_rho.tex
Two-rho rows: G2_I_P (K_C/(b b')), G2_IV_P, G3_III_P (sign as corrected, fix F5).
Rows with P_row = 0: Group I Row II (both parts), Group I Row V (all parts).
"""
import env  # noqa
import sympy as sp
from core import term, symmetrize, conj_swap
from sources import SOURCES, Source, KK
import users as US

I = sp.I
R = sp.Rational
SRC = {s.key: s for s in SOURCES}

K_PB = [(1, [KK("x'", "w'", 'z', 'w'), KK('y', 'z', 'x', 'z')]),
        (-R(1, 2), [KK("x'", "w'", 'x', 'w'), KK('x', 'z', 'y', 'z')])]
K_PBbar = [(1, [KK("x'", 'w', 'z', "w'"), KK('y', 'z', 'x', 'z')]),
           (-R(1, 2), [KK("x'", 'w', 'x', "w'"), KK('x', 'z', 'y', 'z')])]
K_III = [(-1, [KK("x'", "w'", 'x', 'w'), KK('x', 'z', 'y', 'z')]),
         (-2, [KK("x'", "w'", 'z', 'w'), KK('y', 'z', 'x', 'z')])]
K_PC = [(1, [KK("w'", 'z', 'y', 'w'), KK("x'", 'z', 'x', 'z')])]


def yz(sub):
    """U(y') -> U(z) in the colour of a sub-term"""
    out = []
    for s, col, word in sub:
        out.append((s, [('U', u[1], u[2], 'z') if (u[0] == 'U' and u[3] == "y'") else u for u in col], word))
    return out


def mk(key, base, pref, kname, kern, subs=None, form=None, phase=("w'", 'w')):
    b = SRC[base]
    return Source(key, b.title + r' (region P)', 'P', pref, phase, kname, b.kargs, kern,
                  form or b.form, b.cc, subs if subs is not None else b.subterms, b.rho_labels, ["w'", 'w', 'z'])


PSOURCES = [
    mk('GI.III.P', 'GI.III', I, r'K^{P}_{\rm I.III}', K_III),
    mk('GI.IV.P', 'GI.IV', 2 * I, r'K^{P}_{\mathfrak B}', K_PB),
    mk('GII.I1.P', 'GII.I1', -I, r'K^{P}_{\mathfrak B}', K_PB),
    mk('GII.I2.P', 'GII.I2', I, r'\bar K^{P}_{\mathfrak B}', K_PBbar),
    mk('GIII.I.P', 'GIII.I', -2 * I, r'K^{P}_{\mathfrak B}', K_PB),
    mk('GIII.II.P', 'GIII.II', -2 * I, r'K^{P}_{\mathfrak B}', K_PB),
    mk('GIII.III.P', 'GIII.III', -I, r'K^{P}_{\mathfrak C}', K_PC, yz(SRC['GIII.III'].subterms)),
    mk('GIII.V.P', 'GIII.V', -2 * I, r'K^{P}_{\mathfrak C}', K_PC, yz(SRC['GIII.V'].subterms)),
    mk('GIII.VI.P', 'GIII.VI', -2 * I, r'K^{P}_{\mathfrak C}', K_PC, yz(SRC['GIII.VI'].subterms)),
]
PSRC = {s.key: s for s in PSOURCES}


def p_terms(key, keep=(2,)):
    return US.source_terms(PSRC[key], keep)


# ---------------------------------------------------------------- genuine two-rho terms, region P
U = US.U


def G2_I_P():
    """2 N_P [K_C/(b b')] x colour (no c.c.)"""
    out = []
    kerns = [(1.0, [('KK', "x'", 'z', 'x', 'z'), ('KK', 'z', 'w', 'z', "w'")]),
             (-0.5, [('KK', "x'", 'z', 'x', 'z'), ('KK', "x'", "w'", 'z', 'w')]),
             (-0.5, [('KK', 'x', 'w', 'z', "w'"), ('KK', 'x', 'z', "x'", 'z')]),
             (0.25, [('KK', 'x', 'z', "x'", 'z'), ('KK', "x'", "w'", 'x', 'w')])]
    cols = [(+1, [U('a', "c'", "w'"), U('a', 'c', 'w')], [("c'", 'b', "d'"), ('c', 'b', 'd')], [("d'", "x'"), ('d', 'x')]),
            (-1, [U('a', "c'", "w'"), U('b', 'c', 'z'), U('d', 'e', 'x')], [("c'", 'b', "d'"), ('a', 'c', 'd')], [("d'", "x'"), ('e', 'x')]),
            (-1, [U('a', 'c', 'w'), U('b', "c'", 'z'), U("e'", "d'", "x'")], [('a', "c'", "d'"), ('c', 'b', 'd')], [("e'", "x'"), ('d', 'x')]),
            ('Nc', [U("e'", 'd', "x'"), U('d', 'e', 'x')], [], [("e'", "x'"), ('e', 'x')])]
    for kc, kern in kerns:
        for cc, Us, fs, rs in cols:
            if cc == 'Nc':
                t = term(2.0 * kc, Us, fs, rs, kern, 'P.II.I', ncp=1)
            else:
                t = term(2.0 * kc * cc, Us, fs, rs, kern, 'P.II.I')
            out += symmetrize(t)
    return out


def G2_IV_P():
    """2 Nc N_P K(w,w',z) K(x,x',z) [U(z)-U(x')]^{e'c} [U(z)-U(x)]^{ec} rho^{e'}(x') rho^e(x)  (no c.c.)"""
    out = []
    kern = [('KK', 'w', 'z', "w'", 'z'), ('KK', 'x', 'z', "x'", 'z')]
    for s1, p1 in ((1, 'z'), (-1, "x'")):
        for s2, p2 in ((1, 'z'), (-1, 'x')):
            t = term(2.0 * s1 * s2, [U("e'", 'c', p1), U('e', 'c', p2)], [], [("e'", "x'"), ('e', 'x')], kern, 'P.II.IV', ncp=1)
            out += symmetrize(t)
    return out


def G3_III_P(sign=-1.0):
    """sign x 2 N_P K(x',x,z) KK(w'-z).[KK(z-w) - KK(x-w)/2] x colour(y'->z) + c.c.   (sign -1: corrected, fix F5)"""
    out = []
    kerns = [(1.0, [('KK', "x'", 'z', 'x', 'z'), ('KK', "w'", 'z', 'z', 'w')]),
             (-0.5, [('KK', "x'", 'z', 'x', 'z'), ('KK', "w'", 'z', 'x', 'w')])]
    left = [(+1, [U("e'", "c'", 'z'), U("d'", 'b', 'z')]), (-1, [U('b', "d'", 'z'), U("c'", "e'", "x'")])]
    right = [(+1, [U('a', 'c', 'w')], [('c', 'b', 'd')], ('d', 'x')),
             (-1, [U('b', 'c', 'z'), U('d', 'e', 'x')], [('a', 'c', 'd')], ('e', 'x'))]
    for kc, kern in kerns:
        for s1, Ul in left:
            for s2, Ur, fr, rr in right:
                t = term(sign * 2.0 * kc * s1 * s2, Ul + Ur, [("c'", "d'", 'a')] + fr, [("e'", "x'"), rr], kern, 'P.III.III2')
                out += symmetrize(t) + symmetrize(conj_swap(t))
    return out


def genuineP():
    return dict(G2_I=G2_I_P(), G2_IV=G2_IV_P(), G3_III=G3_III_P())
