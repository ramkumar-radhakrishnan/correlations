"""All three-rho terms proportional to log(vee/Lambda), transcribed from
  (i)  the 'Three rho contributions' of the four-rho notes (two and three Wilson lines),
  (ii) Evolution_three_rho.tex (genuine three-rho terms, row by row).
Each source is written as
    prefactor * N * int e^{-ik(P'-P)}  K(rho positions)  sum_n sign_n  C_n  (rho word)
with N = 1/(2pi)^3 * g^4/(16 pi^5) * 1/k^+ * log(vee/Lambda).
"""
import sympy as sp
from engine import I

R = sp.Rational


class Source:
    def __init__(self, key, title, origin, pref, phase, kname, kargs, kexpl, form, cc, subterms,
                 rho_labels, fixed_labels, remark=''):
        self.key = key
        self.title = title
        self.origin = origin
        self.pref = pref            # sympy number, in units of N
        self.phase = phase          # (P', P): e^{-ik(P'-P)}
        self.kname = kname          # LaTeX name of kernel, e.g. 'K_1'
        self.kargs = kargs          # list of rho-position labels the kernel depends on
        self.kexpl = kexpl          # explicit kernel: list of (coeff, [factors])
        self.form = form            # 'plain', 'A{BC}', '{BC}A'
        self.cc = cc                # True: '+ c.c.' (Groups I and III)
        self.subterms = subterms    # list of (sign, [colour factors w/o rho], [(idx,pos) word])
        self.rho_labels = rho_labels
        self.fixed_labels = fixed_labels
        self.remark = remark


def KK(p, q, r, s):
    return ('KK', p, q, r, s)


def raw(fmt, *args):
    return ('raw', fmt, list(args))


SOURCES = []

# ---------------------------------------------------------------- four-rho remainders
K1 = [(1, [KK("x'", "w'", 'y', 'w'), KK('x', 'z', 'y', 'z')])]
K2 = [(1, [KK("x'", "w'", 'y', 'w'), KK('x', 'z', "x'", 'z')])]
K3 = [(1, [KK("x'", "w'", 'x', 'w'), KK('y', 'z', "x'", 'z')])]
RL = ["x'", 'x', 'y']

SOURCES.append(Source(
    'A1', r'Two Wilson lines, first term', 'fourrho', -I, ("w'", 'w'), r'K_1', RL, K1, 'plain', False,
    [(+1, [('f', 'a', 'b', 'c'), ('U', 'd', 'e', "w'"), ('U', 'e', 'b', 'w')], [('d', "x'"), ('a', 'x'), ('c', 'y')]),
     (-1, [('f', 'a', 'b', 'c'), ('U', 'e', 'b', 'w'), ('U', 'e', 'd', "x'")], [('d', "x'"), ('a', 'x'), ('c', 'y')])],
    RL, ["w'", 'w', 'z']))
SOURCES.append(Source(
    'A2', r'Two Wilson lines, second term', 'fourrho', -I, ("w'", 'w'), r'K_2', RL, K2, 'plain', False,
    [(+1, [('f', 'a', 'b', 'c'), ('U', 'e', 'b', "w'"), ('U', 'e', 'd', 'y')], [('c', "x'"), ('a', 'x'), ('d', 'y')]),
     (-1, [('f', 'a', 'b', 'c'), ('U', 'b', 'e', "w'"), ('U', 'e', 'd', 'w')], [('c', "x'"), ('a', 'x'), ('d', 'y')])],
    RL, ["w'", 'w', 'z']))
SOURCES.append(Source(
    'B1', r'Three Wilson lines, kernel $K_1$', 'fourrho', I, ("w'", 'w'), r'K_1', RL, K1, 'plain', False,
    [(-1, [('f', 'a', 'b', 'c'), ('U', 'd', 'b', "w'"), ('U', 'e', 'c', 'y'), ('U', 'a', 'f', 'x')], [('d', "x'"), ('f', 'x'), ('e', 'y')]),
     (+1, [('f', 'a', 'b', 'c'), ('U', 'e', 'c', 'y'), ('U', 'b', 'd', "x'"), ('U', 'a', 'f', 'x')], [('d', "x'"), ('f', 'x'), ('e', 'y')])],
    RL, ["w'", 'w', 'z']))
SOURCES.append(Source(
    'B2', r'Three Wilson lines, kernel $K_3$', 'fourrho', I, ("w'", 'w'), r'K_3', RL, K3, 'plain', False,
    [(+1, [('f', 'a', 'b', 'c'), ('U', 'e', 'c', "x'"), ('U', 'd', 'b', 'w'), ('U', 'a', 'f', 'y')], [('e', "x'"), ('d', 'x'), ('f', 'y')]),
     (-1, [('f', 'a', 'b', 'c'), ('U', 'e', 'c', "x'"), ('U', 'b', 'd', 'x'), ('U', 'a', 'f', 'y')], [('e', "x'"), ('d', 'x'), ('f', 'y')])],
    RL, ["w'", 'w', 'z']))

# ---------------------------------------------------------------- Group I
# conjugate LO amplitude [U^{ab'}(x') - U^{ab'}(w')] rho^{b'}(x')
def lo_conj(sign_list):
    return [(+1, [('U', 'a', "b'", "x'")]), (-1, [('U', 'a', "b'", "w'")])]


def product(left, right, extra=()):
    out = []
    for s1, c1 in left:
        for s2, c2, w2 in right:
            out.append((s1 * s2, list(extra) + c1 + c2, w2))
    return out


# Row II: kernel KK^i(x'-w') F^i_A(w;x,y)
GII = [(1, [raw(r"\KK^{i}(@0-w')\,\mathfrak F^{i}(w;@1,@2)", "x'", 'x', 'y')])]
right_II = [(+1, [('f', 'h', 'b', 'c'), ('U', 'a', 'h', 'w')], [("b'", "x'"), ('b', 'x'), ('c', 'y')]),
            (-1, [('f', 'a', 'b', 'c'), ('U', 'b', 'd', 'x'), ('U', 'c', 'e', 'y')], [("b'", "x'"), ('d', 'x'), ('e', 'y')])]
SOURCES.append(Source(
    'GI.II', r'Group I, Row II (two-charge part of $\mathbb A^{(3)}$)', 'genuine', -I / 2, ("w'", 'w'),
    r'G_{\rm I.II}', RL, GII, 'A{BC}', True,
    product(lo_conj(None), right_II), RL, ["w'", 'w']))

# Row III: f^{abd}, kernel K_{I.III}
KIII = [(1, [KK("x'", "w'", 'x', 'w'), KK('x', 'z', 'y', 'z')]),
        (2, [KK("x'", "w'", 'y', 'z'), raw(r'\Phi(@0)', 'x')])]
right_III = [(+1, [('U', 'd', 'e', 'x'), ('U', 'b', 'c', 'z')], [("b'", "x'"), ('e', 'x'), ('c', 'y')]),
             (-1, [('U', 'd', 'e', 'x'), ('U', 'b', 'c', 'y')], [("b'", "x'"), ('e', 'x'), ('c', 'y')])]
SOURCES.append(Source(
    'GI.III', r'Group I, Row III ($\bar{\mathbb B}_3^\dagger\mathbb A$)', 'genuine', I, ("w'", 'w'),
    r'K_{\rm I.III}', RL, KIII, 'plain', True,
    product(lo_conj(None), right_III, extra=[('f', 'a', 'b', 'd')]), RL, ["w'", 'w', 'z']))

# Row IV: kernel K_B (2 i N)
KB = [(1, [KK('x', 'w', "x'", "w'"), KK('z', 'w', 'y', 'z')]),
      (R(1, 2), [KK('x', 'w', "x'", "w'"), KK('x', 'z', 'y', 'z')])]
right_IV = [(+1, [('f', 'c', 'd', "c'"), ('U', 'a', 'c', 'w'), ('U', 'b', 'd', 'z'), ('U', 'b', 'e', 'y')],
             [("b'", "x'"), ('e', 'y'), ("c'", 'x')]),
            (-1, [('f', 'a', 'b', 'c'), ('U', 'b', 'e', 'y'), ('U', 'c', 'd', 'x')],
             [("b'", "x'"), ('e', 'y'), ('d', 'x')])]
SOURCES.append(Source(
    'GI.IV', r'Group I, Row IV ($\bar{\mathbb A}^\dagger\mathbb B_2$)', 'genuine', 2 * I, ("w'", 'w'),
    r'K_{\mathfrak B}', RL, KB, 'plain', True,
    product(lo_conj(None), right_IV), RL, ["w'", 'w', 'z']))

# Row V: two kernels, rho positions x', y, y'
RLV = ["x'", 'y', "y'"]
KV1 = [(1, [KK("x'", "w'", 'y', 'z'), KK('z', 'x', "y'", 'x')])]
KV2 = [(1, [KK("x'", "w'", "y'", 'x'), KK('z', 'x', 'y', 'z')])]
right_V = [(+1, [('f', 'a', 'd', 'c'), ('U', 'd', 'b', 'x'), ('U', 'c', 'e', 'z')], [("b'", "x'"), ('e', 'y'), ('b', "y'")]),
           (-1, [('f', 'a', 'b', 'c'), ('U', 'c', 'd', 'y'), ('U', 'b', 'e', "y'")], [("b'", "x'"), ('d', 'y'), ('e', "y'")])]
SOURCES.append(Source(
    'GI.V1', r"Group I, Row V, first term (measured gluon at $z$)", 'genuine', -I / 2, ("w'", 'z'),
    r'K_{\rm V1}', RLV, KV1, 'A{BC}', True,
    product(lo_conj(None), right_V), RLV, ["w'", 'x', 'z']))
SOURCES.append(Source(
    'GI.V2', r"Group I, Row V, second term (measured gluon at $x$)", 'genuine', -I / 2, ("w'", 'x'),
    r'K_{\rm V2}', RLV, KV2, 'A{BC}', True,
    product(lo_conj(None), right_V), RLV, ["w'", 'x', 'z']))

# ---------------------------------------------------------------- Group II, Row I
left_G1 = [(+1, [('U', 'a', "c'", "w'")], [("c'", "x'"), ('b', 'y')]),
           (-1, [('U', 'b', "c'", 'z'), ('U', 'a', "d'", "x'"), ('U', "c'", "e'", 'y')], [("d'", "x'"), ("e'", 'y')])]
right_G1 = [(+1, [('f', 'c', 'b', 'd'), ('U', 'a', 'c', 'w')], [('d', 'x')]),
            (-1, [('f', 'a', 'c', 'd'), ('U', 'b', 'c', 'z'), ('U', 'd', 'e', 'x')], [('e', 'x')])]
sub = []
for s1, c1, w1 in left_G1:
    for s2, c2, w2 in right_G1:
        sub.append((s1 * s2, c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GII.I1', r'Group II, Row I, first term', 'genuine', -I, ("w'", 'w'),
    r'K_{\mathfrak B}', RL, KB, '{BC}A', False, sub, RL, ["w'", 'w', 'z']))

KBbar = [(1, [KK('x', "w'", "x'", 'w'), KK('z', "w'", 'y', 'z')]),
         (R(1, 2), [KK('x', "w'", "x'", 'w'), KK('x', 'z', 'y', 'z')])]
left_G2 = [(+1, [('f', "c'", 'b', "d'"), ('U', 'a', "c'", "w'")], [("d'", 'x')]),
           (-1, [('f', 'a', "c'", "d'"), ('U', 'b', "c'", 'z'), ('U', "d'", "e'", 'x')], [("e'", 'x')])]
right_G2 = [(+1, [('U', 'a', 'c', 'w')], [('c', "x'"), ('b', 'y')]),
            (-1, [('U', 'b', 'c', 'z'), ('U', 'a', 'd', "x'"), ('U', 'c', 'e', 'y')], [('d', "x'"), ('e', 'y')])]
sub = []
for s1, c1, w1 in left_G2:
    for s2, c2, w2 in right_G2:
        sub.append((s1 * s2, c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GII.I2', r'Group II, Row I, second term', 'genuine', I, ("w'", 'w'),
    r'\bar K_{\mathfrak B}', RL, KBbar, 'A{BC}', False, sub, RL, ["w'", 'w', 'z']))

# ---------------------------------------------------------------- Group III
right_G3 = [(+1, [('f', 'c', 'b', 'd'), ('U', 'a', 'c', 'w')], [('d', 'x')]),
            (-1, [('f', 'a', 'c', 'd'), ('U', 'b', 'c', 'z'), ('U', 'd', 'e', 'x')], [('e', 'x')])]
left_H = [(+1, [('U', 'a', "b'", "x'"), ('U', "c'", "d'", 'y'), ('U', "c'", 'b', 'z')], [("b'", "x'"), ("d'", 'y')]),
          (-1, [('U', 'b', "d'", 'z'), ('U', "d'", "e'", 'y'), ('U', 'a', "c'", "w'")], [("c'", "x'"), ("e'", 'y')])]
sub = []
for s1, c1, w1 in left_H:
    for s2, c2, w2 in right_G3:
        sub.append((s1 * s2, c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GIII.I', r'Group III, Row I', 'genuine', -2 * I, ("w'", 'w'),
    r'K_{\mathfrak B}', RL, KB, 'plain', True, sub, RL, ["w'", 'w', 'z']))

left_I = [(+1, [('U', "c'", "d'", 'y'), ('U', 'a', "b'", "x'"), ('U', 'b', "c'", 'z')], [("d'", 'y'), ("b'", "x'")]),
          (-1, [('U', 'a', "b'", "x'")], [('b', 'y'), ("b'", "x'")])]
sub = []
for s1, c1, w1 in left_I:
    for s2, c2, w2 in right_G3:
        sub.append((s1 * s2, c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GIII.II', r'Group III, Row II', 'genuine', -2 * I, ("w'", 'w'),
    r'K_{\mathfrak B}', RL, KB, 'plain', True, sub, RL, ["w'", 'w', 'z']))

KC = [(1, [KK("y'", 'z', 'x', 'z'), KK("x'", "y'", 'y', 'w')])]
left_J = [(+1, [('U', "e'", "c'", "y'"), ('U', "d'", 'b', 'z')], [("e'", "x'")]),
          (-1, [('U', 'b', "d'", 'z'), ('U', "c'", "e'", "x'")], [("e'", "x'")])]
right_J = [(+1, [('U', 'a', 'c', 'w')], [('c', 'y'), ('b', 'x')]),
           (-1, [('U', 'b', 'c', 'z'), ('U', 'a', 'd', 'y'), ('U', 'c', 'e', 'x')], [('d', 'y'), ('e', 'x')])]
sub = []
for s1, c1, w1 in left_J:
    for s2, c2, w2 in right_J:
        sub.append((s1 * s2, [('f', "c'", "d'", 'a')] + c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GIII.III', r'Group III, Row III, first part', 'genuine', -I, ("y'", 'w'),
    r'K_{\mathfrak C}', RL, KC, 'A{BC}', True, sub, RL, ["y'", 'w', 'z']))

left_K = [(+1, [('U', "e'", "c'", "y'")], [("e'", "x'")]),
          (-1, [('U', "c'", "e'", "x'")], [("e'", "x'")])]
right_K = [(+1, [('U', 'c', 'e', 'x'), ('U', 'a', 'd', 'y')], [('e', 'x'), ('d', 'y')]),
           (-1, [('U', 'a', 'd', 'w'), ('U', 'c', 'e', 'x')], [('e', 'x'), ('d', 'y')])]
sub = []
for s1, c1, w1 in left_K:
    for s2, c2, w2 in right_K:
        sub.append((s1 * s2, [('f', "c'", 'c', 'a')] + c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GIII.V', r'Group III, Row V', 'genuine', -2 * I, ("y'", 'w'),
    r'K_{\mathfrak C}', RL, KC, 'plain', True, sub, RL, ["y'", 'w', 'z']))

left_L = left_J
right_L = [(+1, [('U', 'b', 'c', 'z'), ('U', 'a', 'd', 'y'), ('U', 'c', 'e', 'x')], [('d', 'y'), ('e', 'x')]),
           (-1, [('U', 'a', 'd', 'y')], [('d', 'y'), ('b', 'x')])]
sub = []
for s1, c1, w1 in left_L:
    for s2, c2, w2 in right_L:
        sub.append((s1 * s2, [('f', "c'", "d'", 'a')] + c1 + c2, w1 + w2))
SOURCES.append(Source(
    'GIII.VI', r'Group III, Row VI', 'genuine', -2 * I, ("y'", 'w'),
    r'K_{\mathfrak C}', RL, KC, 'plain', True, sub, RL, ["y'", 'w', 'z']))
