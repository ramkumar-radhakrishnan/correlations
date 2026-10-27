r"""Your BFKL.tex and Other_terms.tex (input/), cut into blocks and terms with the machinery of ../llsplit.

UB[name]: the blocks of BFKL.tex (terms you already wrote with log(vee/k+)), UO[name]: the blocks of Other_terms.tex.
The names are those of ../llsplit (row of Group_I/II/III.tex)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LL = os.path.normpath(os.path.join(HERE, '..', 'llsplit'))
if LL not in sys.path:
    sys.path.insert(0, LL)
import blocks as BL  # noqa: E402

FILES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n') for g in ('BFKL', 'Other_terms')}

B_RANGES = {
    '4rho a': (34, 39), '4rho b': (40, 45), '4rho 1': (47, 68), '4rho 2': (69, 101),
    'G1 II(2a)': (106, 110), 'G1 III N1': (112, 114), 'G1 III N2': (115, 117), 'G1 IV': (119, 121),
    'G1 V(2b UV)': (126, 130),
    'G2 I(2d)': (136, 138), 'G2 I(2e)': (139, 141), 'G2 I(2f)': (142, 144), 'G2 IV': (150, 155),
    'G3 I': (158, 160), 'G3 II': (162, 164)}
O_RANGES = {
    'G1 II(1)': (40, 42), 'G1 II(2b)': (43, 48),
    'G1 III N1': (57, 59), 'G1 III N2': (60, 62), 'G1 III N3a': (63, 78), 'G1 III N3b': (79, 87),
    'G1 IV': (93, 96), 'G1 V(1)': (99, 101), 'G1 V(2a)': (103, 105), 'G1 V(2b UV)': (106, 110),
    'G1 V(2b fin)': (111, 118),
    'G2 I(1)': (122, 124), 'G2 I(2a)': (127, 129), 'G2 I(2b)': (130, 132), 'G2 I(2c)': (133, 135),
    'G2 I(2d)': (136, 138), 'G2 I(2e)': (139, 141), 'G2 I(2f)': (142, 144), 'G2 IV': (150, 152),
    'G3 I': (155, 157), 'G3 II': (159, 161), 'G3 III(1)': (164, 166), 'G3 III(2)': (168, 172),
    'G3 V': (176, 178), 'G3 VI': (180, 182)}


def align_text(g, a, b):
    t = '\n'.join(FILES[g][a - 1:b])
    i = t.index(r'\begin{align}') + len(r'\begin{align}')
    j = t.rindex(r'\end{align}')
    return t[i:j]


def g2_i1_parts(t):
    import re
    m = re.search(r"&\s*\+\s*(?=\\frac\{1\}\{\(2\\pi\)\^\{3\}\}\\,\\int_\{x,y,z,w,x',w\^\{\\prime\}\})", t)
    a = t[:m.start()]
    b = r'\frac{d^{3}N}{d^{3}k} & = ' + t[m.end():]
    return a, b


def user_blocks():
    UB = {n: BL.Block('B ' + n, align_text('BFKL', a, b)) for n, (a, b) in B_RANGES.items()}
    UO = {}
    for n, (a, b) in O_RANGES.items():
        t = align_text('Other_terms', a, b)
        if n == 'G2 I(1)':
            A, Bt = g2_i1_parts(t)
            UO['G2 I(1) A'] = BL.Block('O G2 I(1) A', A)
            UO['G2 I(1) B'] = BL.Block('O G2 I(1) B', Bt)
        else:
            UO[n] = BL.Block('O ' + n, t)
    return UB, UO


if __name__ == '__main__':
    import re
    UB, UO = user_blocks()
    for tag, D in (('BFKL.tex', UB), ('Other_terms.tex', UO)):
        print('=====', tag)
        for n, B in D.items():
            print('%-14s braced %-5s terms %d' % (n, B.braced, len(B.terms)))
            for k, t in enumerate(B.terms):
                print('     %2d %s' % (k, re.sub(r'\s+', ' ', BL.cell_text(t))[:110]))
