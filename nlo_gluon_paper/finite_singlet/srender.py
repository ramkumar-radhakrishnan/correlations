"""LaTeX of the singlet colour factors (project.pkl)."""
import sympy as sp
from trace import Nc

PT = {"x'": "x'", "y'": "y'", "w'": "w'", 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w'}


def letter(l):
    if l[0] == 'U':
        return 'U(%s)' % PT[l[1]]
    if l[0] == 'Ud':
        return r'U^\dagger(%s)' % PT[l[1]]
    if l[0] == 'D':
        return r'\big(U(%s)-U(%s)\big)' % (PT[l[1]], PT[l[2]])
    if l[0] == 'Dd':
        return r'\big(U^\dagger(%s)-U^\dagger(%s)\big)' % (PT[l[1]], PT[l[2]])
    return 'T^e'


def word(w):
    return r'{\rm Tr}\big\{' + ''.join(letter(l) for l in w) + r'\big\}'


def _poly(s):
    """polynomial in N_c with rational coefficients -> (sign, latex magnitude); '' for 1"""
    s = sp.expand(s)
    terms = sp.Add.make_args(s)
    if len(terms) == 1:
        c, rest = s.as_coeff_Mul()
        sign = '-' if c < 0 else '+'
        c = abs(c)
        n = sp.degree(rest, Nc) if rest != 1 else 0
        nc = '' if n == 0 else ('N_c' if n == 1 else 'N_c^{%d}' % n)
        if c == 1:
            mag = nc
        elif c.q == 1:
            mag = '%d%s' % (c.p, nc)
        else:
            mag = r'\tfrac{%d}{%d}%s' % (c.p, c.q, nc)
        return sign, mag
    f = sp.factor(s)
    if f.could_extract_minus_sign():
        return '-', '(%s)' % sp.latex(sp.expand(-f)).replace('N_{c}', 'N_c')
    return '+', '(%s)' % sp.latex(s).replace('N_{c}', 'N_c')


def coef(c, first):
    c = sp.expand(c)
    re, im = sp.re(c), sp.im(c)
    if im == 0:
        sign, mag = _poly(re)
        iu = ''
    else:
        assert re == 0, c
        sign, mag = _poly(im)
        iu = 'i'
    if mag.startswith(r'\tfrac{1}{'):
        mag = r'\tfrac{%s}{%s' % (iu or '1', mag[len(r'\tfrac{1}{'):]) if iu else mag
        iu = ''
    out = iu + (r'\,' if iu and mag else '') + mag
    if sign == '+':
        sign = '' if first else '+'
    return sign + out


def mono(ws, c, first):
    return coef(c, first) + r'\,' + r'\,'.join(word(w) for w in ws)


def bracket(d, per_line=2, indent=r'&\qquad'):
    """d: {words: coef} -> lines 'a + b \\ + c' (no outer brackets)"""
    items = sorted(d.items(), key=lambda it: (sum(len(w) for w in it[0]), str(it[0])))
    lines, cur, n = [], [], 0
    for k, (ws, c) in enumerate(items):
        t = mono(ws, c, k == 0)
        width = sum(len(w) + (2 if any(l[0] in ('D', 'Dd') for l in w) else 0) for w in ws)
        if cur and (n + width > 9):
            lines.append(''.join(cur))
            cur, n = [], 0
        cur.append(t)
        n += width
    if cur:
        lines.append(''.join(cur))
    return ('\\notag\\\\\n' + indent).join(lines)
