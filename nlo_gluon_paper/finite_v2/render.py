"""LaTeX of the compact two-rho remnants (compact.pkl)."""
from fractions import Fraction

PT = {"x'": "x'", "y'": "y'", "w'": "w'", 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w'}


def coef_tex(c, first):
    """c = i * r with r rational (the remnants are purely imaginary)"""
    assert abs(c.real) < 1e-12, c
    r = Fraction(c.imag).limit_denominator(64)
    sign = '-' if r < 0 else ('' if first else '+')
    r = abs(r)
    if r == 1:
        mag = 'i\\,'
    elif r.denominator == 1:
        mag = '%d\\,i\\,' % r.numerator
    else:
        mag = '\\tfrac{%s}{%d}\\,' % ('i' if r.numerator == 1 else '%d\\,i' % r.numerator, r.denominator)
    return sign + mag


def idx(*ix):
    return ''.join(ix)


def fac_tex(u):
    if u[0] == 'U':
        return 'U^{%s}(%s)' % (idx(u[1], u[2]), PT[u[3]])
    if u[0] == 'D':
        return '\\big[U^{%s}(%s)-U^{%s}(%s)\\big]' % (idx(u[1], u[2]), PT[u[3][0]], idx(u[1], u[2]), PT[u[4][0]])
    raise ValueError(u)


def term_tex(keyn, c, first):
    key, n = keyn if isinstance(keyn[1], int) else (keyn, 0)
    fs, Us, rs = key
    s = coef_tex(c, first) + ('N_c\\,' if n == 1 else ('N_c^{%d}\\,' % n if n > 1 else ''))
    s += '\\,'.join('f^{%s}' % idx(*f) for f in fs) + '\\,'
    s += '\\,'.join(fac_tex(u) for u in Us) + '\\,'
    s += '\\,'.join('\\rho^{%s}(%s)' % (i, PT[p]) for i, p in rs)
    return s


def merge_tex(pair):
    keep, gone = pair
    return '%s\\to %s' % (PT[gone], PT[keep])


def row_block(groups, kname):
    """align body: sum over merges  K|_{gone->keep} [ terms ]"""
    lines = []
    for n, (pair, terms) in enumerate(groups.items()):
        head = ('' if n == 0 else '+') + '%s\\big|_{%s}\\Big[' % (kname, merge_tex(pair))
        items = list(terms.items())
        body = []
        for m, (k, c) in enumerate(items):
            body.append(term_tex(k, c, m == 0))
        # one term per line
        chunks = list(body)
        txt = head + chunks[0]
        for ch in chunks[1:]:
            txt += '\\notag\\\\\n&\\qquad' + ch
        txt += '\\Big]'
        lines.append(txt)
    return '&' + '\\notag\\\\\n&'.join(lines)
