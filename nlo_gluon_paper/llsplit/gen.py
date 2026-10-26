r"""Split every term of your final row expressions into
    J  (the coefficient of log(k+/Lambda), written with log(k+/Lambda)),
    B  (the coefficient of log(vee/k+), written with log(vee/k+)),
    F  (the finite remainder, Lambda -> 0 where possible),
exactly, term by term, using logsplit.py.  A term carries at most one logarithmic factor (a single log or one of
the brackets [log(k+/Lambda) + log(vee/k+ - 1)], [log(vee/k+) + log(k+/Lambda)]); the few terms with a log inside a sum
with a constant are treated in SPECIAL.  check(): for every evaluable term, original = J + B + F_exact numerically."""
import re

import blocks as BL
import logsplit as LS
import texeval as TE

LKL, LVK = LS.LKL, LS.LVK


def rx(s):
    """regex for the LaTeX string s, insensitive to spaces and \\,"""
    toks = re.findall(r'\\[A-Za-z]+|\\.|.', LS.norm(s))
    return re.compile(r'(?:\s|\\,)*'.join(re.escape(t) for t in toks))


def strip_user_colours(t):
    return TE.unwrap_colours(t, ('blue', 'green', 'orange'))


# bracket units (one logarithmic factor)
L1 = [r"\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+\log{\left(\frac{\vee}{k^{+}}-1\right)}\right]",
      r"\left[\log{\left(\frac{\vee}{k^{+}}-1\right)}+\log{\left(\frac{k^{+}}{\Lambda}\right)}\right]"]
L2 = [r"\left[\log{\left[\frac{\vee}{k^{+}}\right]}+\log{\left[\frac{k^{+}}{\Lambda}\right]}\right]",
      r"\left[\log{\left[\frac{\vee}{k^{+}}\right]}+\log{\left[\frac{k^{+}}{\Lambda}\right]}\right]"]
UNITS = [(u, (1, 1, LS.L1MA, LS.L1MA)) for u in L1] + [(u, (1, 1, None, None)) for u in L2]
UNITS += [(s, v) for s, v in []]
UNIT_RX = [(rx(u), v, u) for u, v in UNITS]
LOG_RX = None


def coef(c, L):
    if c == 1:
        return L
    if c == -1:
        return r'\left(-%s\right)' % L
    return r'%d\,%s' % (c, L)


def find_unit(t):
    """the logarithmic factor of a term: (start, end, (cL, cV, R, Rx), text) or None"""
    for r, v, u in UNIT_RX:
        m = r.search(t)
        if m:
            # strip an enclosing {..} left over from the colour command
            return m.start(), m.end(), v, u
    found = []
    for lg in BL.logs_in(t):
        k = LS.norm(lg)
        if k in LS.SPLIT:
            found.append(lg)
        elif k not in LS.FINITE and not re.match(r'\\log\\frac\{\((w|y)-|\\log\\frac\{\(w', k):
            raise KeyError('unknown log: %s' % lg[:120])
    big = [lg for lg in found if LS.SPLIT[LS.norm(lg)][:2] != (0, 0)]
    small = [lg for lg in found if LS.SPLIT[LS.norm(lg)][:2] == (0, 0)]
    if len(big) > 1:
        raise ValueError('two large logs in one term: %s' % big)
    lg = (big or small or [None])[0]
    if lg is None:
        return None
    i = t.index(lg)
    return i, i + len(lg), LS.SPLIT[LS.norm(lg)], lg


def split_term(t, hidden=False):
    """-> dict(J=..., B=..., F=..., Fx=...) LaTeX of the term (None = no such part)"""
    t = strip_user_colours(t)
    u = find_unit(t)
    if u is None:
        return dict(J=None, B=None, F=t, Fx=t)
    i, j, (cL, cV, R, Rx), txt = u
    out = {}
    out['J'] = t[:i] + coef(cL, LKL) + t[j:] if cL else None
    cVe = 0 if hidden else cV
    out['B'] = t[:i] + coef(cVe, LVK) + t[j:] if cVe else None
    if hidden and txt in L1:
        R = Rx = LS.LVK1
    out['F'] = t[:i] + R + t[j:] if R is not None else None
    out['Fx'] = t[:i] + Rx + t[j:] if Rx is not None else None
    # other small logs in the term: Lambda -> 0 in F
    if out['F'] is not None:
        for lg in BL.logs_in(out['F']):
            k = LS.norm(lg)
            if k in LS.SPLIT and LS.SPLIT[k][:2] == (0, 0) and LS.SPLIT[k][2] is not None:
                out['F'] = out['F'].replace(lg, LS.SPLIT[k][2])
    return out


if __name__ == '__main__':
    BLK = BL.blocks()
    for name, B in BLK.items():
        print('==', name)
        for n, t in enumerate(B.terms):
            try:
                s = split_term(t)
                print('   %2d  J:%s B:%s F:%s' % (n, s['J'] is not None, s['B'] is not None, s['F'] is not None))
            except Exception as e:
                print('   %2d  ERROR %s' % (n, str(e)[:150]))


# ---------------------------------------------------------------- special terms
HIDDEN_ROWS = {('G2 IV', 1), ('G3 III(1)', 1), ('G3 III(2)', 1), ('G3 V', 1), ('G3 VI', 1)}
POWER = {('G1 III N1', 0), ('G1 III N1', 3), ('G1 III N2', 0), ('G1 III N3a', 0)}     # 1/Lambda terms: cancel exactly
RED = {('G1 III N2', 2), ('G1 III N2', 3)}                                           # your red terms: cancel each other


def special(name, n, t):
    t = strip_user_colours(t)
    if name == 'G1 II(1)':
        i = t.index(r'\left[\frac{11}{3}')
        j = t.index(r'\Bigg]\times', i)
        J = t[:i] + r'\left(-2\,%s\right)' % LKL + t[j:]
        F = t.replace(r'\frac{11}{3}-2\,\log{\left(\frac{k^{+}}{\Lambda}\right)}', r'\frac{11}{3}')
        assert F != t
        return dict(J=J, B=None, F=F, Fx=F)
    if name == 'G1 V(2b UV)':
        m = rx(r'\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]').search(t)
        J = t[:m.start()] + r'2\,' + LKL + t[m.end():]
        F = t[:m.start()] + r'\left[-\frac{11}{6}\right]' + t[m.end():]
        Fx = t[:m.start()] + r'\left[2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]' + t[m.end():]
        return dict(J=J, B=None, F=F, Fx=Fx, E1=True)
    if name == 'G1 V(2b fin)' and n == 5:
        m = rx(r'\bigg[-\frac{11}{6}+2\log{\left[\frac{k^{+}-\Lambda}{\Lambda}\right]}\bigg]').search(t)
        J = t[:m.start()] + r'2\,' + LKL + t[m.end():]
        F = t[:m.start()] + r'\bigg[-\,\frac{11}{6}\bigg]' + t[m.end():]
        Fx = t[:m.start()] + r'\bigg[-\,\frac{11}{6}+2\,\log{\left(1-\frac{\Lambda}{k^{+}}\right)}\bigg]' + t[m.end():]
        return dict(J=J, B=None, F=F, Fx=Fx)
    return None


def split(name, n, t):
    if (name, n) in RED:
        return dict(J=None, B=None, F=None, Fx=None, red=True)
    if (name, n) in POWER:
        return dict(J=None, B=None, F=None, Fx=None, power=True)
    s = special(name, n, t)
    if s is not None:
        return s
    return split_term(t, hidden=(name, n) in HIDDEN_ROWS)


def all_splits():
    BLK = BL.blocks()
    out = {}
    for name, B in BLK.items():
        out[name] = [split(name, n, t) for n, t in enumerate(B.terms)]
    return BLK, out
