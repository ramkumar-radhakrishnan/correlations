r"""Pointwise numerical evaluation of your LaTeX terms, as written.

A term is translated into Python: (a-b) -> transverse vector, \cdot -> dot product, ^{2} -> square, ^{i} -> component,
\delta -> Kronecker delta, U^{ab}(p) / \rho^{a}(p) / f^{abc} / N_c -> colour tensors, \{\rho,\rho\} -> 2 rho rho
(classical charges), e^{-ik\cdot(a-b)} -> phase, \log -> log|.|, \mathcal{A} -> the arctan function of Group I Row III.
The transverse indices of a term are summed explicitly, the colour indices by einsum.  Points, Wilson lines (adjoint
SU(3), symmetric U = U^T as in your convention U^{bc} = U^{cb}) and charges come from a configuration of ../groups23/fp.py.
Terms with a remaining p+, q+, xi or momentum integral are not evaluated (TexError)."""
import itertools
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'groups23')))
import fp  # noqa: E402

F = fp.F
NCV = 3
NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
PT = {"x": 'x', "y": 'y', "z": 'z', "w": 'w', "x'": 'xp', "y'": 'yp', "w'": 'wp', "z'": 'zp'}


class TexError(Exception):
    pass


# ---------------------------------------------------------------- colour algebra
class Col:
    """sum of monomials (coef, factors); factors: ('U', i, j, p), ('r', i, p), ('f', i, j, k)"""

    def __init__(self, mons):
        self.m = mons

    def __add__(self, o):
        if isinstance(o, Col):
            return Col(self.m + o.m)
        if o == 0:
            return self
        raise TexError('colour + scalar')

    __radd__ = __add__

    def __neg__(self):
        return Col([(-c, f) for c, f in self.m])

    def __pos__(self):
        return self

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, Col):
            return Col([(c1 * c2, f1 + f2) for c1, f1 in self.m for c2, f2 in o.m])
        return Col([(c * o, f) for c, f in self.m])

    __rmul__ = __mul__

    def __truediv__(self, o):
        return Col([(c / o, f) for c, f in self.m])

    def value(self, cfg):
        tot = 0.0
        for c, facs in self.m:
            if c == 0:
                continue
            labels = {}
            ops, subs = [], []
            for fa in facs:
                if fa[0] == 'U':
                    ops.append(cfg.U[fa[3]])
                    idx = fa[1:3]
                elif fa[0] == 'r':
                    ops.append(cfg.rho[fa[2]])
                    idx = fa[1:2]
                else:
                    ops.append(F)
                    idx = fa[1:4]
                s = ''
                for i in idx:
                    if i not in labels:
                        labels[i] = 'abcdefghijklmnopqrstuvwxyz'[len(labels)]
                    s += labels[i]
                subs.append(s)
            cnt = {}
            for s in subs:
                for ch in s:
                    cnt[ch] = cnt.get(ch, 0) + 1
            if any(v != 2 for v in cnt.values()):
                raise TexError('colour index not contracted: %s %s' % (facs, cnt))
            tot += c * np.einsum(','.join(subs) + '->', *ops)
        return tot


def scalar_or_col(v, cfg):
    return v.value(cfg) if isinstance(v, Col) else v


# ---------------------------------------------------------------- translation
def match_brace(t, i, o='{', c='}'):
    d = 0
    for k in range(i, len(t)):
        if t[k] == o:
            d += 1
        elif t[k] == c:
            d -= 1
            if d == 0:
                return k + 1
    raise TexError('unbalanced')


def unwrap_colours(t, cols=('blue', 'red', 'green', 'orange', 'black')):
    """{\\color{c} X} -> X (a colour group is not a mathematical group)"""
    for col in cols:
        cmd = r'{\color{%s}' % col
        while cmd in t:
            i = t.index(cmd)
            j = match_brace(t, i)
            t = t[:i] + t[i + len(cmd):j - 1] + t[j:]
    return t


def strip_colour_cmds(t):
    t = unwrap_colours(t)
    for col in ('blue', 'red', 'green', 'orange', 'black'):
        t = t.replace(r'\color{%s}' % col, '')
    t = t.replace(r'\nc', '').replace(r'\rc', '')
    while r'\new{' in t or r'\old{' in t:
        for cmd, keep in ((r'\new{', True), (r'\old{', False)):
            if cmd in t:
                i = t.index(cmd)
                j = match_brace(t, i + len(cmd) - 1)
                t = t[:i] + (t[i + len(cmd):j - 1] if keep else '') + t[j:]
    return t


def frac(t):
    while r'\frac' in t:
        i = t.index(r'\frac')
        j = i + 5
        while t[j] == ' ':
            j += 1
        if t[j] != '{':      # \frac12 style not used
            raise TexError('frac without braces')
        a_end = match_brace(t, j)
        k = a_end
        while t[k] == ' ':
            k += 1
        b_end = match_brace(t, k)
        t = t[:i] + '((' + t[j + 1:a_end - 1] + ')/(' + t[k + 1:b_end - 1] + '))' + t[b_end:]
    return t


def remove_int(t):
    """drop \\int_{..}^{..}; refuse integrals over p+, q+, xi, momenta"""
    if re.search(r'dp\^\{\+\}|dq\^\{\+\}|d\\xi|d\^\{2\}p', t):
        raise TexError('remaining longitudinal/momentum integral')
    while r'\int' in t:
        i = t.index(r'\int')
        j = i + 4
        for _ in range(2):
            while j < len(t) and t[j] == ' ':
                j += 1
            if j < len(t) and t[j] in '_^':
                j += 1
                if t[j] == '{':
                    j = match_brace(t, j)
                else:
                    j += 1
        t = t[:i] + t[j:]
    return t


PNT = r"(x'|y'|w'|z'|x|y|z|w|w\^\{\\prime\}|w\^\\prime)"


def pname(s):
    s = s.replace(r'^{\prime}', "'").replace(r'^\prime', "'")
    return PT[s]


def idx_tokens(s):
    s = s.replace(r'\,', '').replace(' ', '').replace(r'\dagger', '')
    return re.findall(r"[a-z]'?", s)


def translate(t):
    t = strip_colour_cmds(t)
    # \mathcal{K}(a-b) = (a-b)/(a-b)^2 (notation of 4rho.tex and of the evolution lists)
    t = re.sub(r'\\mathcal\{K\}\s*\(([^()]*)\)', lambda m: r'\frac{(%s)}{(%s)^{2}}' % (m.group(1), m.group(1)), t)
    t = t.replace(r'\frac{i}{', r'1j\frac{1}{')
    t = re.sub(r'\\(notag|nonumber|quad|qquad)', ' ', t)
    t = t.replace('\\\\', ' ').replace('&', ' ')
    t = re.sub(r'\\hspace\{[^}]*\}', ' ', t)
    t = re.sub(r'\\label\{[^}]*\}', ' ', t)
    for s in (r'\,', r'\;', r'\!', r'\:'):
        t = t.replace(s, ' ')
    t = remove_int(t)
    # anticommutator {rho, rho}
    t = re.sub(r'(\\(left|Big|big|bigg|Bigg))?\\\{\s*(?=\\rho)', ' ACOM(', t)
    t = re.sub(r'(\\rho\^\{[^}]*\}\(\s*[^)]*\))\s*(\\(right|Big|big|bigg|Bigg))?\\\}', r'\1)', t)
    t = re.sub(r'\\(left|right|Bigg|bigg|Big|big)\s*(?=[\[\]\(\)]|\\\{|\\\})', '', t)
    t = t.replace(r'\{', '(').replace(r'\}', ')')
    t = re.sub(r'\\mathbf\{(\w)\}', r'\1', t)
    # phases
    def phase(m):
        inner = m.group(1)
        inner = re.sub(r'ik\s*(\\cdot)?\s*\(' + PNT + r'\s*-\s*' + PNT + r'\s*\)',
                       lambda q: '1j*KD("%s","%s")' % (pname(q.group(2)), pname(q.group(3))), inner)
        inner = re.sub(r'ik\s*(\\cdot)?\s*\\left\(\s*' + PNT + r'\s*-\s*' + PNT + r'\s*\\right\)',
                       lambda q: '1j*KD("%s","%s")' % (pname(q.group(2)), pname(q.group(3))), inner)
        return ' EXPI(' + inner + ') '
    while 'e^{' in t:
        i = t.index('e^{')
        j = match_brace(t, i + 2)
        t = t[:i] + phase(re.match(r'(.*)', t[i + 3:j - 1], re.S)) + t[j:]
    # colour tensors
    t = re.sub(r"U\^\{([^}]*)\}\s*\(\s*" + PNT + r"\s*\)",
               lambda m: ' CU("%s","%s","%s") ' % (*(idx_tokens(m.group(1)) if r'\dagger' not in m.group(1)
                                                     else idx_tokens(m.group(1))[::-1]), pname(m.group(2))), t)
    t = re.sub(r"\\rho\^\{([^}]*)\}\s*\(\s*" + PNT + r"\s*\)",
               lambda m: ' CR("%s","%s") ' % (idx_tokens(m.group(1))[0], pname(m.group(2))), t)
    t = re.sub(r"f\^\{([^}]*)\}", lambda m: ' CF("%s","%s","%s") ' % tuple(idx_tokens(m.group(1))), t)
    t = t.replace('N_{c}', ' NCV ')
    # transverse vectors and indices
    t = re.sub(r'\(\s*' + PNT + r'\s*-\s*' + PNT + r'\s*\)', lambda m: ' D("%s","%s") ' % (pname(m.group(1)), pname(m.group(2))), t)
    t = re.sub(r"\\delta[_^]\{\s*(\w'?)\s*(\w'?)\s*\}", lambda m: ' DEL("%s","%s") ' % (m.group(1), m.group(2)), t)
    t = re.sub(r'(D\("[^"]*","[^"]*"\))\s*\^\{\s*(\w\'?)\s*\}', lambda m: ' COMP(%s,"%s") ' % (m.group(1), m.group(2)) if not m.group(2).isdigit() else m.group(0), t)
    # symbols
    reps = [(r'\mathcal{A}_{0}', ' AFN0 '), (r'\mathcal{A}', ' AFN '), (r'\mathrm{arc tan}', ' ATAN '), (r'\mathrm{arc\,tan}', ' ATAN '), (r'\sqrt', ' SQRT '),
            (r'\log^{2}', ' LOG2 '), (r'\log', ' LOG '), ('k^{+ 2}', ' KP**2 '), ('k^{+2}', ' KP**2 '),
            ('k^{+\\,2}', ' KP**2 '), ('k^{+}', ' KP '), (r'\vee', ' V '), (r'\Lambda', ' L '), (r'\epsilon', ' EPS '),
            (r'\gamma_{E}', ' GE '), (r'\gamma', ' GE '), (r'\mu^{2}_{\overline{MS}}', ' MU2 '), (r'\pi', ' PI '),
            ('d_{\\perp}', ' DP '), ('ig^{4}', ' 1j*G**4 '), ('i g^{4}', ' 1j*G**4 '), ('g^{4}', ' G**4 '),
            ('g^{2}', ' G**2 '), (r'\cdot', ' @ '), (r'\times', ' * ')]
    for a, b in reps:
        t = t.replace(a, b)
    t = re.sub(r'(?<![A-Za-z_])i(?=\s*(\*|1j|G\*\*|\(|CF))', ' 1j* ', t)
    t = frac(t)
    t = t.replace('[', '(').replace(']', ')').replace('{', '(').replace('}', ')')
    t = t.replace('^', '**')
    # implicit multiplication
    toks = re.findall(r'1j|\d+\.\d+|\d+|[A-Za-z_][A-Za-z_0-9]*|"[^"]*"|\*\*|[-+*/@(),]', t)
    funcs = {'LOG', 'LOG2', 'D', 'CU', 'CR', 'CF', 'ACOM', 'DEL', 'EXPI', 'KD', 'SQRT', 'ATAN', 'IX', 'COMP'}
    out = []
    for tk in toks:
        if out:
            prev = out[-1]
            prev_end = prev not in ('(', ',', '-', '+', '*', '/', '@', '**') and prev not in funcs
            cur_start = tk not in (')', ',', '-', '+', '*', '/', '@', '**')
            if prev_end and cur_start and not (tk == '(' and prev in funcs) and not prev.startswith('"'):
                out.append('*')
        out.append(tk)
    s = ''.join(out)
    s = s.replace(')(', ')*(').replace('IX*(', 'IX(')
    s = re.sub(r'\)\*\*2', ')**2', s)
    return s


# ---------------------------------------------------------------- evaluation
class Vec(np.ndarray):
    """transverse vector: v**2 = v.v, v**4 = (v.v)^2"""
    def __pow__(self, n):
        d = float(np.asarray(self) @ np.asarray(self))
        if n == 2:
            return d
        if n == 4:
            return d * d
        raise TexError('vector power %s' % n)

    def __matmul__(self, o):
        return float(np.asarray(self) @ np.asarray(o))

class Env:
    def __init__(self, cfg, KP=1.0, V=1e3, L=1e-3, EPS=0.37, MU2=1.7):
        self.cfg, self.KP, self.V, self.L, self.EPS, self.MU2 = cfg, KP, V, L, EPS, MU2

    def names(self, ix):
        c = self.cfg

        def D(a, b):
            return (c.pt[a] - c.pt[b]).view(Vec)

        def AFN(L=None):
            XW, XZ = D('x', 'w'), D('x', 'z')
            s = np.sqrt((XW @ XW) * (XZ @ XZ) - (XW @ XZ)**2)
            kp, L = self.KP, (self.L if L is None else L)
            return (np.arctan(((kp - L) * (XZ @ XZ) - kp * (XW @ XZ)) / (kp * s))
                    - np.arctan((L * (XZ @ XZ) - kp * (XW @ XZ)) / (kp * s))) / s
        return dict(
            D=D, COMP=lambda v, n: float(v[ix[n]]), KD=lambda a, b: c.kv @ (c.pt[a] - c.pt[b]), EXPI=np.exp, IX=lambda n: ix[n],
            DEL=lambda a, b: 1.0 if ix[a] == ix[b] else 0.0,
            CU=lambda i, j, p: Col([(1.0, [('U', i, j, p)])]), CR=lambda i, p: Col([(1.0, [('r', i, p)])]),
            CF=lambda i, j, k: Col([(1.0, [('f', i, j, k)])]), ACOM=lambda a, b: 2 * a * b,
            LOG=lambda v: np.log(np.abs(v)), LOG2=lambda v: np.log(np.abs(v))**2, SQRT=np.sqrt, ATAN=np.arctan,
            AFN=AFN(), AFN0=AFN(0.0), NCV=NCV, KP=self.KP, V=self.V, L=self.L, EPS=self.EPS, GE=np.euler_gamma, MU2=self.MU2,
            PI=np.pi, G=fp.G, DP=2.0 - self.EPS)


def term_value(tex, env):
    code = translate(tex)
    letters = sorted(set(re.findall(r'COMP\(D\("[^"]*","[^"]*"\),"([^"]*)"\)', code)) | set(sum((list(m) for m in re.findall(r'DEL\("([^"]*)","([^"]*)"\)', code)), [])))
    tot = 0.0
    for vals in itertools.product((0, 1), repeat=len(letters)):
        ix = dict(zip(letters, vals))
        try:
            v = eval(code, {'np': np}, env.names(ix))
        except TexError:
            raise
        except Exception as e:
            raise TexError('%s in %s' % (e, code[:300]))
        tot += scalar_or_col(v, env.cfg)
    return tot


def cfg(seed):
    return fp.Cfg(seed, sym=True, names=NAMES)


if __name__ == '__main__':
    print(translate(r"\frac{(x'-w')^{i}\,(y-z)^{j}}{(x'-w')^{2}\,(y-z)^{2}}\,\delta^{ij}\,e^{-ik(w^\prime-w)} \left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\,\rho^{b'}(x')\right]f^{abc}\{\rho^{c}(y),\rho^{b}(x)\}"))
