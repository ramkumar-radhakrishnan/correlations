r"""Pointwise evaluation of your equations with the longitudinal variable kept: p^+, q^+ or xi as a numpy array.

Built on ../llsplit/texeval.py (transverse vectors, colour tensors, classical charges, U = U^T).  Added here:
  * p^{+}, q^{+}, \xi            -> arrays PP, QQ, XI (one evaluation gives the integrand on all quadrature nodes);
  * every phase e^{...}           -> its value, with the longitudinal factor, e.g. e^{-ik\cdot\frac{p^+}{k^+}(y'-z)};
  * [1/p^+]_+, [k^+/p^+]_+, ...   -> a marker; the integrand is A + B * marker, and the plus prescription is
                                     A(p) + [B(p) - B(0)] * marker(p)   (your definition: subtract the value at p = 0);
  * \int_a^b dp^+/(2pi) (or dq^+, dxi) is removed from the text and done numerically (composite Gauss-Legendre in log p).
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'llsplit'))
import texeval as TE  # noqa: E402

LINES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n') for g in ('Group_I', 'Group_II', 'Group_III')}


# ---------------------------------------------------------------- array-safe colour algebra
def _add(self, o):
    if isinstance(o, TE.Col):
        return TE.Col(self.m + o.m)
    if np.isscalar(o) and o == 0:
        return self
    raise TE.TexError('colour + scalar')


def _value(self, cfg):
    tot = 0.0
    for c, facs in self.m:
        if np.isscalar(c) and c == 0:
            continue
        labels, ops, subs = {}, [], []
        for fa in facs:
            if fa[0] == 'U':
                ops.append(cfg.U[fa[3]])
                idx = fa[1:3]
            elif fa[0] == 'r':
                ops.append(cfg.rho[fa[2]])
                idx = fa[1:2]
            else:
                ops.append(TE.F)
                idx = fa[1:4]
            s = ''
            for i in idx:
                if i not in labels:
                    labels[i] = 'abcdefghijklmnopqrstuvwxyz'[len(labels)]
                s += labels[i]
            subs.append(s)
        tot = tot + c * np.einsum(','.join(subs) + '->', *ops)
    return tot


TE.Col.__add__ = _add
TE.Col.__radd__ = _add
TE.Col.value = _value
TE.Col.__array_ufunc__ = None          # ndarray * Col -> Col.__rmul__ (array coefficients)


# ---------------------------------------------------------------- text helpers
def align_rhs(group, a, b):
    """right-hand side of the align environment in lines a..b (1-based), after the first '& ='"""
    t = '\n'.join(LINES[group][a - 1:b])
    i = t.index(r'\begin{align}') + len(r'\begin{align}')
    j = t.rindex(r'\end{align}')
    t = t[i:j]
    t = re.sub(r'(?<!\\)%[^\n]*', '', t)
    k = re.search(r'&\s*=', t)
    return t[k.end():]


def _script(t, j):
    """argument of _ or ^ at t[j] (after the symbol): text, end index"""
    while t[j] == ' ':
        j += 1
    if t[j] == '{':
        je = TE.match_brace(t, j)
        return t[j + 1:je - 1], je
    return t[j], j + 1


def strip_long_int(t, var=r'p\^\{\+\}'):
    r"""remove every \int_{a}^{b} (or \int^{b}_{a}) \frac{dp^{+}}{2\pi} (resp. dq^+, d\xi); return text, [(a, b, measure)]"""
    lims = []
    pos = 0
    while True:
        m = re.compile(r'\\int\s*(?=[_^])').search(t, pos)
        if not m:
            break
        i, j = m.start(), m.end()
        first = t[j]
        a1, j = _script(t, j + 1)
        while t[j] == ' ':
            j += 1
        if t[j] not in '_^' or t[j] == first:
            pos = j
            continue
        a2, ke = _script(t, j + 1)
        lo, hi = (a1, a2) if first == '_' else (a2, a1)
        rest = t[ke:]
        mm = re.match(r'(\s|\\,)*(\\frac\{d%s\}\{2\\pi\}|\\frac\{d%s\}\{%s\}|d%s)' % (var, var, r'\\xi\^\{2\}', var), rest)
        if not mm:
            # the measure may come later in the term (e.g. \int_0^{k^+} \int_{x,..} e^{..} \frac{dq^+}{2\pi})
            far = list(re.finditer(r'\\frac\{d%s\}\{2\\pi\}' % var, rest))
            if len(far) == 1 and not re.match(r'(\s|\\,)*\\frac\{d', rest):
                f = far[0]
                lims.append((lo, hi, f.group(0)))
                t = t[:i] + ' ' + rest[:f.start()] + ' ' + rest[f.end():]
                pos = i
                continue
            pos = ke
            continue
        lims.append((lo, hi, mm.group(2)))
        t = t[:i] + ' ' + rest[mm.end():]
        pos = i
    return t, lims


LONG = [(r'k^{+\,3}', '{KP}^{3}'), (r'k^{+3}', '{KP}^{3}'), (r'p^{+\,2}', '{PP}^{2}'), (r'p^{+2}', '{PP}^{2}'), (r'p^{+}', ' PP '), (r'q^{+}', ' QQ '), (r'\xi', ' XI ')]


def longvars(t):
    for a, b in LONG:
        t = t.replace(a, b)
    return t


PLUS = [(r'\left[\frac{1}{p^{+}}\right]_{+}', r'\frac{1}{p^{+}}'), (r'\left[\frac{k^{+}}{p^{+}}\right]_{+}', r'\frac{k^{+}}{p^{+}}'),
        (r'\left[\frac{1}{q^{+}}\right]_{+}', r'\frac{1}{q^{+}}'), (r'\left[\frac{k^{+}}{q^{+}}\right]_{+}', r'\frac{k^{+}}{q^{+}}'),
        (r'\frac{\xi}{(1-\xi)_{+}}', r'\frac{\xi}{(1-\xi)}')]


def mark_plus(t):
    """replace the plus distribution [m(x)]_+ by the marker PLS; return text and m as LaTeX (or None)"""
    found = None
    for p, m in PLUS:
        if p in t:
            assert found is None or found == m, 'two different plus distributions'
            found = m
            t = t.replace(p, r' PLS ')
    return t, found


PT_RE = r"(x'|y'|w'|z'|x|y|z|w|w\^\{\\prime\}|w\^\\prime)"


class Ctx:
    """one configuration and the longitudinal values"""

    def __init__(self, cfg, KP=1.0, V=20.0, L=1e-3, EPS=0.37):
        self.cfg, self.KP, self.V, self.L, self.EPS = cfg, KP, V, L, EPS

    def env(self, extra):
        base = TE.Env(self.cfg, KP=self.KP, V=self.V, L=self.L, EPS=self.EPS)
        names = base.names

        class E:
            cfg = self.cfg

            def names(_, ix):
                d = names(ix)
                d.update(extra)
                return d
        return E()


def phase_values(t, ctx, long):
    """replace each e^{...} by a placeholder; its value includes the longitudinal factor"""
    vals = {}
    n = 0
    while 'e^{' in t:
        i = t.index('e^{')
        j = TE.match_brace(t, i + 2)
        inner = t[i + 3:j - 1]
        m = re.search(r"(\\left)?\(\s*" + PT_RE + r"\s*-\s*" + PT_RE + r"\s*(\\right)?\)", inner)
        if not m:
            raise TE.TexError('phase ' + inner)
        a, b = TE.pname(m.group(2)), TE.pname(m.group(3))
        fac = inner[:m.start()] + ' ' + inner[m.end():]
        fac = re.sub(r'^\s*-\s*ik\s*(\\cdot)?', '', fac.strip())
        fac = fac.replace(r'\,', ' ').strip()
        if fac:
            code = TE.translate(longvars(fac))
            f = eval(code, {'np': np}, dict(PP=long.get('PP'), QQ=long.get('QQ'), XI=long.get('XI'), KP=ctx.KP, PI=np.pi))
        else:
            f = 1.0
        kd = ctx.cfg.kv @ (ctx.cfg.pt[a] - ctx.cfg.pt[b])
        name = 'PHX%d' % n
        vals[name] = np.exp(-1j * f * kd)
        t = t[:i] + ' ' + name + ' ' + t[j:]
        n += 1
    return t, vals


_cache = {}


# ---------------------------------------------------------------- Einstein-correct evaluation of the translated code
import ast  # noqa: E402
import itertools  # noqa: E402


class V2:
    """transverse vector whose components may be arrays (p+ on the quadrature nodes)"""

    def __init__(self, a, b):
        self.c = (a, b)

    def __add__(self, o):
        return V2(self.c[0] + o.c[0], self.c[1] + o.c[1])

    def __sub__(self, o):
        return V2(self.c[0] - o.c[0], self.c[1] - o.c[1])

    def __neg__(self):
        return V2(-self.c[0], -self.c[1])

    def __mul__(self, o):
        if isinstance(o, (V2, TE.Col)):
            return NotImplemented
        return V2(self.c[0] * o, self.c[1] * o)

    __rmul__ = __mul__
    __array_ufunc__ = None

    def __truediv__(self, o):
        return V2(self.c[0] / o, self.c[1] / o)

    def __matmul__(self, o):
        return self.c[0] * o.c[0] + self.c[1] * o.c[1]

    def __pow__(self, n):
        d = self @ self
        if n == 2:
            return d
        if n == 4:
            return d * d
        raise TE.TexError('vector power %s' % n)

    def __getitem__(self, i):
        return self.c[i]


class T:
    """value with free transverse indices: letters (tuple) and data {assignment tuple: value}"""

    def __init__(self, letters, data):
        self.L, self.d = tuple(letters), data

    @staticmethod
    def s(v):
        return T((), {(): v})


def _combine(a, b, op, contract):
    L = list(a.L) + [x for x in b.L if x not in a.L]
    common = [x for x in a.L if x in b.L] if contract else []
    out = {}
    for vals in itertools.product((0, 1), repeat=len(L)):
        asg = dict(zip(L, vals))
        r = op(a.d[tuple(asg[x] for x in a.L)], b.d[tuple(asg[x] for x in b.L)])
        key = tuple(asg[x] for x in L if x not in common)
        out[key] = out[key] + r if key in out else r
    return T([x for x in L if x not in common], out)


def _unary(a, f):
    return T(a.L, {k: f(v) for k, v in a.d.items()})


class Ev:
    def __init__(self, names):
        self.n = names

    def ev(self, node):
        m = getattr(self, 'e_' + type(node).__name__)
        return m(node)

    def e_Expression(self, node):
        return self.ev(node.body)

    def e_Constant(self, node):
        return T.s(node.value)

    def e_Name(self, node):
        return T.s(self.n[node.id])

    def e_UnaryOp(self, node):
        a = self.ev(node.operand)
        return _unary(a, (lambda v: -v) if isinstance(node.op, ast.USub) else (lambda v: v))

    def e_BinOp(self, node):
        a, b = self.ev(node.left), self.ev(node.right)
        op = node.op
        if isinstance(op, ast.Add):
            return _combine(a, b, lambda x, y: x + y, False)
        if isinstance(op, ast.Sub):
            return _combine(a, b, lambda x, y: x - y, False)
        if isinstance(op, ast.Mult):
            return _combine(a, b, lambda x, y: x * y, True)
        if isinstance(op, ast.Div):
            return _combine(a, b, lambda x, y: x / y, True)
        if isinstance(op, ast.Pow):
            return _combine(a, b, lambda x, y: x ** y, False)
        if isinstance(op, ast.MatMult):
            return _combine(a, b, lambda x, y: x @ y, True)
        raise TE.TexError('op %s' % op)

    def e_Call(self, node):
        fname = node.func.id
        if fname == 'COMP':
            v = self.ev(node.args[0])
            letter = node.args[1].value
            assert v.L == ()
            vec = v.d[()]
            return T((letter,), {(0,): float(vec[0]), (1,): float(vec[1])})
        if fname == 'DEL':
            a, b = node.args[0].value, node.args[1].value
            if a == b:
                return T.s(2.0)
            return T((a, b), {(i, j): (1.0 if i == j else 0.0) for i in (0, 1) for j in (0, 1)})
        if fname in ('D', 'CU', 'CR', 'CF', 'KD'):
            return T.s(self.n[fname](*[x.value for x in node.args]))
        args = [self.ev(x) for x in node.args]
        f = self.n[fname]
        if len(args) == 1:
            return _unary(args[0], f)
        r = args[0]
        for x in args[1:]:
            r = _combine(r, x, lambda p, q: (p, q), False)
        return _unary(r, lambda pq: f(*_flatten(pq)))


def _flatten(pq):
    if isinstance(pq, tuple):
        out = []
        for x in pq:
            out.extend(_flatten(x))
        return out
    return [pq]


def value(tex, ctx, **long):
    """value of a term (no longitudinal integral left); longitudinal variables as arrays; indices contracted per product"""
    t, ph = phase_values(tex, ctx, long)
    t = longvars(t)
    if t not in _cache:
        _cache[t] = ast.parse(TE.translate(t), mode='eval')
    tree = _cache[t]
    extra = dict(ph)
    extra.update({k: v for k, v in long.items() if v is not None})
    for k in ('PP', 'QQ', 'XI', 'PLS'):
        extra.setdefault(k, np.nan)
    names = ctx.env(extra).names({})
    pt = ctx.cfg.pt
    names['D'] = lambda u, v: V2(*(pt[u] - pt[v]))
    r = Ev(names).ev(tree)
    if r.L:
        raise TE.TexError('free indices left: %s' % (r.L,))
    return TE.scalar_or_col(r.d[()], ctx.cfg)


def plus_value(tex, ctx, var, x, x0=0.0):
    """integrand with your plus prescription [m(x)]_+ f(x) = m(x) [f(x) - f(x0)]:  A(x) + m(x) [f(x) - f(x0)]"""
    t, m = mark_plus(tex)
    if m is None:
        return value(t, ctx, **{var: x})
    z = np.zeros_like(x)
    A = value(t, ctx, **{var: x, 'PLS': z})
    f = value(t, ctx, **{var: x, 'PLS': z + 1}) - A
    xs = np.full_like(x, x0)
    A0 = value(t, ctx, **{var: xs, 'PLS': z})
    f0 = value(t, ctx, **{var: xs, 'PLS': z + 1}) - A0
    env = dict(PP=None, QQ=None, XI=None, KP=ctx.KP, PI=np.pi)
    env[var] = x
    mx = eval(TE.translate(longvars(m)), {'np': np}, env)
    return A + mx * (f - f0)


def gl_log(a, b, panels=60, order=10):
    """nodes, weights for int_a^b with u = log x (a > 0)"""
    xg, wg = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(np.log(a), np.log(b), panels + 1)
    X, W = [], []
    for u0, u1 in zip(edges[:-1], edges[1:]):
        u = 0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0)
        X.append(np.exp(u))
        W.append(0.5 * (u1 - u0) * wg * np.exp(u))
    return np.concatenate(X), np.concatenate(W)


def gl_lin(a, b, panels=60, order=10):
    xg, wg = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(a, b, panels + 1)
    X, W = [], []
    for u0, u1 in zip(edges[:-1], edges[1:]):
        X.append(0.5 * (u1 - u0) * xg + 0.5 * (u1 + u0))
        W.append(0.5 * (u1 - u0) * wg)
    return np.concatenate(X), np.concatenate(W)


# ---------------------------------------------------------------- whole equations
def split_top(rhs):
    """split a right-hand side at the + and - signs that are outside every bracket"""
    t = re.sub(r'\\(notag|nonumber)', ' ', rhs).replace('\\\\', ' ').replace('&', ' ')
    terms, cur, depth, k = [], '', 0, 0
    while k < len(t):
        ch = t[k]
        if t.startswith(r'\{', k) or t.startswith(r'\}', k):
            depth += 1 if t[k + 1] == '{' else -1
            cur += t[k:k + 2]
            k += 2
            continue
        if ch in '{[(':
            depth += 1
        elif ch in '}])':
            depth -= 1
        if ch in '+-' and depth == 0 and cur.strip() and not re.search(r'(\\left|\\right|\\Bigg|\\bigg|\\Big|\\big)\s*$', cur):
            terms.append(cur)
            cur = ch
        else:
            cur += ch
        k += 1
    if cur.strip():
        terms.append(cur)
    return terms


def limit(tex, ctx):
    return float(eval(TE.translate(longvars(tex)), {'np': np}, dict(KP=ctx.KP, V=ctx.V, L=ctx.L, PI=np.pi)))


def _gl(a, b, n):
    xg, wg = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * xg + 0.5 * (b + a), 0.5 * (b - a) * wg


def nodes(a, b, order=12):
    """composite Gauss-Legendre on [a, b] (a < b): geometric panels toward both ends, panels of width <= 0.25 between"""
    X, W = [], []
    span = b - a
    near = min(0.5 * span, 0.5)
    g = np.geomspace(1e-13, 1.0, 50)
    edges = [a] + list(a + near * g) + list(b - near * g) + [b]
    mid0, mid1 = a + near, b - near
    k = mid0
    while k < mid1 - 1e-15:
        k = min(k + 0.25, mid1)
        edges.append(k)
    edges = sorted(set(edges))
    for u0, u1 in zip(edges[:-1], edges[1:]):
        if u1 - u0 <= 0:
            continue
        x, w = _gl(u0, u1, order)
        X.append(x)
        W.append(w)
    return np.concatenate(X), np.concatenate(W)


def integrate(f, a, b, pv=()):
    """int_a^b f (f takes an array); principal value at the points pv inside (a, b); orientation as written"""
    sign = 1.0
    if a > b:
        a, b, sign = b, a, -1.0
    pts = sorted(p for p in pv if a < p < b)
    tot = 0.0
    lo = a
    for s in pts:
        h = 0.5 * min(s - lo, b - s)
        x, w = nodes(lo, s - h)
        tot = tot + np.sum(w * f(x))
        # symmetric pairing around the pole; plain panels (nodes too close to s would amplify the rounding of s)
        t, wt = gl_lin(0.0, h, panels=24, order=12)
        tot = tot + np.sum(wt * (f(s + t) + f(s - t)))
        lo = s + h
    x, w = nodes(lo, b)
    tot = tot + np.sum(w * f(x))
    return sign * tot


VARS = [(r'p\^\{\+\}', 'PP'), (r'q\^\{\+\}', 'QQ'), (r'\\xi', 'XI')]


def term_value(term, ctx, pv=None):
    """one top-level term: its longitudinal integral (if any) is done numerically"""
    for vre, var in VARS:
        t, lims = strip_long_int(term, vre)
        if lims:
            assert len(lims) == 1, lims
            lo, hi, meas = lims[0]
            a, b = limit(lo, ctx), limit(hi, ctx)
            jac = 1.0 / (2 * np.pi) if r'2\pi' in meas else 1.0
            xi2 = var == 'XI' and r'\xi^{2}' in meas
            x0 = 1.0 if var == 'XI' else 0.0
            pts = pv(var) if pv else ()
            return jac * integrate(lambda x: plus_value(t, ctx, var, x, x0) / (x**2 if xi2 else 1.0), a, b, pts)
    return value(term, ctx)


def side_value(rhs, ctx, pv=None):
    return sum(term_value(t, ctx, pv) for t in split_top(rhs))


def strip_delta(t):
    r"""remove \delta^{(2)}\left[ ... \right] (and the 'x' that may follow); return text and the argument"""
    m = re.search(r'\\delta\^\{\(2\)\}\s*\\left\[', t)
    if not m:
        return t, None
    i = m.start()
    k = m.end()
    d = 1
    while d:
        if t.startswith(r'\left[', k):
            d += 1
            k += 6
            continue
        if t.startswith(r'\right]', k):
            d -= 1
            k += 7
            continue
        k += 1
    arg = t[m.end():k - 7]
    rest = re.sub(r'^\s*\\times', ' ', t[k:])
    return t[:i] + ' ' + rest, arg


def moved(cfg, point, where):
    """copy of the configuration with one point moved"""
    import copy
    c2 = copy.copy(cfg)
    c2.pt = dict(cfg.pt)
    c2.pt[point] = where
    return c2


def n_long(term):
    return len(re.findall(r'\\frac\{d(p|q)\^\{\+\}\}\{2\\pi\}|\\frac\{d\\xi\}', term))


def expand(term):
    """a term with more than one longitudinal integral: distribute the prefix and suffix over the outermost brace group"""
    if n_long(term) <= 1:
        return [term]
    for o, c in ((r'\Bigg\{', r'\Bigg\}'), (r'\bigg\{', r'\bigg\}'), (r'\Bigg[', r'\Bigg]')):
        i = term.find(o)
        if i < 0:
            continue
        d, k = 0, i
        while k < len(term):
            if term.startswith(o, k):
                d += 1
                k += len(o)
                continue
            if term.startswith(c, k):
                d -= 1
                if d == 0:
                    break
                k += len(c)
                continue
            k += 1
        pre, inner, suf = term[:i], term[i + len(o):k], term[k + len(c):]
        if n_long(pre) or n_long(suf):
            continue
        out = []
        for t in split_top(inner):
            out.extend(expand(pre + ' ( ' + t + ' ) ' + suf))
        return out
    raise TE.TexError('cannot separate the longitudinal integrals: ' + term[:200])


def side_value(rhs, ctx, pv=None):
    tot = 0.0
    for t in split_top(rhs):
        for u in expand(t):
            tot = tot + term_value(u, ctx, pv)
    return tot


def side_integrand(rhs, ctx, var, x, delta_move=None):
    """sum of the integrands (all terms must share the integral); delta_move(cfg, P) moves a point for a delta^2"""
    tot = 0.0
    for t in split_top(rhs):
        for u in expand(t):
            args = []
            while True:
                u, arg = strip_delta(u)
                if arg is None:
                    break
                args.append(arg)
            arg = args or None
            for vre, v in VARS:
                s, lims = strip_long_int(u, vre)
                if lims:
                    break
            xi2 = bool(lims) and var == 'XI' and r'\xi^{2}' in lims[0][2]
            if arg is not None and delta_move is not None:
                vals = []
                for P in np.atleast_1d(x):
                    moves, jac = delta_move(ctx.cfg, P, ctx.KP)
                    cm = ctx.cfg
                    for pt, where in moves:
                        cm = moved(cm, pt, where)
                    c2 = Ctx(cm, ctx.KP, ctx.V, ctx.L, ctx.EPS)
                    vals.append(jac * plus_value(s, c2, var, np.array([P]), 1.0 if var == 'XI' else 0.0)[0] / (P**2 if xi2 else 1.0))
                tot = tot + np.array(vals)
            else:
                xx = np.atleast_1d(x)
                tot = tot + plus_value(s, ctx, var, xx, 1.0 if var == 'XI' else 0.0) / (xx**2 if xi2 else 1.0)
    return tot


def balance(t):
    r"""remove a \Bigg\} (or \bigg\}) that closes nothing; return text and the number removed"""
    out, depth, k, removed = '', 0, 0, 0
    while k < len(t):
        if t.startswith('\\\\', k):            # a line break \\ is not an escaped brace
            out += t[k:k + 2]
            k += 2
            continue
        m = re.match(r'\\(Bigg|bigg|Big|big)?\\\{', t[k:]) if t[k] == '\\' else None
        if t.startswith(r'\{', k) or (m and m.group(0).endswith(r'\{')):
            L = 2 if t.startswith(r'\{', k) else len(m.group(0))
            depth += 1
            out += t[k:k + L]
            k += L
            continue
        m = re.match(r'\\(Bigg|bigg|Big|big)?\\\}', t[k:]) if t[k] == '\\' else None
        if m is None and t.startswith(r'\}', k):
            m = re.match(r'\\\}', t[k:])
        if m:
            if depth == 0:
                removed += 1
                k += len(m.group(0))
                continue
            depth -= 1
            out += m.group(0)
            k += len(m.group(0))
            continue
        out += t[k]
        k += 1
    if depth > 0:                      # an opening \Bigg\{ that is never closed: close it at the end
        out += r' \Bigg\}' * depth
        removed += depth
    return out, removed
