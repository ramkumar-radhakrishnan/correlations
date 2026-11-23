# s16_lib.py -- parser/evaluator for the colour factors as PRINTED in TeX (index tensors and adjoint traces),
# used by s16b (colour checks of the author's Sec. 7 and of the new Sec. 8).  Adjoint Wilson lines are real
# orthogonal, U^dagger = U^T; (T^e)_{bc} = -i f^{ebc}, (D^h)_{bc} = d^{hbc}.
#   tensor expressions : U^{ab}(x), f^{abc}, d^{abc}, \delta^{ab}, N_{c}, numbers, i, \frac{p}{q}, brackets, sums.
#   trace expressions  : {\rm Tr}\{ ... \} with T^{e}, D^{h} (each generator index occurs twice, summed),
#                        U(x), U^{\dagger}(x), (U(x)-U(y)), (U^{\dagger}(x)-U^{\dagger}(y)), products of traces.
import re, numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fe_colour import SUN
PT = {"x'": 'xp', "w'": 'wp', "y'": 'yp', 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w', 'u': 'u', 'v': 'v'}
def clean(s):
    s = s.replace('^{\\prime}', "'").replace('\\prime', "'")
    s = re.sub(r'\\(left|right)(?=[\[\]\(\)\.\\])', '', s)
    s = re.sub(r'\\(Big|big|bigg|Bigg)(?=[\[\]\(\)\\])', '', s)
    s = re.sub(r'\\(Big|big|bigg|Bigg)l?r?(?=[\[\]\(\)])', '', s)
    s = re.sub(r'\\(notag|times|quad|qquad)(?![a-zA-Z])', ' ', s)
    s = re.sub(r'\\[,;!]', ' ', s).replace('\\\\', ' ').replace('&', ' ')
    s = re.sub(r'\\color\{[a-zA-Z]+\}', '', s)
    return s
def split_idx(s):
    out = []; i = 0
    while i < len(s):
        c = s[i]; i += 1
        if c in ' {}': continue
        while i < len(s) and s[i] == "'": c += "'"; i += 1
        out.append(c)
    return out
TOK = re.compile(r"""
 (?P<rmTr>\{\\rm\s*Tr\}\\\{|\\mathrm\{Tr\}\\\{|\\operatorname\{Tr\}\\\{) |
 (?P<Ud>U\^\{\\dagger\}\((?P<udp>[^()]*)\)) |
 (?P<Uidx>U\^\{(?P<ui>[^{}\\]*)\}\((?P<up>[^()]*)\)) |
 (?P<Um>U\((?P<ump>[^()]*)\)) |
 (?P<T>T\^\{?(?P<ti>[a-z]'?)\}?) |
 (?P<D>D\^\{?(?P<di>[a-z]'?)\}?) |
 (?P<f>f\^\{(?P<fi>[^{}]*)\}) |
 (?P<d>d\^\{(?P<dd>[^{}]*)\}) |
 (?P<delta>\\delta\^\{(?P<de>[^{}()]*)\}) |
 (?P<Nc>N_\{?c\}?) |
 (?P<frac>\\frac\{(?P<fn>[0-9]+)\}\{(?P<fd>[0-9]+)\}) |
 (?P<num>[0-9]+(\.[0-9]+)?) |
 (?P<imag>(?<![a-zA-Z\\])i(?![a-zA-Z])) |
 (?P<lbr>\\\{) | (?P<rbr>\\\}) |
 (?P<lb>\[) | (?P<rb>\]) | (?P<lp>\() | (?P<rp>\)) |
 (?P<plus>\+) | (?P<minus>-) |
 (?P<ws>\s+|\{|\})
""", re.X)
def tokenize(s):
    s = clean(s); pos = 0; out = []
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: raise ValueError('cannot tokenize at: %r' % s[pos:pos + 60])
        pos = m.end(); k = m.lastgroup
        if k == 'ws': continue
        g = m.groupdict()
        if k == 'Uidx': out.append(('U', tuple(split_idx(g['ui'])), PT[g['up'].strip()]))
        elif k == 'Ud': out.append(('Ud', PT[g['udp'].strip()]))
        elif k == 'Um': out.append(('Um', PT[g['ump'].strip()]))
        elif k == 'T': out.append(('T', g['ti']))
        elif k == 'D': out.append(('D', g['di']))
        elif k == 'f': out.append(('f', tuple(split_idx(g['fi']))))
        elif k == 'd': out.append(('d', tuple(split_idx(g['dd']))))
        elif k == 'delta': out.append(('delta', tuple(split_idx(g['de']))))
        elif k == 'Nc': out.append(('Nc',))
        elif k == 'frac': out.append(('num', int(g['fn']) / int(g['fd'])))
        elif k == 'num': out.append(('num', float(m.group(0))))
        elif k == 'imag': out.append(('num', 1j))
        elif k == 'rmTr': out.append(('Tr',))
        else: out.append((k,))
    return out
class Ctx:
    def __init__(s, S, U): s.S, s.U = S, U
# ---------- tensor values: (array, [index names]) ----------
LET = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
def tprod(fs):
    fs = [f for f in fs]
    scal = 1
    arrs = []
    for v, idx in fs:
        if not idx: scal = scal * v
        else: arrs.append((v, idx))
    if not arrs: return (scal, [])
    cnt = {}
    for v, idx in arrs:
        for q in idx: cnt[q] = cnt.get(q, 0) + 1
    if any(c > 2 for c in cnt.values()): raise ValueError('index more than twice: %s' % cnt)
    names = sorted(cnt); mp = {q: LET[n] for n, q in enumerate(names)}
    free = [q for q in names if cnt[q] == 1]
    sub = ','.join(''.join(mp[q] for q in idx) for v, idx in arrs) + '->' + ''.join(mp[q] for q in free)
    r = np.einsum(sub, *[v for v, idx in arrs], optimize=True)
    return (scal * r, free)
def tsum(a, b, sgn):
    va, ia = a; vb, ib = b
    if sorted(ia) != sorted(ib): raise ValueError('sum of different free indices %s %s' % (ia, ib))
    if ia != ib: vb = np.transpose(vb, [ib.index(q) for q in ia])
    return (va + sgn * vb, ia)
class P:
    def __init__(s, toks, ctx): s.t, s.i, s.c = toks, 0, ctx
    def peek(s): return s.t[s.i] if s.i < len(s.t) else ('EOF',)
    def eat(s, k):
        t = s.peek()
        if t[0] != k: raise ValueError('expected %s got %s at %d: %s' % (k, t, s.i, s.t[max(0, s.i - 5):s.i + 5]))
        s.i += 1; return t
    def expr(s, stop=('rb', 'rp', 'EOF', 'rbr')):
        sgn = 1
        if s.peek()[0] in ('plus', 'minus'):
            sgn = 1 if s.eat(s.peek()[0])[0] == 'plus' else -1
        acc = s.term(); acc = (sgn * acc[0], acc[1])
        while s.peek()[0] in ('plus', 'minus'):
            sg = 1 if s.eat(s.peek()[0])[0] == 'plus' else -1
            acc = tsum(acc, s.term(), sg)
        return acc
    def term(s):
        fs = []
        while s.peek()[0] not in ('plus', 'minus', 'rb', 'rp', 'EOF', 'rbr'):
            fs.append(s.factor())
        if not fs: raise ValueError('empty term at %d' % s.i)
        return tprod(fs)
    def factor(s):
        t = s.peek(); S, U = s.c.S, s.c.U
        if t[0] == 'lb': s.eat('lb'); v = s.expr(); s.eat('rb'); return v
        if t[0] == 'lp': s.eat('lp'); v = s.expr(); s.eat('rp'); return v
        if t[0] == 'U': s.i += 1; return (U[t[2]], list(t[1]))
        if t[0] == 'f': s.i += 1; return (S.f, list(t[1]))
        if t[0] == 'd': s.i += 1; return (S.d, list(t[1]))
        if t[0] == 'delta': s.i += 1; return (np.eye(S.n), list(t[1]))
        if t[0] == 'Nc': s.i += 1; return (S.N, [])
        if t[0] == 'num': s.i += 1; return (t[1], [])
        if t[0] == 'Tr': return (s.trace(), [])
        raise ValueError('unexpected token %s at %d' % (t, s.i))
    # ----- traces -----
    def trace(s):
        s.eat('Tr'); chain = []
        while s.peek()[0] != 'rbr':
            chain.append(s.melem())
        s.eat('rbr')
        return evaltrace(chain, s.c)
    def melem(s):
        t = s.peek(); S, U = s.c.S, s.c.U
        if t[0] in ('T', 'D'): s.i += 1; return (t[0], t[1])
        if t[0] == 'Um': s.i += 1; return ('M', U[t[1]])
        if t[0] == 'Ud': s.i += 1; return ('M', U[t[1]].T)
        if t[0] == 'lp':
            s.eat('lp'); v = s.mexpr(); s.eat('rp'); return ('M', v)
        raise ValueError('unexpected in trace %s at %d' % (t, s.i))
    def mexpr(s):
        sgn = 1
        if s.peek()[0] in ('plus', 'minus'): sgn = 1 if s.eat(s.peek()[0])[0] == 'plus' else -1
        acc = sgn * s.mprod()
        while s.peek()[0] in ('plus', 'minus'):
            sg = 1 if s.eat(s.peek()[0])[0] == 'plus' else -1
            acc = acc + sg * s.mprod()
        return acc
    def mprod(s):
        m = None; sc = 1
        while s.peek()[0] not in ('plus', 'minus', 'rp'):
            e = s.melem() if s.peek()[0] != 'num' else ('S', s.eat('num')[1])
            if e[0] == 'S': sc = sc * e[1]; continue
            if e[0] != 'M': raise ValueError('generator inside a parenthesis')
            m = e[1] if m is None else m @ e[1]
        return sc * m
def evaltrace(chain, ctx):
    S = ctx.S
    gens = {}
    for n, e in enumerate(chain):
        if e[0] in ('T', 'D'): gens.setdefault(e[1], []).append(n)
    for k, v in gens.items():
        if len(v) != 2: raise ValueError('generator index %s occurs %d times in a trace' % (k, len(v)))
    # einsum: matrix indices m0..; generator letters
    L = iter('abcdefghijklmnopqrstuvwxyz'); gl = {k: next(L) for k in gens}
    nmat = len(chain); mi = [next(L) for _ in range(nmat)]
    subs = []; ops = []
    for n, e in enumerate(chain):
        a, b = mi[n], mi[(n + 1) % nmat]
        if e[0] == 'T': subs.append(gl[e[1]] + a + b); ops.append(S.T)
        elif e[0] == 'D': subs.append(gl[e[1]] + a + b); ops.append(S.D)
        else: subs.append(a + b); ops.append(e[1])
    return np.einsum(','.join(subs) + '->', *ops, optimize=True)
def ev(tex, S, U):
    """evaluate a printed colour expression; returns (value, free index list)"""
    p = P(tokenize(tex), Ctx(S, U)); v = p.expr()
    if p.peek()[0] != 'EOF': raise ValueError('trailing tokens: %s' % (p.t[p.i:p.i + 8],))
    return v
def evs(tex, S, U):
    v, idx = ev(tex, S, U)
    assert not idx, idx
    return v
def bracket_after(s, start):
    """s[start] is at an opening bracket token (\\Big[, \\Bigg[, \\bigg[, [ , \\Big\\{ ...); return (content, end)"""
    m = re.match(r'\\(?:Big|big|bigg|Bigg)?(\[|\\\{|\()', s[start:])
    op = m.group(1); j = start + m.end()
    cl = {'[': ']', '\\{': '\\}', '(': ')'}[op]
    depth = 1; k = j
    while k < len(s):
        if s.startswith(op, k) and not (op == '(' and False): depth += 1; k += len(op); continue
        if s.startswith(cl, k):
            depth -= 1
            if depth == 0:
                # strip a trailing \Big etc. that belongs to the closing bracket
                content = re.sub(r'\\(?:Big|big|bigg|Bigg)$', '', s[j:k])
                return content, k + len(cl)
            k += len(cl); continue
        k += 1
    raise ValueError('unbalanced')
# ---------- random Wilson lines ----------
PTS = ['x', 'xp', 'y', 'yp', 'z', 'w', 'wp']
def randUs(S, rng): return {p: S.adj(S.randV(rng)) for p in PTS}
