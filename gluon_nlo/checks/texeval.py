"""Mechanical evaluation of the author's LaTeX terms on the lattice.

Each align-block is parsed into a sum of monomials (brackets are expanded
distributively), and every monomial is contracted with numpy.einsum.  Nothing is
retyped by hand, so transcription errors cannot creep in.
"""
import re
import string
import numpy as np
import opt_einsum as oe

STRIP = [r'\notag', r'\\', '&', r'\,', r'\;', r'\!', r'\quad', r'\qquad',
         r'\left', r'\right', r'\Bigg', r'\bigg', r'\Big', r'\big']

def clean(s):
    s = s.replace(r'\mathcal K', r'\mathcal{K}').replace(r'\mathcal{K} ', r'\mathcal{K}')
    s = s.replace('N^{2}_{c}', 'N_{c}^{2}')
    for t in STRIP:
        s = s.replace(t, '')
    s = re.sub(r'\s+', '', s)
    return s

TOK = [
    ('KDOT', r'\\mathcal\{K\}\(([^()]*)\)\\cdot\\mathcal\{K\}\(([^()]*)\)'),
    ('KK',   r'K\(([^,()]+),([^,()]+),([^,()]+)\)'),
    ('INV2', r'\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}'),
    ('F',    r'f\^\{([^}]*)\}'),
    ('U',    r'U\^\{([^}]*)\}\(([^()]*)\)'),
    ('RHO',  r'\\rho\^\{([^}]*)\}\(([^()]*)\)'),
    ('NC',   r'N_\{c\}'),
    ('NUM',  r'\d+'),
    ('OP',   r'[+\-]'),
    ('LB',   r'\['),
    ('RB',   r'\]'),
]

def split_idx(s):
    s = s.replace('{', '').replace('}', '').replace(' ', '')
    out = []
    for ch in s:
        if ch == "'":
            out[-1] += "'"
        else:
            out.append(ch)
    return out

def split_pos(s):
    a, b = s.split('-')
    return a, b

def tokenize(s):
    toks, i = [], 0
    while i < len(s):
        for name, pat in TOK:
            m = re.match(pat, s[i:])
            if m:
                toks.append((name, m.groups() if m.groups() else m.group(0)))
                i += m.end()
                break
        else:
            raise ValueError('cannot tokenize at: ' + s[i:i + 60])
    return toks

def parse_sum(toks, i):
    """returns list of monomials (coef, factors) and new index"""
    mons = []
    sign = 1
    while i < len(toks):
        name, val = toks[i]
        if name == 'RB':
            return mons, i
        if name == 'OP':
            sign = 1 if val == '+' else -1
            i += 1
            continue
        term, i = parse_term(toks, i)
        mons += [(sign * c, f) for c, f in term]
        sign = 1
    return mons, i

def parse_term(toks, i):
    cur = [(1.0, [])]
    while i < len(toks):
        name, val = toks[i]
        if name in ('OP', 'RB'):
            break
        if name == 'LB':
            sub, i = parse_sum(toks, i + 1)
            assert toks[i][0] == 'RB'
            i += 1
            cur = [(c1 * c2, f1 + f2) for c1, f1 in cur for c2, f2 in sub]
            continue
        if name == 'NUM':
            cur = [(c * float(val), f) for c, f in cur]
        elif name == 'NC':
            cur = [(c, f + [('NC',)]) for c, f in cur]
        else:
            cur = [(c, f + [(name, val)]) for c, f in cur]
        i += 1
    return cur, i

def parse(tex):
    toks = tokenize(clean(tex))
    mons, i = parse_sum(toks, 0)
    assert i == len(toks), 'unbalanced brackets'
    return mons

class Evaluator:
    def __init__(self, lat, P2=None, P3=None):
        self.lat = lat
        self.P2 = P2      # P2[u,b,v,c]  = < rho^b(u) rho^c(v) >  (ordered)
        self.P3 = P3

    def monomial(self, coef, factors, integ, phase):
        lat = self.lat
        letters = iter(string.ascii_letters)
        names = {}
        def L(name):
            if name not in names:
                names[name] = next(letters)
            return names[name]
        ops, arrs = [], []
        nc = 1.0
        rhos = []
        for fac in factors:
            kind = fac[0]
            if kind == 'NC':
                nc *= lat.N
            elif kind == 'KDOT':
                (a, b), (c, d) = split_pos(fac[1][0]), split_pos(fac[1][1])
                v = L(('vec', len(names)))
                ops += ['p' + a, 'p' + b, v], ['p' + c, 'p' + d, v]
                arrs += [lat.KV, lat.KV]
            elif kind == 'KK':
                ops.append(['p' + fac[1][0], 'p' + fac[1][1], 'p' + fac[1][2]]); arrs.append(lat.KK)
            elif kind == 'INV2':
                a, b = split_pos(fac[1][0])
                ops.append(['p' + a, 'p' + b]); arrs.append(lat.INV2)
            elif kind == 'F':
                ix = split_idx(fac[1][0]); assert len(ix) == 3
                ops.append(['c' + q for q in ix]); arrs.append(lat.f)
            elif kind == 'U':
                ix = split_idx(fac[1][0]); assert len(ix) == 2
                ops.append(['p' + fac[1][1]] + ['c' + q for q in ix]); arrs.append(lat.U)
            elif kind == 'RHO':
                ix = split_idx(fac[1][0]); assert len(ix) == 1
                rhos.append(('p' + fac[1][1], 'c' + ix[0]))
            else:
                raise ValueError(kind)
        if len(rhos) == 2:
            ops.append([rhos[0][0], rhos[0][1], rhos[1][0], rhos[1][1]]); arrs.append(self.P2)
        elif len(rhos) == 3:
            ops.append([r for rr in rhos for r in rr]); arrs.append(self.P3)
        elif len(rhos) != 0:
            raise ValueError('rho count %d' % len(rhos))
        for pos, sgn in phase:
            ops.append(['p' + pos]); arrs.append(lat.ph if sgn > 0 else lat.ph.conj())
        used = {q for o in ops for q in o if q.startswith('p')}
        missing = [p for p in integ if 'p' + p not in used]
        extra = [q for q in used if q[1:] not in integ]
        assert not extra, ('position not integrated', extra)
        sub = ','.join(''.join(L(q) for q in o) for o in ops) + '->'
        val = oe.contract(sub, *arrs, optimize='auto')
        return coef * nc * val * lat.S ** len(missing)

    def total(self, mons, integ, phase):
        return sum(self.monomial(c, f, integ, phase) for c, f in mons)

def blocks(texfile):
    """Return the list of align blocks (raw text) of the file."""
    src = open(texfile).read()
    return re.findall(r'\\begin\{align\}(.*?)\\end\{align\}', src, re.S)

def body(block):
    """Strip the 'lhs = prefactor \\int ... [' head and the final ']' of a block."""
    m = re.search(r'e\^\{-ik\\cdot\s*\(w-w\'\)\}|e\^\{-ik\(w-w\'\)\}', block)
    assert m, 'no phase found'
    rest = block[m.end():]
    return rest

# ---------------------------------------------------------------- traces
TR_PAT = r'\{\\rmTr\}\\\{(.*?)\\\}'
TOK.insert(0, ('TR', TR_PAT))

def tr_chain(content):
    """'T^eU^\\dagger(w)U(w')' -> list of ('T',e) / ('Ud',pos) / ('U',pos)"""
    out, i = [], 0
    pats = [('T', r'T\^\{?([a-z]\'?)\}?'), ('Ud', r'U\^\\dagger\(([^()]*)\)'), ('U', r'U\(([^()]*)\)')]
    while i < len(content):
        for nm, p in pats:
            m = re.match(p, content[i:])
            if m:
                out.append((nm, m.group(1))); i += m.end(); break
        else:
            raise ValueError('bad trace content ' + content[i:])
    return out

_old_monomial = Evaluator.monomial
def _monomial(self, coef, factors, integ, phase):
    lat = self.lat
    newf, extra_ops = [], []
    for fac in factors:
        if fac[0] == 'TR':
            extra_ops.append(tr_chain(fac[1][0]))
        else:
            newf.append(fac)
    if not extra_ops:
        return _old_monomial(self, coef, factors, integ, phase)
    # build einsum by hand: reuse the old code path by turning traces into U/T tensors
    self._tr = extra_ops
    return _old_monomial_tr(self, coef, newf, integ, phase, extra_ops)

def _old_monomial_tr(self, coef, factors, integ, phase, chains):
    lat = self.lat
    letters = iter(string.ascii_letters)
    names = {}
    def L(name):
        if name not in names:
            names[name] = next(letters)
        return names[name]
    ops, arrs = [], []
    # non-trace factors: reuse logic of the original evaluator by a small inline copy
    nc = 1.0; rhos = []
    for fac in factors:
        kind = fac[0]
        if kind == 'NC': nc *= lat.N
        elif kind == 'KDOT':
            (a, b), (c, d) = split_pos(fac[1][0]), split_pos(fac[1][1])
            v = ('vec', len(names)); ops += [['p'+a,'p'+b,v],['p'+c,'p'+d,v]]; arrs += [lat.KV, lat.KV]
        elif kind == 'KK':
            ops.append(['p'+fac[1][0],'p'+fac[1][1],'p'+fac[1][2]]); arrs.append(lat.KK)
        elif kind == 'INV2':
            a, b = split_pos(fac[1][0]); ops.append(['p'+a,'p'+b]); arrs.append(lat.INV2)
        elif kind == 'RHO':
            ix = split_idx(fac[1][0]); rhos.append(('p'+fac[1][1], 'c'+ix[0]))
        else: raise ValueError(kind)
    if len(rhos) == 2:
        ops.append([rhos[0][0], rhos[0][1], rhos[1][0], rhos[1][1]]); arrs.append(self.P2)
    for ci, ch in enumerate(chains):
        n = len(ch)
        for j, (nm, arg) in enumerate(ch):
            i0, i1 = ('m', ci, j), ('m', ci, (j + 1) % n)
            if nm == 'T':
                ops.append(['c' + arg, i0, i1]); arrs.append(lat.T)
            elif nm == 'U':
                ops.append(['p' + arg, i0, i1]); arrs.append(lat.U)
            else:
                ops.append(['p' + arg, i1, i0]); arrs.append(lat.U)   # U^dagger = U^T (real)
    for pos, sgn in phase:
        ops.append(['p' + pos]); arrs.append(lat.ph if sgn > 0 else lat.ph.conj())
    used = {q for o in ops for q in o if isinstance(q, str) and q.startswith('p')}
    missing = [p for p in integ if 'p' + p not in used]
    sub = ','.join(''.join(L(q) for q in o) for o in ops) + '->'
    val = oe.contract(sub, *arrs, optimize='auto')
    return coef * nc * val * lat.S ** len(missing)

Evaluator.monomial = _monomial

def singlet_body(block):
    b = body(block)
    b = b.replace(r'\langle', '').replace(r'\rangle', '')
    return b
