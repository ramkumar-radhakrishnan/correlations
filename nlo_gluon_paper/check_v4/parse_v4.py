"""Parser of every two-rho term of the new Evolution_terms.tex (Secs. 1-4, including the factored
'kernel x [ ... ]' forms of Secs. 3, 4).  Units N = (1/(2pi)^3) g^4/(16 pi^5) (1/k+) log(vee/Lambda).
Returns a list of equations: dict(sec, sub, subsub, n, name, terms) with terms in the format of ../two_rho/core.py.
The colour-singlet (trace) forms are parsed separately by parse_traces()."""
import os
import re
import env  # noqa: F401
from core import term

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'input', 'Evolution_terms.tex')

TOK = re.compile(r"""
 (?P<KK>\\mathcal\{K\}\((?P<k1>\w'?)-(?P<k2>\w'?)\))
|(?P<INV>\\frac\{1\}\{\((?P<i1>\w'?)-(?P<i2>\w'?)\)\^\{2\}\})
|(?P<FR>\\frac\{(?P<fn>\d+)\}\{(?P<fd>\d+)\})
|(?P<F>f\^\{(?P<f1>\w)(?P<f2>\w)(?P<f3>\w)\})
|(?P<U>U\^\{(?P<u1>\w)(?P<u2>\w)\}\((?P<up>\w'?)\))
|(?P<R>\\rho\^\{(?P<r1>\w)\}\((?P<rp>\w'?)\))
|(?P<NC>N_\{?c\}?)
|(?P<BO>\\Big\[)
|(?P<BC>\\Big\])
|(?P<INT>(?<![\w{^])\d+(?![\w}]))
|(?P<S>[+-])
""", re.X)


class Cur:
    def __init__(self, sign=1):
        self.sign, self.coef, self.ncp = sign, 1.0, 0
        self.ks, self.inv, self.Us, self.fs, self.rs = [], [], [], [], []

    def empty(self):
        return not (self.ks or self.inv or self.Us or self.fs or self.rs)


def kern_of(ks, inv):
    assert len(ks) % 2 == 0, ks
    out = [('KK', ks[i][0], ks[i][1], ks[i + 1][0], ks[i + 1][1]) for i in range(0, len(ks), 2)]
    out += [('KK', a, b, a, b) for a, b in inv]
    return out


def parse_body(body, tag):
    """body: text between the outer \\Bigg[ and \\Bigg] -> list of terms"""
    out = []
    pre = None          # (coef, kern) of an open \Big[ block
    cur = Cur()

    def emit():
        nonlocal cur
        assert len(cur.rs) == 2, (tag, cur.rs, cur.Us)
        kern = kern_of(cur.ks, cur.inv)
        c = cur.sign * cur.coef
        if pre is not None:
            c *= pre[0]
            kern = pre[1] + kern
        assert len(kern) == 2, (tag, kern)
        out.append(term(c, [(m['u1'], m['u2'], m['up']) for m in cur.Us], [(m['f1'], m['f2'], m['f3']) for m in cur.fs],
                        [(m['r1'], m['rp']) for m in cur.rs], kern, tag, ncp=cur.ncp))
        cur = Cur()

    for m in TOK.finditer(body):
        k = m.lastgroup
        if k == 'S':
            if len(cur.rs) == 2:
                emit()
            assert cur.empty() and cur.coef == 1.0, (tag, body[max(0, m.start() - 80):m.start() + 20])
            cur.sign = -1 if m.group() == '-' else 1
        elif k == 'FR':
            cur.coef *= int(m['fn']) / int(m['fd'])
        elif k == 'INT':
            cur.coef *= int(m.group())
        elif k == 'NC':
            cur.ncp += 1
        elif k == 'KK':
            if len(cur.rs) == 2:
                emit()
            cur.ks.append((m['k1'], m['k2']))
        elif k == 'INV':
            if len(cur.rs) == 2:
                emit()
            cur.inv.append((m['i1'], m['i2']))
        elif k == 'F':
            cur.fs.append(m)
        elif k == 'U':
            cur.Us.append(m)
        elif k == 'R':
            cur.rs.append(m)
        elif k == 'BO':
            assert not (cur.Us or cur.fs or cur.rs) and pre is None, tag
            pre = (cur.sign * cur.coef, kern_of(cur.ks, cur.inv))
            cur = Cur()
        elif k == 'BC':
            if len(cur.rs) == 2:
                emit()
            assert cur.empty(), tag
            pre = None
    if len(cur.rs) == 2:
        emit()
    assert cur.empty(), (tag, cur.rs)
    return out


HEAD = re.compile(r"\\section\{([^}]*)\}|\\subsection\{((?:[^{}]|\{[^}]*\})*)\}|\\subsubsection\{([^}]*)\}|\\textbf\{([^}]*)\}"
                  r"|\\begin\{align\}(.*?)\\end\{align\}|(Color singlet form|In color singlet form)", re.S)
NAME = re.compile(r"\\frac\{d\^\{?3\}?N_?\{?(.*?)\}?\}\{d\^\{?3\}?k\}", re.S)


def equations(path=SRC):
    tex = open(path).read()
    tex = tex[tex.index(r'\begin{document}'):]
    h = dict(sec='', sub='', subsub='')
    eqs = []
    singlet_next = False
    for m in HEAD.finditer(tex):
        if m.group(1) is not None:
            h = dict(sec=m.group(1), sub='', subsub='')
        elif m.group(2) is not None:
            h['sub'], h['subsub'] = m.group(2), ''
        elif m.group(3) is not None:
            h['subsub'] = m.group(3)
        elif m.group(4) is not None:
            h['sub'], h['subsub'] = m.group(4), ''
        elif m.group(6) is not None:
            singlet_next = True
        else:
            body = m.group(5)
            nm = NAME.search(body)
            name = re.sub(r'\\[,!]|\\rm|\\mathrm|\s|\{|\}', '', nm.group(1)) if nm else ''
            n = len(eqs)
            if r'\Bigg[' not in body:          # e.g. the total in the corrected copy
                eqs.append(dict(h, n=n, name=name, kind='other', text=body, terms=[]))
                continue
            i0 = body.index(r'\Bigg[')
            i1 = body.rindex(r'\Bigg]')
            inner = body[i0 + len(r'\Bigg['):i1]
            if singlet_next:
                eqs.append(dict(h, n=n, name=name, kind='singlet', text=inner, terms=[]))
                singlet_next = False
                continue
            tag = '%s|%s|%s|%d' % (h['sec'], h['sub'], h['subsub'], n)
            eqs.append(dict(h, n=n, name=name, kind='terms', terms=parse_body(inner, tag), text=inner))
    return eqs


# ---------------------------------------------------------------- colour-singlet forms
TTOK = re.compile(r"""
 (?P<KK>\\mathcal\s?K\((?P<k1>\w'?)-(?P<k2>\w'?)\))
|(?P<KF>K\((?P<a>\w'?),(?P<b>\w'?),z\))
|(?P<TR>\{\\rm\s?Tr\}\\\{)
|(?P<BO>\\Big\[)
|(?P<BC>\\Big\])
|(?P<NC>N_\{?c\}?)
|(?P<INT>(?<![\w{^])\d+(?![\w}]))
|(?P<S>[+-])
""", re.X)
LET = re.compile(r"U\^\\dagger\((\w'?)\)|U\((\w'?)\)|T\^e|\(|\)|-")


def trace_word(s):
    """'U^\\dagger(x)U(y)T^e...' or '(U^\\dagger(a)-U^\\dagger(w))(U(b)-U(w'))' -> list of (coef, word)"""
    if '(U' in s:
        m = re.fullmatch(r"\(U\^\\dagger\((\w'?)\)-U\^\\dagger\((\w'?)\)\)\(U\((\w'?)\)-U\((\w'?)\)\)", s)
        assert m, s
        a, w, b, wp = m.groups()
        return [(s1 * s2, (('Ud', p), ('U', q))) for s1, p in ((1, a), (-1, w)) for s2, q in ((1, b), (-1, wp))]
    w = []
    for m in re.finditer(r"U\^\\dagger\((\w'?)\)|U\((\w'?)\)|T\^e", s):
        if m.group(1):
            w.append(('Ud', m.group(1)))
        elif m.group(2):
            w.append(('U', m.group(2)))
        else:
            w.append(('T', 'e'))
    assert ''.join(re.findall(r"U\^\\dagger\(\w'?\)|U\(\w'?\)|T\^e", s)) == s.replace(' ', ''), s
    return [(1, tuple(w))]


def split_traces(text):
    """find every Tr{...} with balanced braces; returns list of (start, end, content)"""
    out = []
    i = 0
    key = r'{\rm Tr}\{'
    while True:
        j = text.find(key, i)
        if j < 0:
            break
        k = j + len(key)
        depth = 1
        while depth:
            if text.startswith(r'\{', k):
                depth += 1
                k += 2
            elif text.startswith(r'\}', k):
                depth -= 1
                k += 2
            else:
                k += 1
        out.append((j, k, text[j + len(key):k - 2]))
        i = k
    return out


def parse_traces(text):
    """-> list of (coef, ncp, kern, [words...]) with words a list of (coef, word) alternatives per trace
    (a product of traces is the product of the sums).  Kernels: ('KK',a,b,c,d) for K(a-b).K(c-d), K(a,b,z) -> ('KK',a,'z',b,'z')."""
    trs = split_traces(text)
    # replace traces by placeholders
    s, last, words = '', 0, []
    for j, k, c in trs:
        s += text[last:j] + ' @%d@ ' % len(words)
        words.append(trace_word(c.replace(' ', '')))
        last = k
    s += text[last:]
    tok = re.compile(TTOK.pattern + r"|(?P<W>@(?P<wi>\d+)@)", re.X)
    out = []
    pre = None
    cur = dict(sign=1, coef=1, ncp=0, ks=[], kf=[], ws=[])

    def emit():
        nonlocal cur
        kern = [('KK', cur['ks'][i][0], cur['ks'][i][1], cur['ks'][i + 1][0], cur['ks'][i + 1][1])
                for i in range(0, len(cur['ks']), 2)] + [('KK', a, 'z', b, 'z') for a, b in cur['kf']]
        c = cur['sign'] * cur['coef']
        if pre is not None:
            c *= pre[0]
            kern = pre[1] + kern
        assert len(kern) == 2, kern
        out.append((c, cur['ncp'], kern, cur['ws']))
        cur = dict(sign=1, coef=1, ncp=0, ks=[], kf=[], ws=[])

    for m in tok.finditer(s):
        k = m.lastgroup
        if k == 'S':
            if cur['ws']:
                emit()
            cur['sign'] = -1 if m.group() == '-' else 1
        elif k == 'INT':
            cur['coef'] *= int(m.group())
        elif k == 'NC':
            cur['ncp'] += 1
        elif k == 'KK':
            if cur['ws']:
                emit()
            cur['ks'].append((m['k1'], m['k2']))
        elif k == 'KF':
            if cur['ws']:
                emit()
            cur['kf'].append((m['a'], m['b']))
        elif k == 'W':
            cur['ws'].append(words[int(m['wi'])])
        elif k == 'BO':
            kern = [('KK', cur['ks'][i][0], cur['ks'][i][1], cur['ks'][i + 1][0], cur['ks'][i + 1][1])
                    for i in range(0, len(cur['ks']), 2)] + [('KK', a, 'z', b, 'z') for a, b in cur['kf']]
            pre = (cur['sign'] * cur['coef'], kern)
            cur = dict(sign=1, coef=1, ncp=0, ks=[], kf=[], ws=[])
        elif k == 'BC':
            if cur['ws']:
                emit()
            pre = None
    return out


if __name__ == '__main__':
    from collections import Counter
    E = equations()
    for e in E:
        print('%-28s %-28s %-22s %-22s %3d %s' % (e['sec'][:28], e['sub'][:28], e['subsub'][:22], e['name'][:22],
                                                 len(e['terms']), e['kind']))
    print('total terms', sum(len(e['terms']) for e in E))
    for e in E:
        if e['kind'] == 'singlet':
            P = parse_traces(e['text'])
            print('singlet form:', e['sec'], len(P), 'monomials')
