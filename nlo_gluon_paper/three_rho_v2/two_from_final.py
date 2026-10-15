"""The two-rho terms of the final three-rho equation (ThreeRho_final.tex), as one equation.
Reads the black and blue lines of ThreeRho_final.tex, reorders each in the order written there, collects equal
terms row by row, writes TwoRho_from_final.tex, then re-parses the equation it wrote and checks it."""
import re
import sys
import sympy as sp
import env
from core import term, symmetrize, Config, pointwise_total
from book import Book
from engine import tex_colour, Nc
import final_eq as FE
sys.path.insert(0, env.ALL2)
import make_doc as MD
import sets as S2
import target as TG

ROWS = [('A,B', ('4rho.a', '4rho.b'), r"From four $\rho$'s ($d^3N_a$, $d^3N_b$)"),
        ('I.II', ('I.II.1', 'I.II.2'), 'Group I, Row II'),
        ('I.III', ('I.III',), 'Group I, Row III'),
        ('D', ('D',), r'Group I, Row III: $\mathcal D$'),
        ('I.IV', ('I.IV',), 'Group I, Row IV'),
        ('I.V', ('I.V.1', 'I.V.2'), 'Group I, Row V'),
        ('II.I', ('II.I',), 'Group II, Row I'),
        ('III.I', ('III.I',), 'Group III, Row I'),
        ('III.II', ('III.II',), 'Group III, Row II'),
        ('III.III', ('III.III',), 'Group III, Row III'),
        ('III.V', ('III.V',), 'Group III, Row V'),
        ('III.VI', ('III.VI',), 'Group III, Row VI')]


def row_entries(final):
    out = []
    for row, blocks, title in ROWS:
        B = Book()
        for r in final:
            if r['block'] in blocks:
                for t in FE.terms_of(r):
                    for s in symmetrize(t, (2,)):
                        B.add(s, row)
        B.finalize()
        bk = B.dump()
        ents = []
        for key in MD.ordered(bk):
            for net in bk[key]:
                c = sp.expand(net['groups'].get(row, 0))
                if c != 0:
                    ents.append((key, c, net['col']))
        out.append((row, title, ents))
    return out


# ------------------------------------------------------------------ rendering
def frac(q):
    q = sp.Rational(q)
    if q.q == 1:
        return '%d' % abs(q.p)
    return r'\tfrac{%d}{%d}' % (abs(q.p), q.q)


def coef_tex(c):
    """c = sum_k a_k Nc^k with rational a_k; explicit signs, parsable by coef_parse"""
    poly = sp.Poly(sp.expand(c), Nc)
    terms = []
    for (k,), a in sorted(poly.terms(), key=lambda t: t[0][0]):
        a = sp.Rational(a)
        s = '-' if a < 0 else '+'
        nc = '' if k == 0 else ('N_c' if k == 1 else 'N_c^{%d}' % k)
        if k > 0 and abs(a) == 1:
            body = nc
        elif k > 0:
            body = frac(a) + r'\,' + nc
        else:
            body = frac(a)
        terms.append(s + body)
    if len(terms) == 1:
        return terms[0][0] if terms[0][1:] == '1' else terms[0]
    inner = ''.join(terms)
    inner = inner[1:] if inner.startswith('+') else inner
    return r'+\big(%s\big)' % inner


def kern_tex(key):
    em, soft = MD.split_key(key)
    out = [e.replace(r'\KK', r'\mathcal{K}') for e in em]
    for s in soft:
        m = re.match(r'K\(([^,]+),([^,]+),z\)', s)
        a, b = m.group(1), m.group(2)
        out.append(r'\tfrac{1}{(%s-z)^{2}}' % a if a == b else r'\mathcal{K}(%s-z)\cdot\mathcal{K}(%s-z)' % (a, b))
    return r'\,'.join(out)


def write(path, R):
    eq = [r'\begin{align}',
          r"\frac{d^{3}N^{2\rho\,(3\rho)}_{\rm NLO}}{d^{3}k}={}&\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log\Big[\frac{\vee}{\Lambda}\Big]\int_{x,x',z,w,w'}e^{-ik(w'-w)}\,\Bigg[\notag\\"]
    body = []
    n = 0
    for row, title, ents in R:
        body.append(r'&\quad\underline{\text{%s}}\notag\\' % title)
        for key, c, col in ents:
            n += 1
            body.append(r'&\quad %s\;%s\;%s\notag\\' % (coef_tex(c), kern_tex(key), tex_colour(col)))
    body[-1] = body[-1].replace(r'\notag\\', r'\Bigg]\label{eq:two}')
    eq += body + [r'\end{align}']
    doc = open('two_final_static.tex').read().replace('@@EQ@@', '\n'.join(eq)).replace('@@N@@', str(n))
    open(path, 'w').write(doc)
    return n


# ------------------------------------------------------------------ re-parse and check
CT = re.compile(r"([+-])(?:\\tfrac\{(\d+)\}\{(\d+)\}|(\d+))?(?:\\,)?(N_c(?:\^\{(\d+)\})?)?")
KRE = re.compile(r"\\mathcal\{K\}\(([^)]*)\)")
INV = re.compile(r"\\tfrac\{1\}\{\((\w'?)-z\)\^\{2\}\}")
FRE = re.compile(r"f\^\{(\w)(\w)(\w)\}")
URE = re.compile(r"U\^\{(\w)(\w)\}\((\w'?)\)")
RRE = re.compile(r"\\rho\^\{(\w)\}\((\w'?)\)")


def coef_parse(s, ncval=3):
    s = s.strip()
    if s.startswith(r'+\big('):
        s = s[len(r'+\big('):-len(r'\big)')]
        if not s.startswith('-'):
            s = '+' + s
    tot = sp.Integer(0)
    pos = 0
    while pos < len(s):
        m = CT.match(s, pos)
        assert m and m.end() > pos, s[pos:]
        sign = -1 if m.group(1) == '-' else 1
        a = sp.Rational(int(m.group(2)), int(m.group(3))) if m.group(2) else (sp.Integer(int(m.group(4))) if m.group(4) else sp.Integer(1))
        k = 0 if not m.group(5) else (int(m.group(6)) if m.group(6) else 1)
        tot += sign * a * ncval ** k
        pos = m.end()
    return tot


def reparse(path):
    tex = open(path).read()
    eq = tex[tex.index(r'\begin{align}'):tex.index(r'\end{align}')]
    row = None
    titles = {t: r for r, _, t in ROWS}
    out = []
    for ln in eq.split('\n'):
        m = re.match(r'^&\\quad\\underline\{\\text\{(.*)\}\}\\notag\\\\$', ln)
        if m:
            row = titles[m.group(1)]
            continue
        m = re.match(r'^&\\quad (.*?)\\;(.*?)\\;(.*?)(?:\\notag\\\\|\\Bigg\]\\label\{eq:two\})$', ln)
        if not m:
            continue
        c = coef_parse(m.group(1))
        ktx = m.group(2)
        kern = []
        ks = [tuple(s.split('-')) for s in KRE.findall(ktx)]
        for i in range(0, len(ks), 2):
            kern.append(('KK', ks[i][0], ks[i][1], ks[i + 1][0], ks[i + 1][1]))
        for a in INV.findall(ktx):
            kern.append(('KK', a, 'z', a, 'z'))
        col = m.group(3)
        out.append((row, term(float(c), URE.findall(col), FRE.findall(col), RRE.findall(col), kern, row)))
    return out


def check(path, final):
    P2 = reparse(path)
    A = S2.A_sets()
    lab = {'A,B': 'A,B', 'D': 'D'}
    other = [t for k in ('g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III', '4rho') for t in A[k]]
    J = TG.jimwlk()
    direct = [s for r in final for t in FE.terms_of(r) for s in symmetrize(t, (2,))]
    ok = True
    for seed in (3, 7, 11):
        cfg = Config(seed, True)
        a = pointwise_total([t for _, t in P2], cfg)[0].real
        b = pointwise_total(direct, cfg)[0].real
        tot = a + pointwise_total(other, cfg)[0].real
        j = pointwise_total(J, cfg)[0].real
        print('seed %2d: written equation %.8f  direct reordering %.8f  | + genuine + 4rho %.8f  JIMWLK x LO %.8f' % (seed, a, b, tot, j))
        ok &= abs(a - b) < 1e-9 * (1 + abs(b)) and abs(tot - j) < 1e-9 * (1 + abs(j))
    cfg = Config(3, True)
    for row, _, _ in ROWS:
        a = pointwise_total([t for r, t in P2 if r == row], cfg)[0].real
        key = row if row in ('A,B', 'D') else '3:' + row
        b = pointwise_total(A[key], cfg)[0].real
        same = abs(a - b) < 1e-9 * (1 + abs(b))
        ok &= same
        print('   %-8s written %.8f  earlier (TwoRho_all) %.8f  %s' % (row, a, b, 'same' if same else 'DIFFERENT'))
    print('ALL CHECKS PASS' if ok else 'CHECK FAILED')
    return ok


if __name__ == '__main__':
    final = [r for r in FE.reparse('ThreeRho_final.tex') if r['color'] != 'red']
    R = row_entries(final)
    n = write('TwoRho_from_final.tex', R)
    print('entries', n, {row: len(e) for row, _, e in R})
    check('TwoRho_from_final.tex', final)
