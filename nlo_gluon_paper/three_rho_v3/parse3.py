"""Parse input/three_rho_all_v3.tex: all three-rho terms proportional to log(vee/Lambda), with c.c., as sent.
Commented lines (%) are dropped.  Every term becomes a dict in the format of ../two_rho/core.py.

Coefficient of a term = (overall prefactor, +-i) x sign x (2 if '2\\,') x (1/2 if '\\frac{1}{2}') x (2 if it contains an
anticommutator {rho, rho}).  Only the fully symmetrised (classical) part of the charges is checked here, so
1/2 {A,B} C counts as ABC and {A,B} C as 2 ABC.  '\\-' at the start of a term is read as a minus sign."""
import os
import re
import env  # noqa
from core import term

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'input', 'three_rho_all_v3.tex')

LABEL = re.compile(r"\{\\color\{(\w+)\}\s*\\(alpha|beta|gamma)\s*('?)\s*\}")
KRE = re.compile(r"\\mathcal\{K\}\(([^)]*)\)")
FRE = re.compile(r"f\^\{\s*(\w)\s*(\w)\s*(\w)\s*\}")
RHO = re.compile(r"\\rho\^\{(\w)\}\((\w'?)\)")
URE = re.compile(r"U\^\{\s*(\w)\s*(\w)\s*\}\((\w'?)\)")

CLASS_KERN = {
    'alpha': [('KK', "x'", "w'", 'x', 'w'), ('KK', 'x', 'z', 'y', 'z')],
    "alpha'": [('KK', "x'", "w'", 'x', 'w'), ('KK', "x'", 'z', 'y', 'z')],
    'beta': [('KK', 'x', 'w', "x'", "w'"), ('KK', 'w', 'z', 'y', 'z')],
    "beta'": [('KK', 'x', 'w', "x'", "w'"), ('KK', "w'", 'z', 'y', 'z')],
}
CLASS_TEX = {'alpha': r'\alpha', "alpha'": r"\alpha'", 'beta': r'\beta', "beta'": r"\beta'"}

BLOCKS = [('4rho.a', r'From four $\rho$\'s, $d^3N_a$'), ('4rho.b', r'From four $\rho$\'s, $d^3N_b$'),
          ('I.II.1', 'Group I, Row II, part 1'), ('I.II.2', 'Group I, Row II, part 2'),
          ('I.III', 'Group I, Row III'), ('I.IV', 'Group I, Row IV'),
          ('I.V.1', 'Group I, Row V, part 1'), ('I.V.2', 'Group I, Row V, part 2'),
          ('II.I', 'Group II, Row I'), ('III.I', 'Group III, Row I'), ('III.II', 'Group III, Row II'),
          ('III.III', 'Group III, Row III'), ('III.V', 'Group III, Row V'), ('III.VI', 'Group III, Row VI')]
SHORT = {'4rho.a': r'$4\rho_a$', '4rho.b': r'$4\rho_b$', 'I.II.1': 'I.II.1', 'I.II.2': 'I.II.2', 'I.III': 'I.III',
         'I.IV': 'I.IV', 'I.V.1': 'I.V.1', 'I.V.2': 'I.V.2', 'II.I': 'II.I', 'III.I': 'III.I', 'III.II': 'III.II',
         'III.III': 'III.III', 'III.V': 'III.V', 'III.VI': 'III.VI'}


def kern_of(text):
    pairs = [tuple(s.replace(' ', '').split('-')) for s in KRE.findall(text)]
    assert len(pairs) == 4, text[:100]
    return [('KK',) + pairs[0] + pairs[1], ('KK',) + pairs[2] + pairs[3]]


def norm_kern(kern):
    out = []
    for f in kern:
        a, b, c, d = f[1:]
        out.append(tuple(sorted([(a, b), (c, d)])))
    return sorted(out)


def parse():
    raw = open(SRC).read()
    tex = '\n'.join(l for l in raw.split('\n') if not l.lstrip().startswith('%'))
    blocks = re.findall(r'\\begin\{align\}(.*?)\\end\{align\}', tex, re.S)
    assert len(blocks) == len(BLOCKS), len(blocks)
    out = []
    for (key, title), b in zip(BLOCKS, blocks):
        head, body = b.split(r'\Bigg[', 1)
        if r'-\,\frac{i}{(2\pi)^{3}}' in head:
            pref = -1j
        elif r'\frac{i}{(2\pi)^{3}}' in head or r'\frac{ig^{4}}{16\pi^{5}}' in head:
            pref = 1j
        else:
            raise ValueError(head)
        pos, n = 0, 0
        for m in LABEL.finditer(body):
            txt = body[pos:m.start()]
            pos = m.end()
            n += 1
            cls = m.group(2) + m.group(3)
            t = txt.replace(r'\notag', '').replace(r'\nonumber', '').replace('&', '').replace('\\\\', '').strip()
            t = re.sub(r'^(\\quad|\\qquad|\s)+', '', t)    # leading spacing commands
            t = re.sub(r'^\\-', '-', t)                     # '\-' typed for a minus sign
            sign = -1 if t.startswith('-') else 1
            t2 = re.sub(r'^(\\,)+', '', t.lstrip('+-').strip()).strip()
            c = 1.0
            if re.match(r'^2\s*\\,', t2):
                c *= 2
            if t2.startswith(r'\frac{1}{2}'):
                c *= 0.5
            anti = r'\big\{' in t
            if anti:
                c *= 2
            kern = kern_of(t)
            fs = [tuple(f) for f in FRE.findall(t)]
            rs = [(i, p) for i, p in RHO.findall(t)]
            Us = [(i, j, p) for i, j, p in URE.findall(t)]
            assert len(fs) == 1 and len(rs) == 3 and Us, (key, n, t[:80])
            out.append(dict(block=key, title=title, n=n, cls=cls, sign=sign, mult=c, pref=pref, anti=anti,
                            kern=kern, kern_ok=norm_kern(kern) == norm_kern(CLASS_KERN[cls]),
                            fs=fs, rs=rs, Us=Us, nU=len(Us), text=t))
    return out


def as_term(e):
    return term(e['pref'] * e['sign'] * e['mult'], e['Us'], e['fs'], e['rs'], e['kern'], '%s.%d' % (e['block'], e['n']))


if __name__ == '__main__':
    from collections import Counter
    E = parse()
    print('terms', len(E))
    print(Counter(e['block'] for e in E))
    print(Counter((e['cls'], e['nU']) for e in E))
    for e in E:
        if not e['kern_ok']:
            print('kernel differs from its label %s: %s.%d %s' % (e['cls'], e['block'], e['n'], e['kern']))
        idx = [x for f in e['fs'] for x in f] + [i for i, _ in e['rs']] + [x for i, j, _ in e['Us'] for x in (i, j)]
        cnt = Counter(idx)
        badidx = {k: v for k, v in cnt.items() if v != 2}
        if badidx:
            print('index count != 2 in %s.%d: %s' % (e['block'], e['n'], badidx))
