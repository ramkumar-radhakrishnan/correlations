"""Parse input/three_rho_all.tex (the three-rho terms proportional to log(vee/Lambda), with c.c., as sent).
Every term becomes a dict in the format of ../two_rho/core.py, coefficient in units of N, charges in the written order."""
import os
import re
import env  # noqa
from core import term

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'input', 'three_rho_all.tex')

LABEL = re.compile(r"\{\\color\{(\w+)\}\s*\\(alpha|beta|gamma)\s*('?)\s*\}")
KRE = re.compile(r"\\mathcal\{K\}\(([^)]*)\)")
PHI = re.compile(r"\\mathcal\{\\phi\}_\{(w'?)\}\((x'?)\)")
FRE = re.compile(r"f\^\{\s*(\w)\s*(\w)\s*(\w)\s*\}")
RHO = re.compile(r"\\rho\^\{(\w)\}\((\w'?)\)")
URE = re.compile(r"U\^\{\s*(\w)\s*(\w)\s*\}\((\w'?)\)")

# canonical kernel of every class (common variables)
CLASS_KERN = {
    'alpha': [('KK', "x'", "w'", 'x', 'w'), ('KK', 'x', 'z', 'y', 'z')],
    "alpha'": [('KK', "x'", "w'", 'x', 'w'), ('KK', "x'", 'z', 'y', 'z')],
    'beta': [('KK', 'x', 'w', "x'", "w'"), ('KK', 'w', 'z', 'y', 'z')],
    "beta'": [('KK', 'x', 'w', "x'", "w'"), ('KK', "w'", 'z', 'y', 'z')],
    'gamma': [('KK', "x'", "w'", 'y', 'z'), ('PHI', 'x', 'z', 'w')],
    "gamma'": [('KK', 'x', 'w', 'y', 'z'), ('PHI', "x'", 'z', "w'")],
}

# blocks in the order of the file: (key, title)
BLOCKS = [('4rho.a', r'From four $\rho$\'s, $d^3N_a$'), ('4rho.b', r'From four $\rho$\'s, $d^3N_b$'),
          ('I.II.1', 'Group I, Row II, part 1'), ('I.II.2', 'Group I, Row II, part 2'),
          ('I.III', 'Group I, Row III'), ('I.IV', 'Group I, Row IV'),
          ('I.V.1', 'Group I, Row V, part 1'), ('I.V.2', 'Group I, Row V, part 2'),
          ('II.I', 'Group II, Row I'), ('III.I', 'Group III, Row I'), ('III.II', 'Group III, Row II'),
          ('III.III', 'Group III, Row III'), ('III.V', 'Group III, Row V'), ('III.VI', 'Group III, Row VI')]


def kern_of(text):
    ks = KRE.findall(text)
    out = []
    pairs = [tuple(s.replace(' ', '').split('-')) for s in ks]
    for n in range(0, len(pairs), 2):
        (a, b), (c, d) = pairs[n], pairs[n + 1]
        out.append(('KK', a, b, c, d))
    for w, x in PHI.findall(text):
        out.append(('PHI', x, 'z', w))
    return out


def parse():
    tex = open(SRC).read()
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
        pos = 0
        n = 0
        for m in LABEL.finditer(body):
            txt = body[pos:m.start()]
            pos = m.end()
            n += 1
            cls = m.group(2) + m.group(3)
            t = txt.replace(r'\notag', '').replace(r'\nonumber', '').replace('&', '').replace('\\\\', '').strip()
            sign = -1 if t.startswith('-') else 1
            t2 = t.lstrip('+-').strip()
            mult = 2 if re.match(r'^2\s*\\,', t2) else 1
            kern = kern_of(t)
            fs = [tuple(f) for f in FRE.findall(t)]
            rs = [(i, p) for i, p in RHO.findall(t)]
            Us = [(i, j, p) for i, j, p in URE.findall(t)]
            assert len(fs) == 1 and len(rs) == 3 and Us, (key, n, t[:80])
            out.append(dict(block=key, title=title, n=n, cls=cls, sign=sign, mult=mult, pref=pref,
                            kern_written=kern, fs=fs, rs=rs, Us=Us, text=t))
    return out


def as_term(e, kern='class'):
    """kern='class': the kernel of the class label; 'written': the kernel as written."""
    k = CLASS_KERN[e['cls']] if kern == 'class' else e['kern_written']
    return term(e['pref'] * e['sign'] * e['mult'], e['Us'], e['fs'], e['rs'], k, '%s.%d' % (e['block'], e['n']))


if __name__ == '__main__':
    from collections import Counter
    E = parse()
    print('terms', len(E), Counter(e['block'] for e in E))
    import acc3  # noqa
    from acc3 import kkey
    def ck(kern):
        parts = []
        sgn = 1
        for f in kern:
            if f[0] == 'PHI':
                parts.append(('PHI', f[1], f[3]))
            else:
                k, s = __import__('core').kkey(*f[1:]) if False else (None, None)
        return None
    for e in E:
        w = e['kern_written']
        c = CLASS_KERN[e['cls']]
        def norm(kern):
            out = []
            for f in kern:
                if f[0] == 'KK':
                    a, b, c_, d = f[1:]
                    p1, p2 = (a, b), (c_, d)
                    out.append(('KK',) + tuple(sorted([p1, p2])))
                else:
                    out.append(f)
            return sorted(out)
        if norm(w) != norm(c):
            print('kernel differs from class %-7s %s.%d: written %s' % (e['cls'], e['block'], e['n'], w))
