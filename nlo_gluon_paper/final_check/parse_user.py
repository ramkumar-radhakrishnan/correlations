"""Parser of the two-rho lists of Evolution_terms.tex (Sec. 2, 'Symmetrizing three and four rho's').
Every term is  coef  K(a-b).K(c-d) [1/(e-z)^2 | K(e-z).K(f-z)] [f f] U..U rho rho, in units of
N = (1/(2pi)^3) g^4/(16 pi^5) (1/k+) log(vee/Lambda).  Returns (label, term) with the term format of ../two_rho/core.py."""
import re
import env  # noqa
from core import term

KRE = re.compile(r"\\mathcal\{K\}\(([^)]*)\)")
INV = re.compile(r"\\frac\{1\}\{\((\w'?)-z\)\^\{2\}\}")
FRE = re.compile(r"f\^\{(\w)(\w)(\w)\}")
URE = re.compile(r"U\^\{(\w)(\w)\}\((\w'?)\)")
RRE = re.compile(r"\\rho\^\{(\w)\}\((\w'?)\)")
COEF = re.compile(r"^\s*([+-])?\s*(?:\\,)?\s*(?:\\frac\{(\d+)\}\{(\d+)\}|(\d+))?\s*(?:\\,)?\s*(N_\{c\})?\s*(?:\\,)?\s*(?=\\mathcal\{K\})")
NC = 3


def parse_chunk(ch, label):
    ch = ch.replace(r'\Bigg[', '').replace(r'\Bigg]', '').strip()
    m = COEF.match(ch)
    assert m, ch[:120]
    sign = -1 if m.group(1) == '-' else 1
    a = (int(m.group(2)) / int(m.group(3))) if m.group(2) else (int(m.group(4)) if m.group(4) else 1)
    ncp = 1 if m.group(5) else 0
    ks = [tuple(s.split('-')) for s in KRE.findall(ch)]
    assert len(ks) in (2, 4), ch
    kern = [('KK', ks[0][0], ks[0][1], ks[1][0], ks[1][1])]
    if len(ks) == 4:
        kern.append(('KK', ks[2][0], ks[2][1], ks[3][0], ks[3][1]))
    for p in INV.findall(ch):
        kern.append(('KK', p, 'z', p, 'z'))
    assert len(kern) == 2, ch
    Us, fs, rs = URE.findall(ch), FRE.findall(ch), RRE.findall(ch)
    assert len(rs) == 2, ch
    return term(sign * a, Us, fs, rs, kern, label, ncp=ncp)


def sec2_terms(path):
    tex = open(path).read()
    tex = tex[tex.index(r'\section{Symmetrizing three and four'):]
    out = []
    head = {'sub': '', 'subsub': ''}
    pos = 0
    pat = re.compile(r"\\subsection\{([^}]*)\}|\\subsubsection\{([^}]*)\}|\\textbf\{([^}]*)\}|\\begin\{align\}(.*?)\\end\{align\}", re.S)
    neq = {}
    for m in pat.finditer(tex):
        if m.group(1) is not None:
            head['sub'], head['subsub'] = m.group(1), ''
        elif m.group(2) is not None:
            head['subsub'] = m.group(2)
        elif m.group(4) is not None:
            body = m.group(4)
            name = re.search(r"\\frac\{d\^\{3\}N_\{?([^}]*?)\}?\}\{d\^\{3\}k\}", body)
            part = 'N1' if 'N_{1' in body[:80] else ('N2' if 'N_{2' in body[:80] else '')
            lab = '%s|%s|%s' % (head['sub'], head['subsub'], part)
            chunks = re.split(r"\\notag|\\\\|&", body)
            for ch in chunks:
                if r'\rho^' in ch and r'\mathcal{K}' in ch:
                    out.append((lab, parse_chunk(ch, lab)))
    return out


if __name__ == '__main__':
    from collections import Counter
    T = sec2_terms('input/Evolution_terms.tex')
    print(len(T))
    for k, v in Counter(l for l, _ in T).items():
        print('%-45s %d' % (k, v))
