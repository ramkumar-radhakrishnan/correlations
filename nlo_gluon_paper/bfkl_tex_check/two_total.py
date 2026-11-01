r"""Every two-rho list of your BFKL.tex (the 'written as' forms, the 'From 4 rho's' lists and the 'BFKL contribution'),
parsed with ../check_v4/parse_v4.py, evaluated pointwise (random points, random classical charges; U = U^T and general
U).  Compared: each list / their sum against BFKL x LO (../bfkl/refP.py) and against your 'BFKL contribution'.
Output: two_total.json"""
import json, os, re, sys
from collections import OrderedDict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'check_v4'))
import env  # noqa
from core import Config, pointwise_total
import refP
import sets as S2
import parse_v4 as P

SRC = os.path.join(HERE, 'input', 'BFKL.tex')
SEEDS = (3, 7, 11)


def equations(path=SRC):
    """parse_v4.equations, but an align that is not a two-rho list (three-rho first forms) is kept as kind 'other'"""
    tex = open(path).read()
    tex = tex[tex.index(r'\begin{document}'):]
    tex = re.sub(r'\\begin\{comment\}.*?\\end\{comment\}', '', tex, flags=re.S)
    h = dict(sec='', sub='', subsub='')
    eqs, singlet_next = [], False
    line_of = lambda pos: tex[:pos].count('\n')
    for m in P.HEAD.finditer(tex):
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
            e = dict(h, n=len(eqs), pos=m.start(), kind='other', terms=[], text=body)
            if singlet_next:
                e['kind'] = 'singlet'
                singlet_next = False
            elif r'\Bigg[' in body and r'\Bigg\{' not in body and r'\bigg[' not in body:
                inner = body[body.index(r'\Bigg[') + 6:body.rindex(r'\Bigg]')]
                try:
                    e['terms'] = P.parse_body(inner, '%s|%s|%s|%d' % (h['sec'], h['sub'], h['subsub'], len(eqs)))
                    e['kind'] = 'terms'
                except AssertionError as ex:
                    e['kind'] = 'parse error: %s' % str(ex)[:120]
            eqs.append(e)
    return eqs


def rel(a, b, sym, seeds=SEEDS):
    r = []
    for s in seeds:
        cfg = Config(s, sym)
        va, sa = pointwise_total(a, cfg)
        vb, sb = pointwise_total(b, cfg)
        r.append(abs(va - vb) / max(sa + sb, 1e-300))
    return max(r)


def both(a, b):
    return dict(sym=rel(a, b, True), general=rel(a, b, False))


if __name__ == '__main__':
    E = equations()
    out = OrderedDict()
    out['equations'] = [dict(sec=e['sec'], sub=e['sub'], subsub=e['subsub'], kind=e['kind'], n_terms=len(e['terms']))
                        for e in E]
    for e in out['equations']:
        print(e)
    lists = [e for e in E if e['kind'] == 'terms' and not e['sec'].startswith('BFKL')]
    final = [e for e in E if e['kind'] == 'terms' and e['sec'].startswith('BFKL')]
    tot = [t for e in lists for t in e['terms']]
    B = refP.bfkl()
    A = S2.A_sets()
    four = [t for e in lists if e['sec'].startswith('From') for t in e['terms']]
    out['From 4 rho two-rho lists vs 4rho.tex lists'] = both(four, A['4rho'])
    out['sum of all two-rho lists vs BFKLxLO'] = both(tot, B)
    for k, e in enumerate(final):
        out['BFKL contribution eq %d vs BFKLxLO' % k] = both(e['terms'], B)
        out['BFKL contribution eq %d vs sum of lists' % k] = both(e['terms'], tot)
    for k, v in out.items():
        if k != 'equations':
            print(k, v)
    json.dump(out, open(os.path.join(HERE, 'two_total.json'), 'w'), indent=1, default=str)
