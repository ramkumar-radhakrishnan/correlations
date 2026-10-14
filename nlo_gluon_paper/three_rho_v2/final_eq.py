"""One equation with all three-rho terms proportional to log(vee/Lambda), final form, in the notation of the list as sent.
black = as sent; red = removed (F3, F4); blue = changed (T1-T3, operator order) or added (D).
Writes ThreeRho_final.tex; then re-parses the equation it wrote and checks it (three-rho cancellation, two-rho terms)."""
import re
import sys
import env
from core import term, symmetrize, Config, pointwise_total
import parse as P
import fixes as FX
import oporder as OP

KTEX = {
    'alpha': r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(y-z)",
    "alpha'": r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(x'-z)\cdot\mathcal{K}(y-z)",
    'beta': r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(w-z)\cdot\mathcal{K}(y-z)",
    "beta'": r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(w'-z)\cdot\mathcal{K}(y-z)",
}
LAB = {'alpha': r'\alpha', "alpha'": r"\alpha'", 'beta': r'\beta', "beta'": r"\beta'", 'gamma': r'\gamma', "gamma'": r"\gamma'"}
HEAD = [('4rho.a', r"From four $\rho$'s, $d^3N_a$ (its overall minus sign is absorbed in the signs of its terms)"),
        ('4rho.b', r"From four $\rho$'s, $d^3N_b$"),
        ('I.II.1', 'Group I, Row II, part 1'), ('I.II.2', 'Group I, Row II, part 2 (c.c.)'),
        ('I.III', 'Group I, Row III'), ('D', r"Group I, Row III: $\mathcal D$, replaces the $\gamma,\gamma'$ terms"),
        ('I.IV', 'Group I, Row IV'), ('I.V.1', 'Group I, Row V, part 1'), ('I.V.2', 'Group I, Row V, part 2'),
        ('II.I', 'Group II, Row I'), ('III.I', 'Group III, Row I'), ('III.II', 'Group III, Row II'),
        ('III.III', 'Group III, Row III'), ('III.V', 'Group III, Row V'), ('III.VI', 'Group III, Row VI')]


def rho(i, p):
    return r'\rho^{%s}(%s)' % (i, p)


def coef(c):
    if c == 1:
        return '+'
    if c == -1:
        return '-'
    if c == 0.5:
        return r'+\tfrac12\,'
    if c == -0.5:
        return r'-\tfrac12\,'
    return ('+%d\\,' % c) if c > 0 else ('%d\\,' % c)


def flip(text):
    t = text.strip()
    if t.startswith('-'):
        r = t[1:].lstrip()
        if r.startswith('\\,'):
            r = r[2:].lstrip()
        return '+' + r
    return '-\\,' + t.lstrip('+').lstrip()


def generated(e, of, sgn):
    byp = {p: (i, p) for i, p in e['rs']}
    o = [byp[p] for p in of['order']]
    c = sgn * e['sign'] * e['mult']
    if of['form'] == 'ABC':
        rp = r'\,'.join(rho(*r) for r in o)
    elif of['form'] == 'A{BC}':
        c = c / 2
        rp = r'%s\big\{%s,%s\big\}' % (rho(*o[0]), rho(*o[1]), rho(*o[2]))
    else:
        c = c / 2
        rp = r'\big\{%s,%s\big\}%s' % (rho(*o[0]), rho(*o[1]), rho(*o[2]))
    us = r'\,'.join(r'U^{%s%s}(%s)' % u for u in e['Us'])
    return r'%s%s\, f^{%s%s%s}\,%s\,%s' % (coef(c), KTEX[e['cls']], *e['fs'][0], rp, us)


def lines():
    sent = {(e['block'], e['n']): e for e in P.parse()}
    fixed = FX.typo_fixed()
    D = FX.d_entries()
    out = []
    for blk, title in HEAD:
        out.append(('head', title))
        es = D if blk == 'D' else [e for e in fixed if e['block'] == blk]
        sgn = -1 if blk == '4rho.a' else 1
        for e in es:
            num = ('D%d' % e['n']) if blk == 'D' else str(e['n'])
            lab = LAB[e['cls']]
            if (e['block'], e['n']) in FX.F4_DROP:
                out.append(('red', num, sent[(blk, e['n'])]['text'], lab, 'removed, F4'))
                continue
            if blk == 'I.III' and e['cls'] in ('gamma', "gamma'"):
                out.append(('red', num, sent[(blk, e['n'])]['text'], lab, 'removed, F3'))
                continue
            of = OP.operator_form(e)
            reordered = of['form'] != 'ABC' or [p for _, p in e['rs']] != of['order']
            tags = []
            if e.get('fix'):
                tags.append(e['fix'])
            if blk == 'D':
                tags.append('added, F3')
            elif reordered:
                tags.append('c.c.' if (of['cc'] and of['form'] == 'ABC') else 'order')
            if tags:
                out.append(('blue', num, generated(e, of, sgn), lab, ', '.join(tags)))
            else:
                t = e['text'] if sgn == 1 else flip(e['text'])
                out.append(('black', num, t, lab, ''))
    return out


def write(path):
    L = lines()
    eq = [r'\begin{align}',
          r"\frac{d^{3}N^{3\rho}_{\rm NLO}}{d^{3}k}={}&\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log\Big[\frac{\vee}{\Lambda}\Big]\int_{x,y,z,x',w,w'}e^{-ik(w'-w)}\,\Bigg[\notag\\"]
    body = []
    for item in L:
        if item[0] == 'head':
            body.append(r'&\quad\underline{\text{%s}}\notag\\' % item[1])
            continue
        kind, num, txt, lab, tag = item
        content = r'%s\;(%s)' % (txt.strip().rstrip('\\,').rstrip(), lab)
        tagtex = (r'\quad\text{\scriptsize[%s]}' % tag) if tag else ''
        if kind == 'red':
            body.append(r'&\quad\text{(%s)}\ \ \textcolor{red}{%s}%s\notag\\' % (num, content, tagtex))
        elif kind == 'blue':
            body.append(r'&\quad\text{(%s)}\ \ \textcolor{blue}{%s}%s\notag\\' % (num, content, tagtex))
        else:
            body.append(r'&\quad\text{(%s)}\ \ %s\notag\\' % (num, content))
    body[-1] = body[-1].replace(r'\notag\\', r'\Bigg]\label{eq:three}')
    eq += body + [r'\end{align}']
    doc = open('final_static.tex').read().replace('@@EQ@@', '\n'.join(eq))
    open(path, 'w').write(doc)
    return L


# ------------------------------------------------------------------ re-parse the written equation and check it
TOK = re.compile(r"\\big\\\{|\\big\\\}|\\rho\^\{(\w)\}\((\w'?)\)")
LINE = re.compile(r"^&\\quad\\text\{\(([^)]*)\)\}\\ \\ (.*)\\notag\\\\$|^&\\quad\\text\{\(([^)]*)\)\}\\ \\ (.*)\\Bigg\]\\label\{eq:three\}$")
CLASS_OF_LAB = {v: k for k, v in LAB.items()}


def reparse(path):
    tex = open(path).read()
    eq = tex[tex.index(r'\begin{align}'):tex.index(r'\end{align}')]
    out, blk = [], None
    hmap = {t: b for b, t in HEAD}
    for ln in eq.split('\n'):
        m = re.match(r'^&\\quad\\underline\{\\text\{(.*)\}\}\\notag\\\\$', ln)
        if m:
            blk = hmap[m.group(1)]
            continue
        m = LINE.match(ln)
        if not m:
            continue
        num, body = (m.group(1), m.group(2)) if m.group(1) is not None else (m.group(3), m.group(4))
        color = 'red' if r'\textcolor{red}' in body else ('blue' if r'\textcolor{blue}' in body else 'black')
        body = re.sub(r'\\quad\\text\{\\scriptsize\[[^]]*\]\}', '', body)
        body = body.replace(r'\textcolor{red}{', '').replace(r'\textcolor{blue}{', '')
        lab = re.search(r"\\;\((\\(?:alpha|beta|gamma)'?)\)", body).group(1)
        t = body.strip()
        sign = -1 if t.startswith('-') else 1
        t2 = t.lstrip('+-').strip()
        if t2.startswith('\\,'):
            t2 = t2[2:].strip()
        mult = 2 if re.match(r'^2\s*\\,', t2) else (0.5 if t2.startswith(r'\tfrac12') else 1)
        toks = [(mm.group(0), mm.group(1), mm.group(2)) for mm in TOK.finditer(t)]
        rs = [(i, p) for tk, i, p in toks if i]
        if toks[0][0].startswith(r'\big\{'):
            form = '{BC}A'
        elif len(toks) > 1 and toks[1][0].startswith(r'\big\{'):
            form = 'A{BC}'
        else:
            form = 'ABC'
        out.append(dict(block=blk, num=num, color=color, cls=CLASS_OF_LAB[lab], sign=sign, mult=mult, form=form,
                        kern=P.kern_of(t), fs=[tuple(f) for f in P.FRE.findall(t)], rs=rs,
                        Us=[tuple(u) for u in P.URE.findall(t)]))
    return out


def terms_of(r):
    c = 1j * r['sign'] * r['mult']
    rs = r['rs']
    if r['form'] == 'ABC':
        words = [rs]
    elif r['form'] == 'A{BC}':
        words = [rs, [rs[0], rs[2], rs[1]]]
    else:
        words = [rs, [rs[1], rs[0], rs[2]]]
    return [term(c, r['Us'], r['fs'], w, r['kern'], '%s.%s' % (r['block'], r['num'])) for w in words]


def check(path):
    from collections import defaultdict
    import acc3 as A3
    import target as TG
    sys.path.insert(0, env.ALL2)
    import sets as S2
    R = reparse(path)
    final = [r for r in R if r['color'] != 'red']
    print('re-parsed lines: %d (black %d, blue %d, red %d)' % (len(R), sum(r['color'] == 'black' for r in R),
                                                              sum(r['color'] == 'blue' for r in R), sum(r['color'] == 'red' for r in R)))
    res = {}
    for seed in (5, 17, 23):
        acc = defaultdict(lambda: defaultdict(complex))
        cfg = A3.Config3(seed, True)
        for r in final:
            for t in terms_of(r):
                A3.accumulate3(symmetrize(t, (3,)), cfg, acc, label=r['block'])
        bad = A3.report3(acc, show=False)
        print('three-rho, seed %d: kernel keys %d, non-cancelling %d' % (seed, len(acc), len(bad)))
        res['three_%d' % seed] = len(bad)
    two = [s for r in final for t in terms_of(r) for s in symmetrize(t, (2,))]
    A = S2.A_sets()
    other = [t for k in ('g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III', '4rho') for t in A[k]]
    earlier = [t for k in S2.THREE_LABELS + ['A,B'] for t in A[k]]
    J = TG.jimwlk()
    for seed in (3, 7, 11):
        cfg = Config(seed, True)
        a = pointwise_total(two, cfg)[0].real
        b = pointwise_total(earlier, cfg)[0].real
        tot = a + pointwise_total(other, cfg)[0].real
        j = pointwise_total(J, cfg)[0].real
        print('two-rho, seed %d: from this equation %.8f, earlier %.8f; + genuine + 4rho = %.8f, JIMWLK x LO = %.8f' % (seed, a, b, tot, j))
        res['two_%d' % seed] = (a, b, tot, j)
    return res


if __name__ == '__main__':
    write('ThreeRho_final.tex')
    check('ThreeRho_final.tex')
