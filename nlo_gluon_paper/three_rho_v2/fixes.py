"""The list as sent, the typo fixes T1-T3 and the fixes F3, F4.
T1: Group III Row I, term 15: rho^h(x') -> rho^c(x')
T2: Group III Row II, term 13: U^{hl}(y) -> U^{hl}(x)
T3: Group I Row II, part 1, alpha terms: KK(x'-w).KK(y-z) -> KK(x-z).KK(y-z)  (applied by using the class kernel)
F4: Group I Row II: remove the alpha terms of the U(w)A^(3) sub-terms (part 1 terms 1,2,5,6; part 2 terms 1,2,5,6)
F3: Group I Row III: remove the gamma, gamma' terms (part 2, edge kernel phi); add D (input/D_terms.tex)."""
import os
import re
import parse as P

HERE = os.path.dirname(os.path.abspath(__file__))
F4_DROP = {('I.II.1', n) for n in (1, 2, 5, 6)} | {('I.II.2', n) for n in (1, 2, 5, 6)}


def d_entries():
    tex = open(os.path.join(HERE, 'input', 'D_terms.tex')).read()
    body = tex.split(r'\Bigg[', 1)[1]
    out = []
    pos, n = 0, 0
    for m in P.LABEL.finditer(body):
        txt = body[pos:m.start()]
        pos = m.end()
        n += 1
        t = txt.replace(r'\notag', '').replace('&', '').replace('\\\\', '').strip()
        sign = -1 if t.startswith('-') else 1
        mult = 2 if re.match(r'^2\s*\\,', t.lstrip('+-').strip()) else 1
        out.append(dict(block='D', title=r'$\mathcal D$', n=n, cls=m.group(2) + m.group(3), sign=sign, mult=mult, pref=1j,
                        kern_written=P.kern_of(t), fs=[tuple(f) for f in P.FRE.findall(t)], rs=P.RHO.findall(t),
                        Us=P.URE.findall(t), text=t))
    return out


def typo_fixed():
    E = P.parse()
    for e in E:
        if e['block'] == 'III.I' and e['n'] == 15:
            e['rs'] = [('l', 'x'), ('e', 'y'), ('c', "x'")]
            e['fix'] = 'T1'
        if e['block'] == 'III.II' and e['n'] == 13:
            e['Us'] = [('d', 'e', 'y'), ('h', 'l', 'x'), ('b', 'd', 'z'), ('h', 'a', "w'")]
            e['fix'] = 'T2'
        if e['block'] == 'I.II.1' and e['cls'] == 'alpha':
            e['fix'] = 'T3'
    return E


def corrected():
    E = [e for e in typo_fixed() if (e['block'], e['n']) not in F4_DROP
         and not (e['block'] == 'I.III' and e['cls'] in ('gamma', "gamma'"))]
    return E + d_entries()
