r"""The split machinery of ../llsplit (blocks.py, gen.py, three_check.py, pole_check.py) on your new files (input/).

Changes for the new files:
  * Group_II.tex: one line fewer before Row I part 2 (a stray display was removed): the blocks after it move up by 1.
  * Group_I.tex Row V 2b UV: your bracket is now [2 log(k+/Lambda) + 2 log(1 - Lambda/k+) - 11/6] (E1 done).
  * 4rho.tex: log(vee/Lambda) is written [log(k+/Lambda) + log(vee/k+)].
Import this module before the llsplit modules are used."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
LL = os.path.normpath(os.path.join(HERE, '..', 'llsplit'))
for sub in ('llsplit', 'bfkl_rows', 'groups23'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.append(p)

import blocks as BL  # noqa: E402

BL.FILES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n')
            for g in ('Group_I', 'Group_II', 'Group_III', '4rho')}
for name, (g, a, b) in list(BL.FINAL.items()):
    if g == 'Group_II' and a > 52:
        BL.FINAL[name] = (g, a - 1, b - 1)
# every block must start and end its align environment
for name, (g, a, b) in BL.FINAL.items():
    assert r'\begin{align}' in BL.FILES[g][a - 1] and r'\end{align}' in BL.FILES[g][b - 1], (name, a, b)

import gen as GN  # noqa: E402

_special = GN.special


def special(name, n, t):
    if name == 'G1 V(2b UV)':
        t = GN.strip_user_colours(t)
        m = GN.rx(r'\left[2\,\log{\left(\frac{k^{+}}{\Lambda}\right)}+2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]').search(t)
        if m is None:                      # the earlier files (without E1): the old rule
            return _special(name, n, t)
        J = t[:m.start()] + r'2\,' + GN.LKL + t[m.end():]
        F = t[:m.start()] + r'\left[-\frac{11}{6}\right]' + t[m.end():]
        Fx = t[:m.start()] + r'\left[2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]' + t[m.end():]
        return dict(J=J, B=None, F=F, Fx=Fx)
    return _special(name, n, t)


GN.special = special

import three_check as T3  # noqa: E402

FOUR = '\n'.join(BL.FILES['4rho'])
LOGB = r"{\color{blue}\left[\log{\left[\frac{k^{+}}{\Lambda}\right]}+\log{\left[\frac{\vee}{k^{+}}\right]}\right]}"


def four_rho_three(part):
    """three-rho lists of 4rho.tex ('We rewrite them', 'Rewriting with some simplification'): their log(k+/Lambda)
    part ('J') or log(vee/k+) part ('B'); both logs have coefficient 1"""
    big = T3.LKL if part == 'J' else T3.LVK
    out = []
    for head in (r'We rewrite them (by some simplification)', r'Rewriting with some simplification, we get'):
        i = FOUR.index(head)
        t = FOUR[FOUR.index(r'\begin{align}', i) + len(r'\begin{align}'):FOUR.index(r'\end{align}', i)]
        t = t[re.search(r'&\s*=', t).end():]
        assert t.count(LOGB) >= 1, head
        t = t.replace(LOGB, big)
        out.append(('4rho.tex ' + head.split()[0] + ' ' + head.split()[1], t))
    return out


T3.four_rho_three = four_rho_three
