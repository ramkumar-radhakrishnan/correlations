r"""The equations of your current Group_I/II/III.tex against the versions that were checked against Calculation_main.tex
(../group1/input/Group_I_v2.tex, ../groups23/input/Group_II.tex, Group_III.tex), equation by equation.
Normalised: no colours, spaces, \, , \notag; the integration measures \int_{...} (lists of points) and a leading sign are
reported separately.  Output: eq_diff.json"""
import difflib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rows_check'))
import pe  # noqa: E402
pe.LINES = {g: open(os.path.join(HERE, '..', 'full_check', 'input', '%s.tex' % g)).read().split('\n') for g in ('Group_I', 'Group_II', 'Group_III')}
import eqnum as EN  # noqa: E402

OLD = {'Group_I': '../group1/input/Group_I_v2.tex', 'Group_II': '../groups23/input/Group_II.tex',
       'Group_III': '../groups23/input/Group_III.tex'}


def envs(path):
    L = open(os.path.join(HERE, path)).read().split('\n')
    out, i = [], 0
    while i < len(L):
        if r'\begin{align}' in L[i] and not L[i].lstrip().startswith('%'):
            j = i
            while r'\end{align}' not in L[j]:
                j += 1
            out.append((i + 1, j + 1, '\n'.join(L[i:j + 1])))
            i = j + 1
        else:
            i += 1
    return out


def colours(t):
    for c in ('blue', 'red', 'green', 'orange', 'violet', 'magenta', 'teal', 'olive', 'brown', 'purple', 'black'):
        while '{\\color{%s}' % c in t:
            i = t.index('{\\color{%s}' % c)
            d, k = 1, i + 1
            while d:
                d += {'{': 1, '}': -1}.get(t[k], 0)
                k += 1
            t = t[:i] + t[i + len('{\\color{%s}' % c):k - 1] + t[k:]
    return t


def norm(t):
    t = colours(t)
    t = re.sub(r'\\label\{[^}]*\}', '', t)
    t = re.sub(r'\s+|\\,|\\notag|\\nonumber|\\!|\;', '', t)
    t = t.replace('^\\prime', "'").replace('^{\\prime}', "'")
    meas = re.findall(r'\\int_\{[^{}]*\}', t)
    t = re.sub(r'\\int_\{([a-z,\']*)\}', r'\\int_{M}', t)
    t = re.sub(r'(&=)-', r'\1', t)
    return t, meas


out = {}
for g in ('Group_I', 'Group_II', 'Group_III'):
    new, old = envs('../full_check/input/%s.tex' % g), envs(OLD[g])
    on = [norm(t)[0] for _, _, t in old]
    res = []
    for a, b, t in new:
        nt, nm = norm(t)
        num = EN.eq(g, a, b)
        if nt in on:
            ot = old[on.index(nt)][2]
            sign = (re.search(r'&\s*=\s*-', colours(t)) is None) != (re.search(r'&\s*=\s*-', colours(ot)) is None)
            meas = nm != norm(ot)[1]
            res.append(dict(eq=num, lines=[a, b], status='same' + (' up to the overall sign' if sign else '') +
                            (' and the measure' if meas else '')))
        else:
            best = max(range(len(on)), key=lambda k: difflib.SequenceMatcher(None, on[k], nt, autojunk=False).ratio())
            sm = difflib.SequenceMatcher(None, on[best], nt, autojunk=False)
            ch = [(on[best][i1:i2][:120], nt[j1:j2][:120]) for tag, i1, i2, j1, j2 in sm.get_opcodes() if tag != 'equal']
            res.append(dict(eq=num, lines=[a, b], status='changed', closest_old_lines=old[best][:2],
                            ratio=round(sm.ratio(), 4), changes=ch[:6]))
    out[g] = res
    for r in res:
        print(g, r['eq'], r['lines'], r['status'], r.get('changes', '')[:3] if r['status'] == 'changed' else '')
json.dump(out, open(os.path.join(HERE, 'eq_diff.json'), 'w'), indent=1)
