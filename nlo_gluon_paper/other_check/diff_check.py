r"""Your Other_terms.tex (input/) against the remaining terms of ../full_check/Remaining_RC_DGLAP_finite.tex, equation by
equation, after removing markup (red lines dropped, blue kept, spacing and \left/\right ignored).  -> diff_check.txt"""
import difflib
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def strip_cmd(s, cmd, keep):
    out, i = '', 0
    while True:
        j = s.find('\\' + cmd + '{', i)
        if j < 0:
            return out + s[i:]
        out += s[i:j]
        k = j + len(cmd) + 2
        d, st = 1, k
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0)
            k += 1
        if keep:
            out += s[st:k - 1]
        i = k


def drop_rc(s):
    """remove {\rc ...} groups and lines starting with \rc (removed lines of the earlier file)"""
    s = re.sub(r'^.*\\rc\b.*$', '', s, flags=re.M)
    return s


def norm(s):
    s = strip_cmd(s, 'old', False)
    s = strip_cmd(s, 'new', True)
    s = re.sub(r'\\label\{[^}]*\}', '', s)
    s = re.sub(r'\\(nc)\b', '', s)
    s = re.sub(r'\{\\color\{\w+\}', '{', s)
    s = re.sub(r'\\notag|\\\\|&|\\,|\;|\\!|\\quad|\\qquad|\\left|\\right|\\Big|\\big|\\Bigg|\\bigg', '', s)
    s = re.sub(r'\s+', '', s)
    return s.replace('{', '').replace('}', '')


def eqs(t):
    return [norm(m.group(1)) for m in re.finditer(r'\\begin\{align\}(.*?)\\end\{align\}', t, re.S)]


u = eqs(open(os.path.join(HERE, 'input', 'Other_terms.tex')).read())
r = eqs(drop_rc(open(os.path.join(HERE, '..', 'full_check', 'Remaining_RC_DGLAP_finite.tex')).read()))
lines = []
for i, a in enumerate(u):
    best = max(range(len(r)), key=lambda j: difflib.SequenceMatcher(None, a[:2500], r[j][:2500], autojunk=False).ratio())
    sm = difflib.SequenceMatcher(None, a, r[best], autojunk=False)
    ops = [(op, a[i1:i2][:200], r[best][j1:j2][:200]) for op, i1, i2, j1, j2 in sm.get_opcodes() if op != 'equal']
    lines.append('yours #%d -> earlier #%d  ratio %.4f  %d differences' % (i, best, sm.ratio(), len(ops)))
    lines += ['    %s | yours: %s | earlier: %s' % o for o in ops]
open(os.path.join(HERE, 'diff_check.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
