"""All logarithms in the four files (outside comments), with the line where they occur."""
import re, collections, json
FILES = ['Group_I', 'Group_II', 'Group_III', '4rho']
def strip_comments(t):
    return '\n'.join(l.split('%')[0] if not l.lstrip().startswith('%') else '' for l in t.split('\n'))
def grab(t, i):
    """argument of \\log at position i: {..}, \\left(..\\right), \\left[..\\right], \\Big(..\\Big), (..)"""
    j = i + 4
    while t[j] == ' ' or t[j:j+2] == '\\,':
        j += 1 if t[j] == ' ' else 2
    if t[j:j+3] == '^{2':
        j += 4
    if t[j] == '{':
        d, k = 1, j + 1
        while d:
            d += {'{': 1, '}': -1}.get(t[k], 0); k += 1
        return t[i:k]
    for o, c in ((r'\left(', r'\right)'), (r'\left[', r'\right]'), (r'\Big(', r'\Big)'), (r'\big(', r'\big)'), ('(', ')')):
        if t.startswith(o, j):
            d, k = 1, j + len(o)
            while d:
                if t.startswith(o, k): d += 1; k += len(o); continue
                if t.startswith(c, k): d -= 1; k += len(c); continue
                k += 1
            return t[i:k]
    return t[i:j + 30]
out = collections.OrderedDict()
for f in FILES:
    t = strip_comments(open('input/%s.tex' % f).read())
    for m in re.finditer(r'\\log', t):
        s = re.sub(r'\s+', '', grab(t, m.start()))
        line = t[:m.start()].count('\n') + 1
        out.setdefault(s, []).append('%s:%d' % (f, line))
for s, where in out.items():
    print(len(where), s[:160], ' @', ','.join(sorted(set(where)))[:120])
json.dump(out, open('logscan.json', 'w'), indent=1)
