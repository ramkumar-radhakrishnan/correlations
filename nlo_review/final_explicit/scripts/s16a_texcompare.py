# s16a_texcompare.py -- TeX-level comparison of the author's file Final_terms_explicit_2_author.tex (2026-11-23)
# with (i) the previous author version Final_terms_explicit_1_author.tex and (ii) the verified file
# Final_terms_explicit.tex (colour tensors C^{abc}_{row} substituted, WITH and WITHOUT enclosing brackets).
# Run: python3 -I s16a_texcompare.py -> s16a_texcompare.log
import os, re, difflib
HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE)
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
def strip_cmd(s, cmd):
    o = []; i = 0; key = '\\' + cmd + '{'
    while True:
        j = s.find(key, i)
        if j < 0: o.append(s[i:]); break
        o.append(s[i:j]); k = j + len(key); d = 1; st = k
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0); k += 1
        o.append(strip_cmd(s[st:k - 1], cmd)); i = k
    return ''.join(o)
def norm(s, keepbr=False):
    s = strip_cmd(s, 'ph'); s = re.sub(r'\\color\{[a-zA-Z]+\}', '', s)
    s = re.sub(r'\\(quad|qquad|notag|left|right|big|Big|bigg|Bigg)(?![a-zA-Z])', '', s); s = re.sub(r'\\[,;!]', '', s)
    s = s.replace('\\\\', '').replace('&', ''); s = re.sub(r'\s+', '', s)
    for _ in range(5): s = re.sub(r'\{\{([^{}]*)\}\}', r'{\1}', s)
    s = re.sub(r'\\cdot\{(\([^{}]*\))\}', r'\\cdot\1', s); s = s.replace("^{\\prime}", "'")
    s = s.replace('\\times', '')
    s = s.replace('{', '').replace('}', '')
    return s
def eqs(t):
    res = {}
    for m in re.finditer(r'\\begin\{align\}(.*?)\\end\{align\}', t, re.S):
        body = m.group(1); labs = re.findall(r'\\label\{([^}]*)\}', body); parts = re.split(r'\\label\{[^}]*\}', body)
        for l, p in zip(labs, parts): res[l] = p
    return res
def heads(s): return re.sub(r'^\\frac\{d\^\{3\}(\\Delta)?N\^\{\\rm[^}]*\}(_\{[^}]*\})?\}\{d\^\{3\}k\}(\\bigg\|[^=]*)?=', '', s)
ver = open(os.path.join(TOP, 'Final_terms_explicit.tex')).read()
a1 = open(os.path.join(TOP, 'Final_terms_explicit_1_author.tex')).read()
a2 = open(os.path.join(TOP, 'Final_terms_explicit_2_author.tex')).read()
defs = {m.group(1): m.group(2).strip() for m in re.finditer(r'(C\^\{abc\}_\{\{\\rm [^}]*\}\})&\\equiv (.*?)\\notag', ver)}
vs = ver; vb = ver
for k, v in defs.items():
    vs = vs.replace(k, v); vb = vb.replace(k, r'\Big[' + v + r'\Big]')
V, VB, A1, A2 = eqs(vs), eqs(vb), eqs(a1), eqs(a2)
def cmp(A, B, nameA, nameB):
    pr('== %s vs %s ==' % (nameA, nameB))
    for l in sorted(set(A) | set(B), key=lambda l: (a2.find('label{%s}' % l) if l in A else 10 ** 7)):
        if l not in A: pr('   ONLY in %s : %s' % (nameB, l)); continue
        if l not in B: pr('   ONLY in %s : %s' % (nameA, l)); continue
        a, b = norm(heads(A[l].strip())), norm(heads(B[l].strip()))
        if a == b: pr('   same : %s' % l); continue
        sm = difflib.SequenceMatcher(None, b, a, autojunk=False)
        d = ['%s %r -> %r' % (op, b[i1:i2][:70], a[j1:j2][:70]) for op, i1, i2, j1, j2 in sm.get_opcodes() if op != 'equal']
        pr('   DIFF : %s : %s' % (l, ' | '.join(d)[:900]))
cmp(A2, A1, 'author2', 'author1')
cmp(A2, V, 'author2', 'verified(C substituted, no brackets)')
cmp(A2, VB, 'author2', 'verified(C substituted IN brackets)')
open(os.path.join(HERE, 's16a_texcompare.log'), 'w').write('\n'.join(out) + '\n')
