"""Your new Finite_terms.tex against the remaining terms of ../bfkl_split/Finite_terms_split.tex (Secs. 7-9),
equation by equation (whitespace, \\notag, &, \\, removed; \\new{} kept, \\old{} dropped).
The remaining terms were built by exact substitutions, each checked numerically (../bfkl_split/split_logs.json).
Output: compare_file.json"""
import difflib
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
U = open(os.path.join(HERE, 'input', 'Finite_terms.tex')).read()
S = open(os.path.join(HERE, '..', 'bfkl_split', 'Finite_terms_split.tex')).read()
S = S[S.index(r'\section{Remaining terms: Group I'):]


def strip_cmd(t, cmd, keep):
    while cmd + '{' in t:
        i = t.index(cmd + '{')
        j = i + len(cmd) + 1
        d, k = 1, j
        while d:
            d += {'{': 1, '}': -1}.get(t[k], 0)
            k += 1
        t = t[:i] + (t[j:k - 1] if keep else '') + t[k:]
    return t


S = strip_cmd(S, r'\old', False)
S = strip_cmd(S, r'\new', True).replace(r'\nc', '').replace(r'\rc', '')
while r'{\color{red}' in S:                                   # {\color{red} removed text}
    i = S.index(r'{\color{red}')
    d, k = 1, i + 1
    while d:
        d += {'{': 1, '}': -1}.get(S[k], 0)
        k += 1
    S = S[:i] + S[k:]
blocks = lambda t: re.findall(r'\\begin\{align\}(.*?)\\end\{align\}', t, re.S)
norm = lambda b: re.sub(r'\\notag|\\\\|&|\\,|\\!|\\quad|\\qquad|\s+', '', b)
NAMES = ['G1 II(1)', 'G1 II(2b)', 'G1 II R^i', 'G1 III N1', 'G1 III N2', 'G1 III N3a', 'G1 III N3b', 'G1 III A',
         'G1 IV', 'G1 V(1)', 'G1 V(2a)', 'G1 V(2b UV)', 'G1 V(2b fin)', 'G2 I(1)', 'G2 I(2a)', 'G2 I(2b)', 'G2 I(2c)',
         'G2 I(2d)', 'G2 I(2e)', 'G2 I(2f)', 'G2 IV', 'G3 I', 'G3 II', 'G3 III(1)', 'G3 III(2)', 'G3 V', 'G3 VI']

if __name__ == '__main__':
    bu, bs = blocks(U), blocks(S)
    assert len(bu) == len(bs) == len(NAMES), (len(bu), len(bs))
    out = {}
    for name, a, b in zip(NAMES, bu, bs):
        na, nb = norm(a), norm(b)
        if na == nb:
            out[name] = 'identical'
            continue
        diffs = []
        for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, na, nb, autojunk=False).get_opcodes():
            if op != 'equal':
                diffs.append(dict(op=op, yours=na[i1:i2][:200], mine=nb[j1:j2][:200]))
        out[name] = diffs
    for k, v in out.items():
        print('%-12s %s' % (k, v if v == 'identical' else '%d difference(s): ' % len(v) + '; '.join(
            ('yours "%s" / mine "%s"' % (d['yours'][:60], d['mine'][:60])) for d in v)))
    json.dump(out, open(os.path.join(HERE, 'compare_file.json'), 'w'), indent=1)
