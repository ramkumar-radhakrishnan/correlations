r"""Your Evolution_two_rho.tex and Evolution_three_rho.tex (input/): the explicit equations, row by row.
rows(fname) -> [(group, row, part, first_line, rhs_text)]  (amplitude-level lines with \mathbb are skipped)."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))


def rows(fname, stop=None):
    base = os.environ.get('EVOL_DIR', os.path.join(HERE, 'input'))   # EVOL_DIR=corrected: check the corrected copies
    L = open(fname if os.path.isabs(fname) else os.path.join(base, fname)).read().split('\n')
    g = r = p = None
    out, i = [], 0
    while i < len(L):
        s = L[i]
        if stop and stop in s:
            break
        m = re.match(r'\s*\\section\{(.*)\}', s)
        if m:
            g = m.group(1); r = p = None
        m = re.match(r'\s*\\subsection\{(.*)\}', s)
        if m:
            r = m.group(1); p = None
        m = re.match(r'\s*\\subsubsection\{(.*)\}', s)
        if m:
            p = m.group(1)
        if r'\begin{align}' in s and not s.lstrip().startswith('%'):
            j = i
            while r'\end{align}' not in L[j]:
                j += 1
            t = '\n'.join(L[i:j + 1])
            t = t[t.index(r'\begin{align}') + 13:t.rindex(r'\end{align}')]
            t = re.sub(r'\\label\{[^}]*\}', '', t)
            if r'\mf' not in t:
                rhs = t[re.search(r'&\s*=|=\s*&', t).end():]
                out.append((g, r, p, i + 1, rhs))
            i = j + 1
        else:
            i += 1
    return out


if __name__ == '__main__':
    for f in ('Evolution_two_rho.tex', 'Evolution_three_rho.tex'):
        for e in rows(f, stop=r'\section{All three'):
            print(f[10:13], e[:4])
