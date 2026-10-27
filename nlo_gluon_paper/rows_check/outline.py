"""Outline of an input file: every align environment with its line range and the start of its right-hand side."""
import re
import sys
g = sys.argv[1]
L = open('input/%s.tex' % g).read().split('\n')
i = 0
while i < len(L):
    l = L[i]
    if r'\begin{align}' in l and not l.strip().startswith('%'):
        j = i
        while r'\end{align}' not in L[j]:
            j += 1
        txt = ' '.join(L[i:j + 1])
        m = re.search(r'&\s*=', txt)
        rhs = txt[m.end():m.end() + 130] if m else txt[:130]
        print('%3d-%3d ALIGN %s' % (i + 1, j + 1, re.sub(r'\s+', ' ', rhs)))
        i = j
    elif l.strip() and not l.strip().startswith('%') and i > 30 and not l.strip().startswith(r'\end{document}'):
        print('%3d     %s' % (i + 1, l.strip()[:130]))
    i += 1
