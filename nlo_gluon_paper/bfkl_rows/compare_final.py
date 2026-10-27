r"""Your (commented-out) 'BFKL contribution' of BFKL.tex, both forms (term list and collected), against the sum of all
two-rho log(vee/k+) terms that was checked against BFKL x LO (../final_check/Final_evolution_check.tex, Eq. Bsum).
Pointwise, charges symmetrized, with c.c.  Output: compare_final.json"""
import json
import os
import re
import compare_4rho as C4
import rows as R
import texeval as TE

L = open(os.path.join(R.HERE, 'input', 'BFKL.tex')).read()
L = L[L.index(r'\section{BFKL contribution}'):L.index(r'In color singlet form')]
FORMS = []
for m in re.finditer(r'\\begin\{align\}(.*?)\\end\{align\}', L, re.S):
    t = m.group(1)
    FORMS.append(t[re.search(r'&\s*=', t).end():])

if __name__ == '__main__':
    ref = C4.fec('eq:Bsum')
    out = {}
    for seed in (1, 2, 3):
        c = TE.cfg(seed)
        r = C4.sym(ref, c)
        vals = [C4.sym(t, c) for t in FORMS]
        out[str(seed)] = [float(abs(v - r) / abs(r)) for v in vals]
        print(seed, out[str(seed)])
    json.dump(out, open(os.path.join(R.HERE, 'compare_final.json'), 'w'), indent=1)
