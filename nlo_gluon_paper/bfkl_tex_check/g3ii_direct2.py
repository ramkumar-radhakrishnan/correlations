"""Group III Row II: engine vs old LaTeX row, classically, at fixed labels (no relabeling frame), c.c. included."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_tex as ET
C, T3, R, TE = ET.C, ET.T3, ET.R, ET.TE
from rows_user import rows
OLD = {(r[0].split('\\')[0].strip(), r[1], r[2]): r[4] for r in rows(os.path.join(HERE, '..', 'bfkl_final', 'input', 'bfkl.tex'))}
old = C.prep(OLD[('Group III', 'Row II', None)])
ts = R.THREE['G3 II']()
import itertools
for s in (5, 6):
    c = T3.configs(s)
    tot_e = tot_o = 0
    for perm in itertools.permutations(T3.CH):
        cp = T3.permuted(c, perm)
        tot_e += sum(t['c'] * 3 ** t.get('ncp', 0) * TE.term_value(ET.tex(t), T3.env(cp)) for t in ts)
        tot_o += TE.term_value(old, T3.env(cp)) + np.conj(TE.term_value(old, T3.env(T3.swapped_w(cp))))
    print(s, tot_e, tot_o, tot_e / tot_o)
