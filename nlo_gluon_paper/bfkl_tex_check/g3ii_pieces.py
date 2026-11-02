import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_tex as ET
C, T3, R, TE = ET.C, ET.T3, ET.R, ET.TE
from rows_user import rows
OLD = {(r[0].split('\\')[0].strip(), r[1], r[2]): r[4] for r in rows(os.path.join(HERE, '..', 'bfkl_final', 'input', 'bfkl.tex'))}
old = C.prep(OLD[('Group III', 'Row II', None)])
i = old.index(r'\notag \\ & +\frac{1}{(2\pi)^{3}}')
p1, p2 = old[:i], old[i + len(r'\notag \\ & '):]
c = T3.configs(5)
env = T3.env(c)
v1, v2 = TE.term_value(p1, env), TE.term_value(p2, env)
ts = [t for t in R.THREE['G3 II']() if not t['tag'].endswith('*')]
ve = sum(t['c'] * TE.term_value(ET.tex(t), env) for t in ts)
print('old line1', v1, 'line2', v2, 'sum', v1 + v2, 'engine', ve)
print(TE.translate(p1)[:300])
import numpy as np
tcc = [t for t in R.THREE['G3 II']() if t['tag'].endswith('*')]
vecc = sum(t['c'] * TE.term_value(ET.tex(t), env) for t in tcc)
vocc = np.conj(TE.term_value(old, T3.env(T3.swapped_w(c))))
print('c.c.: engine', vecc, 'mine', vocc)
