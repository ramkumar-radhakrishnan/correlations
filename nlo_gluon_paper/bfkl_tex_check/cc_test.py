import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_tex as ET
C, T3, R, TE = ET.C, ET.T3, ET.R, ET.TE
ts = R.THREE['G3 II']()
d = [t for t in ts if not t['tag'].endswith('*')]
s = [t for t in ts if t['tag'].endswith('*')]
c = T3.configs(5)
for a, b in list(zip(d, s))[:4]:
    va = np.conj(a['c'] * TE.term_value(ET.tex(a), T3.env(T3.swapped_w(c))))
    vb = b['c'] * TE.term_value(ET.tex(b), T3.env(c))
    print(va, vb)
from rows_user import rows
OLD = {(r[0].split('\\')[0].strip(), r[1], r[2]): r[4] for r in rows(os.path.join(HERE, '..', 'bfkl_final', 'input', 'bfkl.tex'))}
old = C.prep(OLD[('Group III', 'Row II', None)])
for cfg, lab in ((c, 'c'), (T3.swapped_w(c), 'swapped'), (T3.configs(6), 'seed6')):
    vo = TE.term_value(old, T3.env(cfg))
    ve = sum(t['c'] * TE.term_value(ET.tex(t), T3.env(cfg)) for t in d)
    print(lab, vo, ve)
