"""Group III Row II: the engine's three-rho operator (../bfkl_final/rows.py) against the old LaTeX row
(../bfkl_final/input/bfkl.tex), classically (symmetrized over the charges; c.c. included on both sides)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_tex as ET
C, T3, R = ET.C, ET.T3, ET.R
from rows_user import rows
OLD = {(r[0].split('\\')[0].strip(), r[1], r[2]): r[4] for r in rows(os.path.join(HERE, '..', 'bfkl_final', 'input', 'bfkl.tex'))}
old = C.prep(OLD[('Group III', 'Row II', None)])
ts = R.THREE['G3 II']()
c = T3.configs(5)
e = sum(t['c'] * 3 ** t.get('ncp', 0) * T3.sym_value(ET.tex(t), c, False) for t in ts)
o = T3.sym_value(old, c, True)
print('engine', e, 'old tex', o, 'ratio', e / o)
