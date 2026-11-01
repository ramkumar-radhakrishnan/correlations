"""cmp_reorder on every three-rho row, Group I Row II with its correct first form (eq:B:G1II2a of Checks_JIMWLK_BFKL)."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmp_reorder as CR
import cmp_first as C
o = CR.run(first_override={'G1 II': [C.ours('G1II2a')]})
json.dump(o, open(os.path.join(HERE, 'cmp_reorder.json'), 'w'), indent=1)
