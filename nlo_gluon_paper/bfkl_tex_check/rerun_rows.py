"""The four rows whose first forms mix indexed and unindexed terms, with the corrected index sums (texfix.py)."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder  # noqa  (patches texeval.term_value)
import cmp_reorder as CR
o = CR.run(only=['G1 IV', 'G2 I first part', 'G3 I', 'G3 II'])
json.dump(o, open(os.path.join(HERE, 'cmp_reorder_fixed.json'), 'w'), indent=1)
