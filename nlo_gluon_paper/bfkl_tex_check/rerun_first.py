"""cmp_first.py with the corrected index sums (texfix.py) -> cmp_first.json"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder  # noqa  (patches texeval.term_value)
import runpy
runpy.run_path(os.path.join(HERE, 'cmp_first.py'), run_name='__main__')
