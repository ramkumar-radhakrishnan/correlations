"""Run ../bfkl/reduceP.py (basis reduction, exact coefficients) on a book of this folder:
python3 reduce_run.py bookP2.pkl reducedP2.pkl"""
import env
import os
import runpy
import sys
sys.argv = ['reduceP.py'] + sys.argv[1:]
runpy.run_path(os.path.join(env.BF, 'reduceP.py'), run_name='__main__')
