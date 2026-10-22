"""paths to the earlier check modules"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
for sub in ('groups23', 'group1', 'three_rho'):
    p = os.path.join(ROOT, sub)
    if p not in sys.path:
        sys.path.insert(0, p)
