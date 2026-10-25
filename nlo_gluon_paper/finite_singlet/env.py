"""paths: ../finite_v2 (decomposition, compact), ../finite, ../groups23, ../group1, ../three_rho, ../two_rho, ../singlet"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ('singlet', 'two_rho', 'finite', 'groups23', 'group1', 'three_rho', 'finite_v2'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.insert(0, p)
