"""paths: ../finite_v2 (C_row), ../finite (transcription ft.py), ../bfkl_check (your earlier list, corrected list),
../groups23, ../three_rho, ../two_rho (reordering engine)"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ('two_rho', 'three_rho', 'group1', 'groups23', 'finite', 'finite_v2', 'bfkl_check'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.insert(0, p)
