"""paths: ../finite_v2 (C_row = exact log(vee/k+) coefficients), ../finite (transcription ft.py), ../groups23, ../three_rho"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ('three_rho', 'group1', 'groups23', 'finite', 'finite_v2'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.insert(0, p)
