"""paths: ../finite (closed forms, transcription), and through it ../groups23, ../group1, ../three_rho"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIN = os.path.normpath(os.path.join(HERE, '..', 'finite'))
for sub in ('groups23', 'group1', 'three_rho'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.insert(0, p)
if FIN not in sys.path:
    sys.path.insert(0, FIN)
