"""Paths: ../two_rho (term format, users, JIMWLK target), ../three_rho (engine, sources), ../bfkl (region-P rows, BFKL)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
TWO = os.path.join(ROOT, 'two_rho')
BF = os.path.join(ROOT, 'bfkl')
THREE = os.path.join(ROOT, 'three_rho')
if TWO not in sys.path:
    sys.path.insert(0, TWO)
import core  # noqa: F401,E402  (sets ../three_rho paths)
for p in (os.path.join(BF, 'reference'), BF, HERE):
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)
