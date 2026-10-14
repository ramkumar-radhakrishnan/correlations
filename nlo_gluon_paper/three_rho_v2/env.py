"""Paths: ../two_rho (term format, symmetrization), ../three_rho (engine, sources), ../bfkl (acc3), ../combined (allsets)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
TWO = os.path.join(ROOT, 'two_rho')
BF = os.path.join(ROOT, 'bfkl')
COMB = os.path.join(ROOT, 'combined')
ALL2 = os.path.join(ROOT, 'two_rho_all')
if TWO not in sys.path:
    sys.path.insert(0, TWO)
import core  # noqa: F401,E402
for p in (os.path.join(BF, 'reference'), BF, HERE):
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)
