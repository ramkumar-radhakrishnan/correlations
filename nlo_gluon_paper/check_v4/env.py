"""Paths: ../two_rho (term format, symmetrization, JIMWLK target), ../three_rho (engine, sources),
../two_rho_all (all earlier sets, after the fixes), ../three_rho_v2 (final three-rho equation), ../bfkl (BFKL)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
TWO = os.path.join(ROOT, 'two_rho')
BF = os.path.join(ROOT, 'bfkl')
ALL2 = os.path.join(ROOT, 'two_rho_all')
V2 = os.path.join(ROOT, 'three_rho_v2')
if TWO not in sys.path:
    sys.path.insert(0, TWO)
import core  # noqa: F401,E402  (sets ../three_rho paths)
for p in (os.path.join(BF, 'reference'), BF, ALL2, HERE):
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)
