"""Paths: ../two_rho (term format, users, targets), ../three_rho (engine, sources), reference/ (region-P rows)."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
TWO = os.path.normpath(os.path.join(HERE, '..', 'two_rho'))
for p in (TWO,):
    if p not in sys.path:
        sys.path.insert(0, p)
import core  # noqa: F401,E402  (sets ../three_rho paths)
for p in (os.path.join(HERE, 'reference'), HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
