"""paths to the earlier modules: ../two_rho (term format), ../three_rho (colour numerics), ../final_check (lists, targets)"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FC = os.path.normpath(os.path.join(HERE, '..', 'final_check'))
_spec = importlib.util.spec_from_file_location('fc_env', os.path.join(FC, 'env.py'))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)          # puts ../two_rho, ../three_rho, ../bfkl, ../two_rho_all, ../final_check on sys.path
for p in (FC, HERE):
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)
