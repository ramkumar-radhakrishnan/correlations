"""Region-T reference rows (rows/rowsT_sym.py), converted and symmetrized to two-rho terms."""
from core import term, symmetrize
import rowsT_sym as RT  # ../three_rho/reference/rowsT_sym.py

PMAP = {'wb': "w'"}


def conv(t, tag, factor=2.0):
    """reference rows are normalized as 2N x [Row]"""
    t = RT.fix_kern(t)
    g = lambda p: PMAP.get(p, p)
    kern = [('KK', g(a), g(b), g(c), g(d)) for (a, b), (c, d) in t.kern]
    return term(factor * t.c, [(i, j, g(p)) for i, j, p in t.Us], t.fs, [(i, g(p)) for i, p in t.rs], kern, tag)


def rows(keep=(2,)):
    R = {**RT.virt_rows(), **RT.real_rows()}
    out = {}
    for name, L in R.items():
        ts = []
        for t in L:
            ts += symmetrize(conv(t, 'REF:' + name), keep)
        out[name] = ts
    return out
