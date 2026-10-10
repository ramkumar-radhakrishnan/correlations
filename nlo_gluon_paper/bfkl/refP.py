"""Region-P reference rows and BFKL target, in units of N_P (log(vee/k+) in place of log(vee/Lambda))."""
import env  # noqa
from core import term, symmetrize
import rowsT_sym as RT
import rowsP_sym as RP
import bsym

PMAP = {'wb': "w'", 'xb': "x'"}


def conv(t, tag, factor=2.0, pmap=None):
    pm = pmap or {'wb': "w'"}
    t = RT.fix_kern(t)
    g = lambda p: pm.get(p, p)
    kern = [('KK', g(a), g(b), g(c), g(d)) for (a, b), (c, d) in t.kern]
    return term(factor * t.c, [(i, j, g(p)) for i, j, p in t.Us], t.fs, [(i, g(p)) for i, p in t.rs], kern, tag)


def rows(keep=(2,)):
    R = {**RP.virt_rows(), **RP.real_rows()}
    out = {}
    for name, L in R.items():
        ts = []
        for t in L:
            ts += symmetrize(conv(t, 'REFP:' + name), keep)
        out[name] = ts
    return out


def bfkl(keep=(2,), factor=2.0):
    ts = []
    for t in bsym.bfkl_target():
        ts += symmetrize(conv(t, 'BFKL', factor, PMAP), keep)
    return ts
