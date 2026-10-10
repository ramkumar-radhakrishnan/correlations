"""JIMWLK x LO target: -N Htilde O_LO (operator level, units of N)."""
from core import term
import hsym  # reference/hsym.py

PMAP = {'wb': "w'", 'xb': "x'"}


def conv(t, tag):
    g = lambda p: PMAP.get(p, p)
    kern = [('KK', g(a), g(b), g(c), g(d)) for (a, b), (c, d) in t.kern]
    return term(t.c, [(i, j, g(p)) for i, j, p in t.Us], t.fs, [(i, g(p)) for i, p in t.rs], kern, tag)


def jimwlk(factor=-1.0):
    return [conv(t, 'JIMWLK') for t in hsym.H_on_LO(factor)]


def lo():
    return [conv(t, 'LO') for t in hsym.LO_terms()]
