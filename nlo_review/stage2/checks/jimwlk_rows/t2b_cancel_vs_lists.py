# do the "=0" blocks contain exactly the list terms (same operator orderings)?  ordered comparison, correct kernel substituted for beta, beta'
import sys, io, contextlib; sys.path.insert(0,".")
with contextlib.redirect_stdout(io.StringIO()):
    from t2_cancel3 import *
from t4_brackets import odiff
KER = {"alpha": ("x'-w'", "x-w", "x-z", "y-z"), "alpha'": ("x'-w'", "x-w", "x'-z", "y-z"), "beta": ("x'-w'", "x-w", "w-z", "y-z"), "beta'": ("x'-w'", "x-w", "w'-z", "y-z")}
order = ["alpha"] * 3 + ["alpha'"] * 3 + ["beta"] * 3 + ["beta'"] * 3
W = World(81)
for (l, hk, ms), tg in zip(allc, order):
    nU = sorted(set(len(m.Us) for m in ms)); assert len(nU) == 1
    for m in ms: m.dots = [(nc_vec(KER[tg][0]), nc_vec(KER[tg][1])), (nc_vec(KER[tg][2]), nc_vec(KER[tg][3]))] if False else m.dots
    import nc
    for m in ms: m.dots = [(nc._vec(KER[tg][0]), nc._vec(KER[tg][1])), (nc._vec(KER[tg][2]), nc._vec(KER[tg][3]))]
    sel = [m for m in allm if m.tag == tg and len(m.Us) == nU[0]]
    d, sz = odiff(ms, sel, W)
    print("block line %d (%s, %dU): header kernel typed %s ; ordered max|block - list class| = %.2e (size %.2g)" % (l, tg, nU[0], "OK" if tuple(hk) == KER[tg] else "K(%s).K(%s)K(%s).K(%s) [expected K(%s).K(%s)]" % (*hk, KER[tg][2], KER[tg][3]), d, sz))
