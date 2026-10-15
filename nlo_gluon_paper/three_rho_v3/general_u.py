"""General (non-symmetric) adjoint U: (a) the terms as written, with each U^{ij} contracted as written;
(b) every Wilson-line chain read as U(p1)...U(pn) rho (U oriented from f towards the charge).
Prints, per class, the largest residual of the sum over terms of the colour factor (kernel fixed within a class)."""
import numpy as np
import env  # noqa
import numerics as nm
from core import NC
import parse3 as P
import networks3 as N

ts, f = nm.structure_constants(NC)
E = P.parse()
for seed in (2, 7):
    rng = np.random.default_rng(seed)
    U = {p: nm.random_adjoint(NC, ts, rng, False) for p in N.PTS}
    r = {p: rng.normal(size=NC * NC - 1) for p in N.PTS}
    for cls in ('alpha', "alpha'", 'beta', "beta'"):
        es = [e for e in E if e['cls'] == cls]
        written = sum((e['pref'] / 1j) * e['sign'] * e['mult'] *
                      nm.eval_colour(P.as_term(e) and ([('U', i, j, p) for i, j, p in e['Us']] + [('f',) + e['fs'][0]] +
                                     [('rho', i, p) for i, p in e['rs']]), f, U, r, NC * NC - 1) for e in es)
        chain = sum((e['pref'] / 1j) * e['sign'] * e['mult'] * N.fingerprint(N.legs_of(e), [(U, r)], f)[0] for e in es)
        scale = max(abs(N.fingerprint(N.legs_of(e), [(U, r)], f)[0]) for e in es)
        print('seed %d %-7s as written: |sum| = %.2e   chains U(p1)..U(pn)rho: |sum| = %.2e   (largest term %.2e)' % (
            seed, cls, abs(written), abs(chain), scale))
