"""Smooth random adjoint Wilson-line field U(x) and colour-charge field rho(x) (SU(3)),
and the colour factors of Group I rows III, V and Group III, typed from the TeX."""
import numpy as np
from scipy.linalg import expm
from lattice import su_n_generators, structure_constants
t = su_n_generators(3); f = structure_constants(t); N = 3; NA = 8
F = f  # f[a,b,c]
rng = np.random.default_rng(123)
_q = rng.normal(size=(6, 2)) * 1.2; _ph = rng.uniform(0, 2*np.pi, (6,)); _c = rng.normal(size=(6, NA)) * 0.9
_qr = rng.normal(size=(5, 2)); _phr = rng.uniform(0, 2*np.pi, (5,)); _cr = rng.normal(size=(5, NA))
def alpha(x): return (_c * np.cos(_q @ x + _ph)[:, None]).sum(0)
def U(x):
    a = alpha(x); M = np.einsum('a,abc->bc', a, F)          # (i a.T)_{bc} = a^a f^{abc}
    return expm(M)                                         # real orthogonal
def rho(x): return (_cr * np.cos(_qr @ x + _phr)[:, None]).sum(0)
def AC(r1, r2): return 2 * np.outer(r1, r2)                # {rho^a, rho^b} for classical charges

def colour(group, P):
    x, y, z, xp, yp, w, wp = (P[k] for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp'))
    Ux, Uy, Uz, Uxp, Uyp, Uw, Uwp = U(x), U(y), U(z), U(xp), U(yp), U(w), U(wp)
    rx, ry, rxp, ryp = rho(x), rho(y), rho(xp), rho(yp)
    if group == 'C69':
        L = (Uxp - Uwp) @ rxp                                                    # [a]
        R = np.einsum('adc,db,ce,eb->a', F, Ux, Uz, AC(ry, ryp)) - np.einsum('abc,cd,be,de->a', F, Uy, Uyp, AC(ry, ryp))
        return L @ R
    if group == 'C610':
        M = -np.einsum('adc,fbe,db,ce->af', F, F, Ux, Uz) + N * Ux
        return np.einsum('ab,b,af,f->', Uxp - Uwp, rxp, M, ry)
    if group == 'C611':
        return np.einsum('ab,b,af,f->', Uxp - Uwp, rxp, Uy - Ux, ry)
    if group == 'R3':      # Group I row III, i f^{abd}[..]_a [..]_{bd}  (i kept in the kernel coefficient)
        L = (Uxp - Uwp) @ rxp
        R = np.einsum('de,bc,e,c->bd', Ux, Uz, rx, ry) - np.einsum('de,bc,e,c->bd', Ux, Uy, rx, ry)
        return np.einsum('abd,a,bd->', F, L, R)
    # Group III, first brackets
    def first(zz, kind):
        Uq = U(zz)
        if kind == 'fUU':   # f^{c'd'a}[U^{c'e'}(q)U^{d'b}(z) - U^{d'b}(z)U^{c'e'}(x')]rho^{e'}(x')  -> [a,b]
            return np.einsum('cda,ce,db,e->ab', F, Uq - Uxp, Uz, rxp)
        if kind == 'U':     # [U^{c'e'}(q) - U^{c'e'}(x')] rho^{e'}(x') -> [c']
            return (Uq - Uxp) @ rxp
    if group in ('X1', 'X2'):
        Lab = first(yp if group == 'X1' else z, 'fUU')
        Rab = np.einsum('ac,cb->ab', Uw, AC(ry, rx)) - np.einsum('cb,ad,ce,de->ab', Uz, Uy, Ux, AC(ry, rx))
        return np.sum(Lab * Rab)
    if group in ('Y1', 'Y2'):
        Lab = first(yp if group == 'Y1' else z, 'fUU')
        Rab = np.einsum('cbd,ac,d->ab', F, Uw, rx) - np.einsum('acd,cb,de,e->ab', F, Uz, Ux, rx)
        return np.sum(Lab * Rab)
    if group in ('W1', 'W2'):
        Lc = first(yp if group == 'W1' else z, 'U')                     # [c']
        Rca = np.einsum('ce,ad,e,d->ca', Ux, Uy, rx, ry) - np.einsum('ad,ce,e,d->ca', Uw, Ux, rx, ry)
        return np.einsum('kca,k,ca->', F, Lc, Rca)                      # f^{c'ca}
    if group in ('Z1', 'Z2'):
        Lab = first(yp if group == 'Z1' else z, 'fUU')
        Rab = np.einsum('cb,ad,ce,d,e->ab', Uz, Uy, Ux, ry, rx) - np.einsum('ad,d,b->ab', Uy, ry, rx)
        return np.sum(Lab * Rab)
    raise ValueError(group)
