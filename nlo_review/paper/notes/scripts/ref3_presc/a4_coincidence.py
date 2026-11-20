# Ref.3 / Major 5(a).  The only non-integrable point of the on-shell kernel of Row I.III is z -> w (integrated gluon at
# the observed one), where Brace^{ij} -> delta^{ij}/(2 (z-w)^2).  It is regulated only by the window |p+ - k+| > Lambda of
# B_3, i.e. |A - B| > lambda B (lambda = Lambda/k+, A=(x-w)^2, B=(x-z)^2), which near w is a band of half-width lambda|x-w|/2
# around the circle B = A.  Claim:  int_{|z-w|<rho0} d^2z  theta_window / (z-w)^2 = 2 pi ln(rho0/(lambda |x-w|)) + O(lambda, rho0/|x-w|).
# The ln(1/lambda) multiplies  A_LO^dagger (f.Qbar) A_LO  and cancels against the complex conjugate (forward Coulomb phase).
import numpy as np
from scipy.integrate import quad

def I_window(lam, rho0, xw=1.0):
    # polar coordinates around w, x at distance xw along the 1-axis: a = x-w = (xw,0), b = x-z = a - c
    a = np.array([xw, 0.0])
    def radial(phi):
        e = np.array([np.cos(phi), np.sin(phi)])
        # A - B = 2 a.c - c^2 ;  B = A - 2a.c + c^2 ;  allowed: |2 a.c - c^2| > lam (A - 2a.c + c^2)
        rs = np.geomspace(1e-12, rho0, 200001)
        ac = rs * (a @ e)
        ok = (np.abs(2 * ac - rs * rs) > lam * (xw * xw - 2 * ac + rs * rs)).astype(float)
        lr = np.log(rs)
        return np.sum(0.5 * (ok[1:] + ok[:-1]) * np.diff(lr))
    # below r ~ lam*xw everything is excluded, so the lower grid end is irrelevant
    val = quad(radial, 0, 2 * np.pi, limit=400, points=[np.pi / 2, 3 * np.pi / 2])[0]
    return val

for lam in (1e-2, 1e-3, 1e-4):
    for rho0 in (0.05, 0.01):
        num = I_window(lam, rho0)
        print(f"lambda={lam:.0e} rho0={rho0}:  numerical {num: .5f}   2pi ln(rho0/lambda) = {2*np.pi*np.log(rho0/lam): .5f}   diff {num-2*np.pi*np.log(rho0/lam): .2e}")
