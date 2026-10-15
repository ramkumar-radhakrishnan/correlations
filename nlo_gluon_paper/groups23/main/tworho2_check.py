"""xi integral of the {rho,rho} part of Eq. 3grhorho, compared with Eq. tworho2.
The colour factor f^{abc}{rho^b(-p), rho^c(-k+p)} is odd under p -> k-p, so only the odd part of the kernel counts:
we compare V(p)/(k-p)^2 - V(k-p)/p^2.  Exact: int_{delta}^{1-delta} dxi B(xi)/(xi D(xi)), delta = Lambda/k+.
Notes (tworho2): -(1/p^2)[L k.p k - (L k^2 + 2E (k-p)^2) p],  L = ln delta,  E = 1/delta - 1 - L.
The exact closed form is checked in exact_forms.py."""
import numpy as np
import mpmath as mp

rng = np.random.default_rng(5)


def V_num(k, p, d):
    q = k - p

    def B(x):
        D = (1 - x) * (p @ p) + x * (q @ q)
        b = (2 * p - 2 * x * k) * (p @ q) + (2 * (k @ p) - 2 / x * (p @ p)) * q + ((1 + x) / (1 - x) * (q @ q) - k @ k + p @ p) * p
        return b / (x * D)
    return np.array([float(mp.quad(lambda x: B(float(x))[i], [d, 0.5, 1 - d])) for i in range(2)])


def V_notes(k, p, d):
    q = k - p
    L = np.log(d)
    E = 1 / d - 1 - L
    return -(L * (k @ p) * k - (L * (k @ k) + 2 * E * (q @ q)) * p) / (p @ p)


if __name__ == '__main__':
    for d in [1e-2, 1e-4, 1e-6]:
        k, p = rng.normal(size=2), rng.normal(size=2)
        q = k - p
        # the factor 1/(k-p)^2 outside the xi integral is not symmetric: kernel = V(p)/(k-p)^2
        v = V_num(k, p, d) / (q @ q) - V_num(k, q, d) / (p @ p)
        n = V_notes(k, p, d) / (q @ q) - V_notes(k, q, d) / (p @ p)
        print('delta=%g  exact %s  notes %s  diff %s' % (d, v, n, v - n))
