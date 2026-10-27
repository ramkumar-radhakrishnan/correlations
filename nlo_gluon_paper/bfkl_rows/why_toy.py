r"""Toy model for Why_log_vee.tex: the p+ integrals of the three terms of the bracket of Group III Row VI, line 1, as vee grows.

In the y' integral everything except the singular factor (y'-z)^m/(y'-z)^2 and the phase is replaced by e^{-(y'-z)^2/2}
(|k| = 1).  Then, with t = (p+ + k+)/k+,
    int d^2y' e^{-i t k.(y'-z)} (y'-z)^m/(y'-z)^2 e^{-(y'-z)^2/2} = -2 pi i k^m/k^2 (1 - e^{-t^2/2})/t,
so every F_n(p+) is F^z f(t), f(t) = (1 - e^{-t^2/2})/t, with F^z the y' -> z coefficient (f -> 1/t for t -> infinity).
In units of C F^z/(2 pi), lambda = vee/k+:
    a      =  int_0^{vee-k+} dp+ F(p+)/(p+ + k+)        =  int_1^lambda dt f(t)/t
    b      = -int_0^{vee-k+} dp+ [F(p+) - F(0)]/p+       = -int_1^lambda dt [f(t) - f(1)]/(t - 1)     ([1/p+]_+)
    line2  = -f(1) log(lambda - 1)                                                                 (log(vee/k+ - 1) part)
    c      = -(1/k+) int_0^{vee-k+} dp+ F(p+)            = -int_1^lambda dt f(t)  (exact: -log lambda + [E1(1/2) - E1(lambda^2/2)]/2)
Output: why_toy.json"""
import json
import os
import numpy as np
from scipy.integrate import quad
from scipy.special import exp1

f = lambda t: (1 - np.exp(-t * t / 2)) / t
f1 = f(1.0)


def integral(fun, a, b):
    pts = np.unique(np.concatenate([[a], np.geomspace(a, b, 80), [b]]))
    return sum(quad(fun, pts[i], pts[i + 1], limit=200, epsabs=1e-13, epsrel=1e-12)[0] for i in range(len(pts) - 1))


if __name__ == '__main__':
    out = dict(f1=f1, rows=[])
    for lam in (1e1, 1e2, 1e3, 1e4, 1e6):
        a = integral(lambda t: f(t) / t, 1, lam)
        b = -integral(lambda t: (f(t) - f1) / (t - 1) if t > 1 else 0.0, 1, lam)
        l2 = -f1 * np.log(lam - 1)
        c = -integral(f, 1, lam)
        c_exact = -np.log(lam) + 0.5 * (exp1(0.5) - exp1(lam * lam / 2))
        assert abs(c - c_exact) < 1e-9
        out['rows'].append(dict(lam=lam, a=a, b=b, line2=l2, b_line2=b + l2, c=c, c_plus_log=c + np.log(lam)))
        print('%8.0e  a=%.5f  b=%.5f  line2=%.5f  b+line2=%.5f  c=%.5f  c+log=%.5f' % (lam, a, b, l2, b + l2, c, c + np.log(lam)))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'why_toy.json'), 'w'), indent=1)
