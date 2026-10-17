"""principal-value p+ integrals of the Row III / Row IV integrands"""
import numpy as np
from scipy.integrate import quad


def pv(g, a, b, pstar, A, n=200):
    """PV int_a^b g(p)/(A (p - p*)) dp for smooth complex g; returns complex"""
    def part(fun):
        if a < pstar < b:
            r1 = quad(fun, a, b, weight='cauchy', wvar=pstar, limit=n, epsabs=0, epsrel=1e-12)[0]
            return r1 / A
        return quad(lambda p: fun(p) / (p - pstar), a, b, limit=n, epsabs=0, epsrel=1e-12,
                    points=None)[0] / A
    return part(lambda p: g(p).real) + 1j * part(lambda p: g(p).imag)


def plain(f, a, b, n=200):
    return (quad(lambda p: f(p).real, a, b, limit=n, epsabs=0, epsrel=1e-12)[0]
            + 1j * quad(lambda p: f(p).imag, a, b, limit=n, epsabs=0, epsrel=1e-12)[0])
