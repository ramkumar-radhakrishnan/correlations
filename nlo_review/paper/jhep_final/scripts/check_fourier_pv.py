#!/usr/bin/env python3
"""Numerical test of the principal-value Fourier transforms (FT_i3), (FT_i4), (FT_i5) of the appendix
"Fourier transforms":

  int d^2k/(2pi) d^2p/(2pi) e^{i k.X + i p.Z} f(k,p) = h(X,Z),

  (i3) f = P 1/(xi k^2 - p^2)                                   h = 1/(xi Z^2 - X^2)
  (i4) f = (k-p)^m (p-xi k)^l / [(k-p)^2 P(p^2 - xi k^2)]       h = -(X+xi Z)^m (X+Z)^l / [(X^2-xi Z^2)(X+Z)^2]
  (i5) f = (k-p)^m (p-xi k)^l / [P(p^2 - xi k^2)(p-xi k)^2]     h = -(X+Z)^l (X+xi Z)^m / [(X^2-xi Z^2)(X+xi Z)^2]

Method.  Both sides are distributions; they are compared after smearing with a Gaussian test function centred at
(X0,Z0):  int d^2X d^2Z h(X,Z) phi(X,Z) = int d^2k d^2p/(2pi)^2 f(k,p) Phi(k,p),
phi = exp(-(X-X0)^2/2a - (Z-Z0)^2/2b),  Phi = (2 pi a)(2 pi b) exp(i k.X0 + i p.Z0 - a k^2/2 - b p^2/2).
With P(1/x) = (1/2i) int_0^infty dt (e^{itx} - e^{-itx}) and 1/q^2 = int_0^infty ds e^{-s q^2} the four-dimensional
transverse integrals are Gaussian and are done in closed form; the remaining integrals over (s,t) are done
numerically.  The two sides are computed independently (momentum space / coordinate space).
Run: python3 -I scripts/check_fourier_pv.py
"""
import numpy as np
from scipy import integrate

def gauss(M, J, real_source, al, be, with_tensor):
    """int d^2u d^2v (al.w)^m (be.w)^l exp(-1/2 w^T M w + c J.w), w=(u,v), c = 1 (real_source) or i.
    M: complex symmetric 2x2 (acts on the pair index), J: 2x2 array (pair index, transverse component).
    Returns the 2x2 tensor in (m,l) (or the scalar integral if with_tensor is False)."""
    Mi = np.linalg.inv(M); det = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
    c = 1.0 if real_source else 1j
    mu = c*(Mi @ J)                                   # mean, shape (2 pair, 2 transverse)
    Z0 = (2*np.pi)**2/det*np.exp(0.5*c*c*np.sum(J*(Mi @ J)))
    if not with_tensor: return Z0
    am = al @ mu; bm = be @ mu
    return Z0*((al @ Mi @ be)*np.eye(2) + np.outer(am, bm))

def lhs(case, xi, a, b, X0, Z0, s, t):
    J = np.array([X0, Z0], float)
    def I(tt):
        M = np.array([[a + 2j*tt*xi, 0], [0, b - 2j*tt]], complex)      # exp(i tt (p^2 - xi k^2))
        if case == 'i3':
            # f = P 1/(xi k^2 - p^2) = - P 1/(p^2 - xi k^2)
            return -gauss(M, J, False, None, None, False)
        al = np.array([1.0, -1.0]); be = np.array([-xi, 1.0])
        q = al if case == 'i4' else be                                  # extra denominator (k-p)^2 or (p-xi k)^2
        M = M + 2*s*np.outer(q, q)
        return gauss(M, J, False, al, be, True)
    return (2*np.pi*a)*(2*np.pi*b)/(2*np.pi)**2*(I(t) - I(-t))/(2j)

def rhs(case, xi, a, b, X0, Z0, s, t):
    J = np.array([X0/a, Z0/b], float); c0 = np.exp(-X0 @ X0/(2*a) - Z0 @ Z0/(2*b))
    def I(tt):
        M = np.array([[1/a - 2j*tt, 0], [0, 1/b + 2j*tt*xi]], complex)  # exp(i tt (X^2 - xi Z^2))
        if case == 'i3':
            return -c0*gauss(M, J, True, None, None, False)
        A = np.array([1.0, xi]); B = np.array([1.0, 1.0])               # (X + xi Z)^m, (X + Z)^l
        q = B if case == 'i4' else A
        M = M + 2*s*np.outer(q, q)
        return -c0*gauss(M, J, True, A, B, True)
    return (I(t) - I(-t))/(2j)

def integrate_st(fun, case):
    """int_0^infty dt [int_0^infty ds] fun(s,t); the half lines are mapped to (0,1)."""
    def comp(sel):
        if case == 'i3':
            g = lambda y: sel(fun(0.0, y/(1 - y)))/(1 - y)**2
            return integrate.quad(g, 0, 1, limit=400, epsabs=1e-11, epsrel=1e-11)[0]
        g = lambda x, y: sel(fun(x/(1 - x), y/(1 - y)))/((1 - x)**2*(1 - y)**2)
        return integrate.dblquad(g, 0, 1, 0, 1, epsabs=1e-9, epsrel=1e-9)[0]
    if case == 'i3':
        return comp(lambda z: z.real) + 1j*comp(lambda z: z.imag)
    out = np.zeros((2, 2), complex)
    for m in range(2):
        for l in range(2):
            out[m, l] = comp(lambda z: z[m, l].real) + 1j*comp(lambda z: z[m, l].imag)
    return out

if __name__ == '__main__':
    worst = 0.0
    configs = [(0.37, 0.8, 1.3, np.array([0.9, -0.4]), np.array([0.3, 1.1])),
               (2.60, 1.1, 0.7, np.array([-0.5, 0.8]), np.array([1.2, 0.2]))]   # xi < 1 and xi > 1
    for case in ('i3', 'i4', 'i5'):
        for (xi, a, b, X0, Z0) in configs:
            L = integrate_st(lambda s, t: lhs(case, xi, a, b, X0, Z0, s, t), case)
            R = integrate_st(lambda s, t: rhs(case, xi, a, b, X0, Z0, s, t), case)
            dev = np.max(np.abs(L - R))/np.max(np.abs(R)); worst = max(worst, dev)
            print('(FT_%s) xi=%.2f: |momentum-space side| = %.6f, |coordinate-space side| = %.6f, rel. dev. %.1e'
                  % (case, xi, np.max(np.abs(L)), np.max(np.abs(R)), dev))
    print('largest relative deviation: %.1e' % worst)
    assert worst < 1e-5
    print('ALL CHECKS PASSED')
