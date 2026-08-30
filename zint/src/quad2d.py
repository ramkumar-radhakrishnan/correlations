"""Direct numerical evaluation of  I[F] = \int d^2 r e^{i q.r} F(r)  in 2D.

Polar coordinates about the origin, q rotated onto the +x axis.
 * angular integral: periodic trapezoid (spectrally accurate)
 * radial integral: adaptive quad on [0,R_n] with R_n = R0 + n*pi/q, followed by
   repeated averaging (Euler transform) of the partial sums -- this cancels the
   slowly decaying oscillatory tail e^{+-i q r}/r^{3/2} to high accuracy.
"""
import numpy as np
from scipy.integrate import quad

def integrate(F, q, R0=None, nlev=14, ntheta=512, split=()):
    """F(rvec) -> complex, rvec shape (...,2).  q: 2-vector."""
    qn = np.hypot(*q)
    c, sn = q[0]/qn, q[1]/qn                      # rotate q -> (qn,0)
    Rot = np.array([[c, -sn], [sn, c]])           # maps (x,y) -> lab frame
    th = 2*np.pi*np.arange(ntheta)/ntheta
    ph = np.exp(1j*th)

    def g(r):                                     # angular integral at radius r
        v = np.stack([r*np.cos(th), r*np.sin(th)], axis=-1) @ Rot.T
        return np.trapezoid(np.exp(1j*qn*r*np.cos(th))*F(v), th, axis=0) \
               if False else (2*np.pi/ntheta)*np.sum(np.exp(1j*qn*r*np.cos(th))*F(v))

    def radial(a, b):
        re = quad(lambda r: (r*g(r)).real, a, b, limit=400, epsabs=1e-12, epsrel=1e-11)[0]
        im = quad(lambda r: (r*g(r)).imag, a, b, limit=400, epsabs=1e-12, epsrel=1e-11)[0]
        return re + 1j*im

    if R0 is None:
        R0 = 6.0/qn
    edges = [0.0] + sorted([x for x in split if 0 < x < R0]) + [R0]
    S = [sum(radial(edges[i], edges[i+1]) for i in range(len(edges)-1))]
    step = np.pi/qn
    for n in range(nlev):
        S.append(S[-1] + radial(R0 + n*step, R0 + (n+1)*step))
    S = np.array(S)
    for _ in range(nlev):                          # repeated averaging
        S = 0.5*(S[:-1] + S[1:])
    return S[0]
