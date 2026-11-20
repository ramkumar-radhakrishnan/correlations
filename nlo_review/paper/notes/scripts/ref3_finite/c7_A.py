# Order A of c7_order.py recomputed with the radial grid starting at r = 0 (the integrand K^m f r = e^m f is bounded
# there) and a finer grid; tail beyond t = 30 from a fit H = a/q^2 + b/q^3.
import numpy as np
src = open(__file__.replace('c7_A.py', 'c7_order.py')).read()
src = src.replace("rr = np.linspace(0, 6.5, 3251)[1:]", "rr = np.linspace(0, 6.5, 6501)")
exec(src.split('for m in (0, 1):')[0])
for m in (0, 1):
    ts = np.concatenate([np.linspace(0, 10, 801), np.linspace(10.025, 30, 800)])
    Hs = np.array([H(m, 1 + tt) for tt in ts])
    SA = np.trapezoid(Hs, ts)
    q = 1 + ts[-200:]; Mfit = np.vstack([1 / q ** 2, 1 / q ** 3]).T
    coef = np.linalg.lstsq(Mfit, Hs[-200:], rcond=None)[0]
    SA += coef[0] / (1 + ts[-1]) + coef[1] / (2 * (1 + ts[-1]) ** 2)
    print('m=%d  order A (refined): %.6f%+.6fi   tail part %.2e' % (m, SA.real, SA.imag, abs(coef[0] / (1 + ts[-1]))), flush=True)
