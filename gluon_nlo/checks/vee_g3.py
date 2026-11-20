"""Group III: large-p+ (vee) behaviour.  Compare the contour-rotated value (= distributional vee->inf limit)
with direct real-axis integration up to Xi = (vee-k+)/k+, and extract the non-decaying coefficient c of
exp(-i A xi) at large xi (A = k.(z-y'))."""
import numpy as np
from scipy import integrate
import pnum
rng = np.random.default_rng(1)
P = {k: rng.normal(size=2) for k in ('x', 'y', 'z', 'xp', 'yp', 'w', 'wp')}
A = pnum.KV @ (P['z'] - P['yp'])
def real_axis(lab, Xi):
    g, c, meas, var, f = pnum.FUN[lab]
    mfac = {'1/K^2': 1, '1/K': 1, '1': 1}[meas]
    F0 = lambda p: f(**dict(P, p=p, K=1.0, **{n: 0.0 for n in pnum.PLNAMES}))
    FG = lambda p: f(**dict(P, p=p, K=1.0, **{n: (1.0 if n == 'PLG' else 0.0) for n in pnum.PLNAMES})) - F0(p)
    g0 = FG(0.0)
    integrand = lambda p: F0(p) + (FG(p) - g0) / (p * (p + 1)) if p > 0 else F0(p)
    re_ = integrate.quad(lambda p: np.real(integrand(p)), 0, Xi, limit=4000)[0]
    im_ = integrate.quad(lambda p: np.imag(integrand(p)), 0, Xi, limit=4000)[0]
    return c * mfac * (re_ + 1j * im_) / (2 * np.pi)
print(f'A = k.(z-y\') = {A:+.3f}')
for lab in ['8.3c', '8.3a', '8.2a', '8.2b1', '8.3b1']:
    rot = pnum.term_value(lab, P)
    vals = [real_axis(lab, Xi) for Xi in (50, 50 + np.pi / abs(A), 400, 400 + np.pi / abs(A))]
    # large-xi coefficient c: integrand * exp(+iA xi) at very large xi
    g, c, meas, var, f = pnum.FUN[lab]
    big = [f(**dict(P, p=X, K=1.0, **{n: 0.0 for n in pnum.PLNAMES})) * np.exp(1j * A * X) for X in (1e4, 1e6)]
    print(f'{lab:6s} rotated {rot:+.5f} | real-axis Xi=50,50+pi/A,400,400+pi/A: ' +
          ' '.join(f'{v:+.5f}' for v in vals) + f' | large-xi coefficient of e^(-iA xi): {big[1]:+.3e}')
