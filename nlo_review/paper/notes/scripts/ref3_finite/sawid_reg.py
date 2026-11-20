# Eq. (6.15) (eq:sawid): the two Fourier integrals over r = x - w, each "conditionally convergent", evaluated with
# different regulators: isotropic Abel e^{-eps r}, two anisotropic Abel regulators e^{-eps c(alpha) r} (non-circular
# cut-off regions) and a Gaussian-in-r regulator.  Lemma (Sec. 4 of the note): the limit eps -> 0 is regulator
# independent when the Fourier transform is continuous at q = k, which holds for cos(theta) != 0.
# Radial integrals analytic, angular integral by adaptive quadrature.  k = (1,0), z-direction theta.
import numpy as np
from scipy.integrate import quad
from scipy.special import erfcx
GE = np.euler_gamma
th = 0.7
def chi(t): return np.log(abs(2 * np.cos(t))) + 0.5j * np.pi * np.sign(np.cos(t))
def c_theta(t): return np.log(abs(np.cos(t))) + 0.5j * np.pi * np.sign(np.cos(t))
def phi_of(al):   # angle from r = x - w to z, in (-pi, pi]
    return (th - al + np.pi) % (2 * np.pi) - np.pi
# radial integrals  R0(s) = int_0^inf dr e^{-s r},  R1(s) = int_0^inf dr ln r e^{-s r},  Re s > 0
def R0(s): return 1 / s
def R1(s): return -(GE + np.log(s)) / s
regs = {
    'Abel isotropic': lambda al: 1.0,
    'Abel ellipse-like c=1+0.6cos(2al+0.3)': lambda al: 1 + 0.6 * np.cos(2 * al + 0.3),
    'Abel egg c=1+0.5cos(al-1.1)': lambda al: 1 + 0.5 * np.cos(al - 1.1),
}
def integrals(eps, creg):
    # I1^i = int d^2r e^{ik.r} K^i(r) [gE + ln(|k| r) + c(theta)],  K^i = rhat^i / r  ->  int dal rhat^i int dr (...)
    # I2^i = int d^2r e^{ik.r} phi e_perp^i / r ,  e_perp = rhat(x-w) rotated by +pi/2, r = x - w
    out1 = np.zeros(2, complex); out2 = np.zeros(2, complex)
    for i in range(2):
        def f1(al, part):
            s = eps * creg(al) - 1j * np.cos(al)          # e^{i k.r} = e^{i r cos al}
            v = (GE + c_theta(th)) * R0(s) + R1(s)
            e = np.array([np.cos(al), np.sin(al)])[i]
            return (v * e).real if part == 0 else (v * e).imag
        def f2(al, part):
            s = eps * creg(al) - 1j * np.cos(al)
            ep = np.array([-np.sin(al), np.cos(al)])[i]
            v = phi_of(al) * ep * R0(s)
            return v.real if part == 0 else v.imag
        pts = [np.pi / 2, 3 * np.pi / 2, (th + np.pi) % (2 * np.pi)]
        for f, o in ((f1, 'o1'), (f2, 'o2')):
            val = sum((1 if p == 0 else 1j) * quad(f, 0, 2 * np.pi, args=(p,), points=pts, limit=4000, epsabs=1e-11, epsrel=1e-11)[0] for p in (0, 1))
            if o == 'o1': out1[i] = val
            else: out2[i] = val
    return out1, out2
pred = 2j * np.pi * np.array([1.0, 0.0]) * chi(th)
print('prediction  I1 = +2 pi i khat chi(theta) = %s ,  I2 = -I1   (theta = %.2f)' % (np.round(pred, 6), th))
for name, creg in regs.items():
    for eps in (1e-2, 1e-3, 1e-4):
        I1, I2 = integrals(eps, creg)
        print('%-40s eps=%.0e  I1=(%.5f%+.5fi, %.5f%+.5fi)  I2=(%.5f%+.5fi, %.5f%+.5fi)  |I1+I2|=%.2e  |I1-pred|=%.2e' % (
            name, eps, I1[0].real, I1[0].imag, I1[1].real, I1[1].imag, I2[0].real, I2[0].imag, I2[1].real, I2[1].imag,
            np.abs(I1 + I2).max(), np.abs(I1 - pred).max()), flush=True)
