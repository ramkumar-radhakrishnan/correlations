# s16d_xi.py -- the longitudinal integrals of Sec. 8 written with xi (author's definitions of 2026-11-23):
#   xi = p+/k+ in Groups I and III (Group I: 0<xi<1, Group III: 0<xi<inf), xi = p+/(p+ + k+) in Group II Row IV.
# For every row the p+-dependent part of the kernel (random complex transverse coefficients A, B, C, ...; random
# k+, recoil argument s = k.(...), and the energy-denominator coefficients a, b > 0) is integrated by quadrature in
# the p+ form of Sec. 6 and in the xi form of Sec. 8; the plus prescriptions are implemented by their definitions
# (subtraction of the whole p+-dependent function at the subtraction point).  Group III (upper limit infinity) with a
# convergence factor e^{-eps p+/k+} = e^{-eps xi} (the change of variables is exact for every eps; eps -> 0 checked
# against the closed forms).  Run: python3 -I s16d_xi.py -> s16d_xi.log
import os, numpy as np
from scipy import integrate, special
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
def cquad(f, a, b, **kw):
    kw.setdefault('limit', 2000); kw.setdefault('epsabs', 1e-13); kw.setdefault('epsrel', 1e-12)
    re_ = integrate.quad(lambda t: f(t).real, a, b, **kw)[0]; im_ = integrate.quad(lambda t: f(t).imag, a, b, **kw)[0]
    return re_ + 1j * im_
rng = np.random.default_rng(2026)
def rc(): return complex(rng.normal(), rng.normal())
res = {}
def rec(name, a, b):
    d = abs(a - b) / max(1, abs(a)); res[name] = max(res.get(name, 0), d)
for trial in range(6):
    kp = rng.uniform(0.4, 3.0); s = rng.uniform(-6, 6); A, B, C, A1, A2, A3, A4 = [rc() for _ in range(7)]
    a_, b_ = rng.uniform(0.1, 3, 2)
    # (a) Group I Row III: int_0^{k+} dp+ [1/p+]_+ e^{-is p+/k+}  ==  int_0^1 dxi [1/xi]_+ e^{-is xi}
    Jp = cquad(lambda p: (np.exp(-1j * s * p / kp) - 1) / p, 0, kp)
    Jx = cquad(lambda x: (np.exp(-1j * s * x) - 1) / x, 0, 1)
    rec('I.III  int_0^{k+}dp+[1/p+]_+e^{-ik.(z-w)p+/k+} = int_0^1 dxi[1/xi]_+e^{-ik.(z-w)xi}', Jp, Jx)
    rec('I.III  closed form Ci|s|-gamma_E-log|s|-i Si(s)', Jx, special.sici(abs(s))[1] - np.euler_gamma - np.log(abs(s)) - 1j * special.sici(s)[0])
    # (b) Group I Row V, three-rho kernel (prefactor 1/k+^2 int dp+/(2pi) -> 1/k+ int dxi/(2pi))
    def Vp(p):
        f = lambda q: np.exp(-1j * s * q / kp)
        g = lambda q: np.exp(-1j * s) * np.exp(1j * s * q / kp)
        return (f(p) * A - kp / p * (f(p) - f(0)) * B - kp / p * (g(p) - g(0)) * C) / (2 * np.pi)
    def Vx(x):
        f = lambda q: np.exp(-1j * s * q); g = lambda q: np.exp(-1j * s) * np.exp(1j * s * q)
        return (f(x) * A - (f(x) - f(0)) / x * B - (g(x) - g(0)) / x * C) / (2 * np.pi)
    rec('I.V (3rho)  (1/k+^2) int_0^{k+} dp+/(2pi) {..} = (1/k+) int_0^1 dxi/(2pi) {..}', cquad(Vp, 0, kp) / kp ** 2, cquad(Vx, 0, 1) / kp)
    # (c) Group I Row V, genuine two-rho kernel
    def V2p(p):
        f = lambda q: np.exp(-1j * s * q / kp)
        reg = kp * f(p) / (p * a_ + (kp - p) * b_) * (p * (kp - p) / kp ** 2 * A1 + (kp - p) / kp * A2 + p / kp * A3 + A4)
        return (reg - kp / p * (f(p) - f(0)) * B - kp / (kp - p) * (f(p) - f(kp)) * C) / (2 * np.pi)
    def V2x(x):
        f = lambda q: np.exp(-1j * s * q)
        reg = f(x) / (x * a_ + (1 - x) * b_) * (x * (1 - x) * A1 + (1 - x) * A2 + x * A3 + A4)
        return (reg - (f(x) - f(0)) / x * B - (f(x) - f(1)) / (1 - x) * C) / (2 * np.pi)
    rec('I.V (2rho)  (1/k+^2) int_0^{k+} dp+/(2pi) {..} = (1/k+) int_0^1 dxi/(2pi) {..}', cquad(V2p, 0, kp) / kp ** 2, cquad(V2x, 0, 1) / kp)
    # (d) Group II Row IV, xi = p+/(p+ + k+), finite vee
    for vr in (4.0, 12.0):
        vee = vr * kp
        Ip = cquad(lambda p: p / (p + kp) ** 2 * np.exp(-1j * s) * np.exp(-1j * s * p / kp) / (2 * np.pi), 0, vee - kp) / kp
        Ix = cquad(lambda x: x / (1 - x) * np.exp(-1j * s / (1 - x)) / (2 * np.pi), 0, 1 - kp / vee) / kp
        rec('II.IV  int_0^{vee-k+} dp+ p+/(p++k+)^2 e^{-ik.(y-y\')(1+p+/k+)} = int_0^{1-k+/vee} dxi xi/(1-xi) e^{-ik.(y-y\')/(1-xi)}', Ip, Ix)
        # the same with the DGLAP fraction xi_D = k+/(p+ + k+) = 1 - xi
        Id = cquad(lambda u: (1 - u) / u * np.exp(-1j * s / u) / (2 * np.pi), kp / vee, 1) / kp
        rec('II.IV  same with xi_D = 1 - xi (DGLAP fraction): int_{k+/vee}^1 dxi_D (1-xi_D)/xi_D e^{-ik.(y-y\')/xi_D}', Ip, Id)
    # (e) Group III Rows III (first part), V, VI: (1/k+) int_0^inf dp+/(2pi) e^{-is p+/k+} {A/(p++k+) - [k+/(p+(p++k+))]_+ B - C/k+}
    for eps in (0.5, 0.05):
        def G3p(p):
            f = lambda q: np.exp(-(1j * s + eps) * q / kp)
            return (f(p) * A / (p + kp) - kp / (p * (p + kp)) * (f(p) - f(0)) * B - f(p) * C / kp) / (2 * np.pi)
        def G3x(x):
            f = lambda q: np.exp(-(1j * s + eps) * q)
            return (f(x) * A / (1 + x) - 1 / (x * (1 + x)) * (f(x) - f(0)) * B - f(x) * C) / (2 * np.pi)
        lim = 60 / eps
        rec('III.III/V/VI  (1/k+) int_0^inf dp+/(2pi){..} = (1/k+) int_0^inf dxi/(2pi){..}  [conv. factor e^{-eps xi}]',
            cquad(G3p, 0, lim * kp) / kp, cquad(G3x, 0, lim) / kp)
        # (f) Group III Row III, genuine two-rho kernel: pieces 1-3
        def H1p(p):
            f = np.exp(-(1j * s + eps) * p / kp); den = 1 / (p * a_ + kp * b_)
            return f * den * (p * A1 / (p + kp) ** 2 + p / (kp * (p + kp)) * A2 + A3 / (p + kp) + A4 / (p + kp) - B / kp) / (2 * np.pi)
        def H1x(x):
            f = np.exp(-(1j * s + eps) * x); den = 1 / (x * a_ + b_)
            return f * den * (x * A1 / (1 + x) ** 2 + x / (1 + x) * A2 + A3 / (1 + x) + A4 / (1 + x) - B) / (2 * np.pi)
        rec('III.III (2rho) piece 1:  int_0^inf dp+/(2pi) [..] = (1/k+) int_0^inf dxi/(2pi) [..]', cquad(H1p, 0, lim * kp), cquad(H1x, 0, lim) / kp)
        def H2p(p):
            f = lambda q: np.exp(-(1j * s + eps) * q / kp)
            return kp / (p * (p + kp)) * (f(p) - f(0)) * C / (2 * np.pi)
        def H2x(x):
            f = lambda q: np.exp(-(1j * s + eps) * q)
            return 1 / (x * (1 + x)) * (f(x) - f(0)) * C / (2 * np.pi)
        rec('III.III (2rho) piece 2:  (1/k+) int dp+/(2pi)[k+/(p+(p++k+))]_+ = (1/k+) int dxi/(2pi)[1/(xi(1+xi))]_+', cquad(H2p, 0, lim * kp) / kp, cquad(H2x, 0, lim) / kp)
        def H3p(p): return np.exp(-(1j * s + eps) * p / kp) * p / (p * a_ + kp * b_) * C / (2 * np.pi)
        def H3x(x): return np.exp(-(1j * s + eps) * x) * x / (x * a_ + b_) * C / (2 * np.pi)
        rec('III.III (2rho) piece 3:  (1/k+^2) int dp+/(2pi) p+/(p+a+k+b) = (1/k+) int dxi/(2pi) xi/(xi a+b)', cquad(H3p, 0, lim * kp) / kp ** 2, cquad(H3x, 0, lim) / kp)
    # eps = 0 against the closed forms (Fourier integrals on [0, inf) by QAWF; the distributions of the conventions)
    def four(g, lo):   # int_lo^inf g(xi) e^{-is xi} dxi
        w = abs(s); sg = np.sign(s)
        c = integrate.quad(lambda t: g(t + lo), 0, np.inf, weight='cos', wvar=w, limlst=200)[0]
        si = integrate.quad(lambda t: g(t + lo), 0, np.inf, weight='sin', wvar=w, limlst=200)[0]
        return np.exp(-1j * s * lo) * (c - 1j * sg * si)
    Ix = four(lambda t: 1 / (1 + t), 0)
    rec('III  int_0^inf dxi e^{-is xi}/(1+xi) = e^{is}E_1(is)  (eps = 0, QAWF)', Ix, np.exp(1j * s) * special.exp1(1j * s))
    Ip = cquad(lambda x: (np.exp(-1j * s * x) - 1) / (x * (1 + x)), 0, 1) + four(lambda t: 1 / (t * (1 + t)), 1) - np.log(2)
    rec('III  int_0^inf dxi [1/(xi(1+xi))]_+ e^{-is xi} = -gamma_E - log|s| - i pi/2 sgn(s) - e^{is}E_1(is)  (eps = 0)', Ip,
        -np.euler_gamma - np.log(abs(s)) - 0.5j * np.pi * np.sign(s) - np.exp(1j * s) * special.exp1(1j * s))
for k, v in res.items(): pr('%-150s max rel dev %.1e' % (k, v))
# control: a wrong Jacobian (dp+ -> dxi without the factor k+) must be detected
kp, s = 1.9, 2.3
Vp_ = cquad(lambda p: np.exp(-1j * s * p / kp), 0, kp) / kp ** 2; Vx_ = cquad(lambda x: np.exp(-1j * s * x), 0, 1) / kp ** 2
pr('control: I.V with the prefactor 1/k+^2 kept after dp+ -> dxi (missing Jacobian k+): rel dev %.2f' % (abs(Vp_ - Vx_) / abs(Vp_)))
open(os.path.join(HERE, 's16d_xi.log'), 'w').write('\n'.join(out) + '\n')
