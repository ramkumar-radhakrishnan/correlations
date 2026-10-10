"""Group I Row III, the B3-bar term integrated over k+ + Lambda < p+ < vee (your part 2):
(1) principal-value numerical integration of your integrand vs your integrated result (part 2);
(2) the large-p+ limit of the integrand = -K_h, and the vee-dependence of part 2 has the same slope.
Common factor i g^4 f^{abd}/(8 pi^5 k+) x colour dropped; k+ = 1."""
import warnings
import numpy as np
warnings.filterwarnings('ignore')
from scipy.integrate import quad

rng = np.random.default_rng(3)
P = {l: rng.normal(size=2) for l in ['x', 'y', 'z', "x'", 'w', "w'"]}
x, y, z, xp, w, wp = (P[l] for l in ['x', 'y', 'z', "x'", 'w', "w'"])
X, Y, Zw = xp - wp, y - z, z - w
A, B = x - w, x - z
a, b = A @ A, B @ B
d = lambda u, v: u @ v
kp = 1.0
pre = 1.0 / (d(X, X) * d(Y, Y) * d(Zw, Zw))
Kh = d(X, Zw) * d(Y, B) / b * pre


def term4(p):
    V = (kp * A - p * B) / (kp * a - p * b)
    br = d(X, Zw) * d(Y, V) - kp / p * d(Y, Zw) * d(X, V) - kp / (kp - p) * d(X, Y) * d(Zw, V)
    return br / (kp - p) * pre


def part2(V_, L):
    Kc = d(X, A) * d(Y, Zw) / a * pre
    Kd = d(X, Y) * d(Zw, B) * a / (a - b) ** 2 * pre
    Ke = d(X, Y) * d(Zw, A) * b / (a - b) ** 2 * pre
    Lo = np.log(abs((a * kp - b * V_) / (a * kp - b * (kp + L))))
    blue = np.log(V_ / L) + np.log(1 - kp / V_)
    green = d(X, Y) / (d(X, X) * d(Y, Y) * (a - b)) * (1 / (V_ - kp) - 1 / L) * kp
    g = [d(X, Zw) * d(Y, A) / (a - b), -d(X, A) * d(Y, Zw) * b / ((a - b) * a), d(X, Y) * d(Zw, A) * b / (a - b) ** 2,
         -d(X, Zw) * d(Y, B) * a / ((a - b) * b), d(X, B) * d(Y, Zw) / (a - b), -d(X, Y) * d(Zw, B) * a / (a - b) ** 2]
    return green - np.log(V_ / (L + kp)) * Kc + blue * (Kd - Ke) + Lo * pre * sum(g)


def pv_integral(lo, hi):
    """principal value of int_lo^hi term4 (simple pole at p0 = k+ a/b if inside)"""
    p0 = kp * a / b
    if not (lo < p0 < hi):
        return quad(term4, lo, hi, limit=2000, epsabs=1e-13, epsrel=1e-12)[0]
    h = 0.5 * min(p0 - lo, hi - p0)
    def g(p):
        if abs(p - p0) < 1e-9:
            p = p0 + 1e-7
        return term4(p) * (p - p0)
    mid = quad(g, p0 - h, p0 + h, weight='cauchy', wvar=p0, limit=500, epsabs=1e-13, epsrel=1e-12)[0]
    left = quad(term4, lo, p0 - h, limit=2000, epsabs=1e-13, epsrel=1e-12)[0]
    right = quad(term4, p0 + h, hi, limit=2000, epsabs=1e-13, epsrel=1e-12)[0]
    return left + mid + right


print('pole at p+ = %.4f' % (kp * a / b))
for V_, L in [(10.0, 1e-2), (100.0, 1e-3), (1e3, 1e-3), (1e3, 1e-4)]:
    num = pv_integral(kp + L, V_)
    print('vee=%-8g Lambda=%-8g  PV integral %.10f   your part 2 %.10f   difference %.1e' % (V_, L, num, part2(V_, L), part2(V_, L) - num))
for p in (1e3, 1e4, 1e5):
    print('p+ = %-8g  p+ x integrand = %.8f   (-K_h = %.8f)' % (p, p * term4(p), -Kh))
s = (part2(1e7, 1e-3) - part2(1e6, 1e-3)) / np.log(10)
print('d(part 2)/d log(vee) = %.8f   (-K_h = %.8f): your vee dependence equals the large-p+ limit' % (s, -Kh))
