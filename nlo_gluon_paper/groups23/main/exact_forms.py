"""Exact (no Lambda -> 0, no vee -> infinity) forms of the longitudinal integrals of Calculation_main.tex, Sec. 2.2,
each checked against direct quadrature of the integrand as written in the notes (the integrands themselves are
verified in lcpt_vac.py and onerho_dimreg.py).  delta = Lambda/k+.  Results -> exact_forms.json."""
import json
import os
import numpy as np
import mpmath as mp

mp.mp.dps = 30
OUT = {}


def rec(name, exact, quad, notes=None):
    OUT[name] = {'exact': float(exact), 'quadrature': float(quad), 'rel_diff': float(abs(exact - quad) / max(abs(quad), 1e-30))}
    if notes is not None:
        OUT[name]['notes_Lambda0_form'] = float(notes)
    print('%-38s exact %.15g  quad %.15g  notes %s' % (name, exact, quad, 'n/a' if notes is None else '%.10g' % notes))


kp, Lam, vee = mp.mpf('1.7'), mp.mpf('0.11'), mp.mpf('23')
d = Lam / kp
L, l1 = mp.log(d), mp.log((1 - d) / d)

# ---------------------------------------------------------------- gtwoloop1: -2 I0 (pole coefficient) and -2 I1 (finite)
w = lambda x: -2 + x * (1 - x) + 1 / x + 1 / (1 - x)
I0q = mp.quad(w, [d, mp.mpf(1) / 2, 1 - d])
I1q = mp.quad(lambda x: w(x) * (mp.log(x) + mp.log(1 - x)), [d, mp.mpf(1) / 2, 1 - d])
I0e = -mp.mpf(11) / 6 + 4 * d - d**2 + mp.mpf(2) / 3 * d**3 + 2 * l1
M = mp.log(1 - d)
I1e = (mp.mpf(67) / 18 - mp.pi**2 / 3 - mp.mpf(23) / 3 * d + mp.mpf(2) / 3 * d**2 - mp.mpf(4) / 9 * d**3
       + (-mp.mpf(11) / 3 + 4 * d - d**2 + mp.mpf(2) / 3 * d**3) * M + (4 * d - d**2 + mp.mpf(2) / 3 * d**3) * L
       + M**2 - L**2 + 2 * L * M + 4 * mp.polylog(2, d))
rec('gtwoloop1 pole coeff (-2 I0)', -2 * I0e, -2 * I0q, mp.mpf(11) / 3 + 4 * L)
rec('gtwoloop1 finite (-2 I1)', -2 * I1e, -2 * I1q, 2 * L**2 - mp.mpf(67) / 9 + 2 * mp.pi**2 / 3)

# ---------------------------------------------------------------- onerho2: int F_a, F_a = [(2-xi)ln xi + (1+xi)ln(1-xi) + (1+xi) P]/(xi(1-xi))
Fa_pole = mp.quad(lambda x: (1 + x) / (x * (1 - x)), [d, mp.mpf(1) / 2, 1 - d])
Fa_fin = mp.quad(lambda x: ((2 - x) * mp.log(x) + (1 + x) * mp.log(1 - x)) / (x * (1 - x)), [d, mp.mpf(1) / 2, 1 - d])
rec('onerho2 pole coeff', 3 * l1, Fa_pole, -3 * L)
rec('onerho2 finite', 2 * M**2 - 2 * L**2 + 4 * mp.polylog(2, d) - mp.pi**2 / 3 + 2 * L * M, Fa_fin, -mp.pi**2 / 3 - 2 * L**2)

# ---------------------------------------------------------------- instgg: int_Lam^{k+-Lam} dp+ (p++k+)/(p+(k+-p+)^2) * k+/2
Eq = mp.quad(lambda p: (p + kp) / (p * (kp - p)**2), [Lam, kp / 2, kp - Lam]) * kp / 2
Ee = kp / Lam - kp / (kp - Lam) + mp.log((kp - Lam) / Lam)
rec('instgg factor E', Ee, Eq, kp / Lam - 1 - L)

# ---------------------------------------------------------------- 6c + 7: the 1/Lambda terms
q6c = mp.quad(lambda p: 1 / p**2, [Lam, 1, vee])
rec('6c: int dp+/p+^2', 1 / Lam - 1 / vee, q6c, 1 / Lam)
q7 = mp.quad(lambda x: 2 / x**2, [Lam / kp, 1, vee / kp]) / kp
rec('7: (1/k+) int 2 dxi/xi^2 (cancels 6c x2)', 2 * (1 / Lam - 1 / vee), q7, 2 / Lam)

# ---------------------------------------------------------------- Psi^B (4+5onerho) converted to the final momentum k+:
#   int_Lam^vee dK/sqrt(K) int_{Lam/K}^{1-Lam/K} dxi f(xi) |g(xi K)>  =  int_Lam^{vee-Lam} dk/sqrt(k) int_{k/vee}^{k/(k+Lam)} dxi xi^{-1/2} f(xi) |g(k)>
x0, x1 = kp / vee, kp / (kp + Lam)
fb_pole = lambda x: (2 - x) / (x * (1 - x))
fb_fin = lambda x: ((x - 2) * mp.log(x) + (x + 1) * mp.log(1 - x)) / (x * (1 - x))
a, V, W = mp.log(1 + d), mp.log(vee / kp), mp.log(1 - kp / vee)
rec('Psi^B pole coeff', mp.log(vee / (kp + Lam)) + mp.log((vee - kp) / Lam), mp.quad(fb_pole, [x0, x1]),
    V + mp.log((vee - kp) / Lam))
rec('Psi^B finite', -a**2 - L**2 + a * L + V**2 + W**2 + V * W, mp.quad(fb_fin, [x0, x1]), -L**2 + V**2 + W**2 + V * W)
# check of the change of variables with a test function
tf = lambda x: mp.cos(3 * x) * x
lhs = mp.quad(lambda K: mp.quad(lambda x: tf(x) * mp.exp(-x * K), [Lam / K, 1 - Lam / K]) / mp.sqrt(K), [2 * Lam, vee])
rhs = mp.quad(lambda k: mp.exp(-k) / mp.sqrt(k) * mp.quad(lambda x: x**-0.5 * tf(x), [k / vee, k / (k + Lam)]), [Lam, vee - Lam])
rec('change of variables (test function)', rhs, lhs)

# ---------------------------------------------------------------- total one-rho bracket (grhoposition), Lambda < k+ < vee - Lambda
polA = -2 * I0q + Fa_pole - 2 * mp.log(vee / Lam)           # gtwoloop1 + onerho2 + onerhopart, by quadrature
polB = mp.quad(fb_pole, [x0, x1])
pol_exact = mp.mpf(11) / 3 + mp.log(Lam**2 / (kp**2 - Lam**2)) + W - 8 * d + 2 * d**2 - mp.mpf(4) / 3 * d**3
rec('total pole coeff', pol_exact, polA + polB, mp.mpf(11) / 3 + L + mp.log(vee / kp - 1) - mp.log(vee / Lam))
# Psi^A finite = gtwoloop1 finite + onerho2 finite (with the notes' overall sign) + onerhopart finite
finA = -2 * I1q + Fa_fin - (V**2 - L**2)
finB = mp.quad(fb_fin, [x0, x1])
fin_exact = (mp.pi**2 / 3 - mp.mpf(67) / 9 + W**2 + V * W
             + mp.mpf(46) / 3 * d - mp.mpf(4) / 3 * d**2 + mp.mpf(8) / 9 * d**3
             + (mp.mpf(22) / 3 - 8 * d + 2 * d**2 - mp.mpf(4) / 3 * d**3) * M - (8 * d - 2 * d**2 + mp.mpf(4) / 3 * d**3) * L
             - 2 * L * M - 4 * mp.polylog(2, d) - a**2 + a * L)
rec('total finite', fin_exact, finA + finB, mp.pi**2 / 3 - mp.mpf(67) / 9 + W**2 + V * W)


# ---------------------------------------------------------------- two-rho: grho3 + tworho2 (antisymmetrised kernel)
def Vnum(k, p, dd):
    q = k - p
    def B(x):
        D = (1 - x) * (p @ p) + x * (q @ q)
        b = (2 * p - 2 * x * k) * (p @ q) + (2 * (k @ p) - 2 / x * (p @ p)) * q + ((1 + x) / (1 - x) * (q @ q) - k @ k + p @ p) * p
        return b / (x * D)
    return np.array([float(mp.quad(lambda x: B(float(x))[i], [dd, 0.5, 1 - dd])) for i in range(2)])


def Vexact(k, p, dd):
    q = k - p
    P, Q = p @ p, q @ q
    Dl = Q - P
    Dlam = np.log((Q - dd * Dl) / (P + dd * Dl))
    lam = Dlam / Dl
    l = np.log((1 - dd) / dd)
    i2 = 1 / dd - 1 / (1 - dd)
    return q * (2 * (k @ q) * (l - Dlam) / P - 2 * i2) - 2 * (p @ q) * lam * k + 2 * (l + Dlam) * p


rng = np.random.default_rng(7)
dd = float(d)
for t in range(2):
    k, p = rng.normal(size=2), rng.normal(size=2)
    q = k - p
    P, Q = p @ p, q @ q
    ve, vn = Vexact(k, p, dd), Vnum(k, p, dd)
    rec('tworho2 V exact vs quad [%d]' % t, np.linalg.norm(ve), np.linalg.norm(vn))
    OUT['tworho2 V exact vs quad [%d]' % t]['max_component_diff'] = float(np.abs(ve - vn).max())
    # grho3 (exact E) in the same normalisation 1/(32 pi^{9/2} (k-p)^2 k^2): -2 E (k-p)^2/p^2 p
    E = 1 / dd - 1 / (1 - dd) + np.log((1 - dd) / dd)
    S = lambda kk, pp: Vexact(kk, pp, dd) - 2 * E * ((kk - pp) @ (kk - pp)) / (pp @ pp) * pp
    T = S(k, p) / Q - S(k, q) / P                      # antisymmetrised kernel (x2)
    l = np.log((1 - dd) / dd)
    T_closed = 2 * l / (P * Q) * ((p @ q) * (q - p) + P * q - Q * p)
    # the notes' combined first term, with ln(Lambda/k+) -> ln(Lambda/(k+-Lambda)):  L'[k^2 p - (k.p) k]/(p^2 (k-p)^2)
    Lp = np.log(dd / (1 - dd))
    Nn = lambda kk, pp: Lp * ((kk @ kk) * pp - (kk @ pp) * kk) / ((pp @ pp) * ((kk - pp) @ (kk - pp)))
    T_notes = Nn(k, p) - Nn(k, q)
    rec('grho3+tworho2 antisym vs closed [%d]' % t, np.linalg.norm(T), np.linalg.norm(T_closed))
    OUT['grho3+tworho2 antisym vs closed [%d]' % t]['max_component_diff'] = float(np.abs(T - T_closed).max())
    OUT['grho3+tworho2 antisym vs closed [%d]' % t]['vs_notes_with_ln(Lambda/(k+-Lambda))'] = float(np.abs(T - T_notes).max())
    print('    |T - notes(ln Lambda/(k+-Lambda))| =', np.abs(T - T_notes).max())

# ---------------------------------------------------------------- combinedgrhorho: the xi integral of the remainder R (4+5tworhosimpl)
for t in range(2):
    k, p = rng.normal(size=2), rng.normal(size=2)
    q = k - p
    P, Q, K, kpd = p @ p, q @ q, k @ k, k @ p
    n0 = 2 * (P - kpd) * k + K * p
    n1 = (2 * kpd - P) * k
    Rq = np.array([float(mp.quad(lambda x: -(n0[i] + x * n1[i]) / ((x * Q + (1 - x) * K) * (1 - x)), [float(x0), float(x1)])) for i in range(2)])
    al = (P * k + K * p) / Q
    bb = ((2 * P - 2 * kpd + K) * k + K * p) / Q
    kpf, Lf, vf = float(kp), float(Lam), float(vee)
    Re = (-al * np.log((vf - kpf) / Lf) - k * np.log(vf / (kpf + Lf)) + bb * np.log(1 + (vf / kpf - 1) * K / Q)
          - bb * np.log(1 + Lf * K / (kpf * Q)))
    Rn = (-al * np.log((vf - kpf) / Lf) - k * np.log(vf / kpf) + bb * np.log(1 + (vf / kpf - 1) * K / Q))
    rec('R xi-integral exact vs quad [%d]' % t, np.linalg.norm(Re), np.linalg.norm(Rq), np.linalg.norm(Rn))
    OUT['R xi-integral exact vs quad [%d]' % t]['max_component_diff'] = float(np.abs(Re - Rq).max())

json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'exact_forms.json'), 'w'), indent=1)
