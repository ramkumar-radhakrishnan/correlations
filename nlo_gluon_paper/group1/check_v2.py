"""Checks of the revised Group_I.tex (input/Group_I_v2.tex) -> check_v2.json.

 1. Row II part 2: the new bracket (with -k^i log(vee/k+)) equals A^(3)_2 of Calculation_main.tex
 2. Row III: the new Eq. (1.10) (fourth line with a minus sign) against first principles with E3
 3. Row III: the new N_2 (overall minus) against the principal-value p+ integral of the new fourth line
 4. Row III: the green terms (1/Lambda, 1/(vee-k+), 1/(k+-Lambda)) now cancel
 5. Row V 2b finite: the 'Simplifying' equation and the final plus-prescribed equation, as written, against the exact p+ integral
 6. running coupling: the beta_0 terms of Row II part 1 and Row V 2b (UV) cancel per ordering (random U, rho; sympy, exact in Lambda)
"""
import json
import os

import numpy as np
import sympy as sp

import fp1
import user1
from fp1 import Cfg, PI, G
from integ import pv

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
out = {}
rng = np.random.default_rng(11)


def rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


# 1. Row II part 2 bracket
dev = []
for t in range(8):
    kv, pv_ = rng.normal(size=2), rng.normal(size=2)
    kp, Lam = rng.uniform(0.5, 2), rng.uniform(1e-3, 0.4)
    vee = kp + rng.uniform(0.1, 50)
    qv = kv - pv_
    K2, P2, Q2, KP = kv @ kv, pv_ @ pv_, qv @ qv, kv @ pv_
    E = (K2 * pv_ - KP * kv) / Q2
    C = ((2 * P2 - 2 * KP + K2) * kv + K2 * pv_) / Q2
    D = (P2 * kv + K2 * pv_) / Q2
    LL, LK = np.log(vee / Lam), np.log(vee / kp)
    main = (np.log(Lam / kp) * E - 2 * kv * LL - 2 * kv * np.log(kp * P2 / (vee * K2 + kp * P2))
            + C * np.log(1 + (vee / kp - 1) * K2 / Q2) - D * np.log((vee - kp) / Lam) - kv * LK)
    new = ((-LL + LK) * E - 2 * kv * LL - 2 * kv * np.log(kp * P2 / (vee * K2 + kp * P2))
           + C * np.log(1 + (vee / kp - 1) * K2 / Q2) - D * (LL + np.log(1 - kp / vee)) - kv * LK)
    dev.append(float(np.abs(new - main).max() / np.abs(main).max()))
out['RowII_part2_new_bracket_vs_main'] = max(dev)

# 2. Row III integrand with the new sign of the fourth line
pts = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        kp, Lam = 1.3, 0.07
        for pp in (0.4, 0.9, 2.1, 3.7):
            L = user1.r3_lines(c, kp, pp, Lam)
            u = L['L1'] + L['L2'] + L['L3'] - L['L4'] + L['L5']
            pts.append(rel(u, fp1.row3(c, kp, pp, Lam, sgt=-1)))
out['RowIII_new_integrand_vs_first_principles_E3'] = max(pts)

# 3./4. N_2 with the overall minus; green terms
n2, green = [], []
for sym, seed, kp, Lam, vee in ((True, 1, 1.3, 0.07, 9.0), (True, 2, 1.3, 0.07, 9.0), (False, 1, 0.8, 0.2, 3.0), (False, 6, 2.0, 0.01, 40.0)):
    c = Cfg(seed, sym=sym, names=NAMES)
    q = user1.r3_kin(c)
    A_, B_ = q['A'], q['B']
    ps = kp * B_ / A_
    N = user1.r3_int(c, kp, Lam, vee)
    N2_new = -sum(N[k] for k in N if k.startswith('N2'))
    g = lambda p: -user1.r3_lines(c, kp, p, Lam, strip=True)['L4']
    n2.append(rel(N2_new, pv(g, kp + Lam, vee, ps, A_) / (2 * PI)))
    rat_new = N['N1.1'] + N['N1.4'] - N['N2.1'] + N['N3a.1']
    green.append(float(abs(rat_new) / abs(N['N1.1'])))
out['RowIII_new_N2_vs_PV_integral'] = max(n2)
out['RowIII_green_terms_sum_over_one_term'] = max(green)

# 5. Row V 2b finite: 'Simplifying' (B2 pieces with U(z), as written) and the final form (as written)
src = open(os.path.join(HERE, '..', 'rowV', 'part2_user_check.py')).read()
exec(src[:src.index('def run(')])            # cq, sq, setup, contractions (rng is re-seeded there)
r5 = []
for Lam_over_kp in (1e-3, 1e-5):
    for _ in range(4):
        x, y, z, xp, wp, k = setup()
        kp = 1.7
        Lam = Lam_over_kp * kp
        Lg = np.log((kp - Lam) / Lam)
        b = contractions(x, y, z, xp, wp)
        Ox, Oz = complex(*rng.normal(size=2)), complex(*rng.normal(size=2))
        ph_wz, ph_wx = np.exp(-1j * k @ (wp - z)), np.exp(-1j * k @ (wp - x))
        ph = lambda p: ph_wz * np.exp(-1j * k @ (z - x) * p / kp)
        phq = lambda q_: ph_wx * np.exp(-1j * k @ (x - z) * q_ / kp)
        D = lambda p: p * sq(y - x) + (kp - p) * sq(y - z)
        Dq = lambda q_: (kp - q_) * sq(y - x) + q_ * sq(y - z)
        c5 = lambda p: [p * (kp - p) / kp**2, -p / (kp - p), (kp - p) / kp, -(kp - p) / p, p / kp, -1.0]
        Kyz = ((xp - wp) @ (y - z)) / (sq(xp - wp) * sq(y - z))
        S = ((y - x) @ (y - z)) / (sq(y - x) * sq(x - z))
        Phi = lambda p: p * (kp - p) / kp**2 + p / (kp - p) + (kp - p) / p
        exact = cq(lambda p: ph(p) / D(p) * sum(ci * bi for ci, bi in zip(c5(p), b)) * Ox
                   - (1 / kp) * Phi(p) * ph_wz * Kyz * S * Oz, Lam, kp - Lam)
        g2 = lambda q_: phq(q_) / Dq(q_) * b[1]
        f4 = lambda p: ph(p) / D(p) * b[3]
        reg_noB2 = lambda p: ph(p) / D(p) * (p * (kp - p) / kp**2 * b[0] + (kp - p) / kp * b[2] - (kp - p) / p * b[3]
                                             + p / kp * b[4] - b[5])
        CTfull = cq(lambda p: -(1 / kp) * Phi(p) * ph_wz * Kyz * S * Oz, Lam, kp - Lam)
        # 'Simplifying' as written: B2 pieces with O_z
        B2_written = -cq(lambda q_: (kp / q_ - 1) * g2(q_), Lam, kp - Lam)
        simp_written = cq(lambda p: reg_noB2(p) * Ox, Lam, kp - Lam) + CTfull + B2_written * Oz
        simp_fixed = cq(lambda p: reg_noB2(p) * Ox, Lam, kp - Lam) + CTfull + B2_written * Ox
        # final form as written in the new file
        regb = lambda p: ph(p) / D(p) * (p * (kp - p) / kp**2 * b[0] + (kp - p) / kp * b[2] + b[3] + p / kp * b[4] - b[5])
        P1 = cq(lambda p: (regb(p) - kp / p * (f4(p) - f4(0))) * Ox, Lam, kp - Lam)
        P2 = -cq(lambda q_: kp / q_ * (g2(q_) - g2(0)) * Ox, 0, kp)
        P3 = cq(lambda q_: g2(q_) * Ox, Lam, kp - Lam)
        E2 = -Lg * ph_wx * b[1] / sq(y - x) * Ox
        E4 = -Lg * ph_wz * b[3] / sq(y - z) * Ox
        CT = -(2 * Lg - 11 / 6) * ph_wz * Kyz * S * Oz
        final = P1 + P2 + P3 + E2 + E4 + CT
        r5.append(dict(Lam_over_kp=Lam_over_kp, simplifying_written=rel(simp_written, exact),
                       simplifying_fixed=rel(simp_fixed, exact), final_written=rel(final, exact)))
out['RowV_2b_finite'] = r5

# 6. beta_0: Row II part 1 (11/3 term) + Row V 2b UV (-11/6 term), per ordering, after y -> x, z -> w
b0 = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        U, rho = c.U, c.rho
        # colour factors at the same points (x: source, w: gluon, x', w')
        UX, UW, UXP, UWP = U['x'], U['w'], U['xp'], U['wp']
        rx, rxp = rho['x'], rho['xp']
        C_II = ((UXP - UWP) @ rxp) @ ((UW - UX) @ rx)            # [U(x')-U(w')] rho(x') [U(w)-U(x)] rho(x)
        C_V = ((UXP - UWP) @ rxp) @ ((UX - UW) @ rx)             # [U(x')-U(w')] [U(y)-U(z)] rho(x') rho(y), y->x, z->w
        E = rng.normal()                                          # 2/eps + 2 gamma + log((x-w)^2 mu^2/4), any value
        II = G**4 * fp1.NC / (32 * PI**4) * (-E) * (11 / 3) * C_II
        V = G**4 * fp1.NC / (16 * PI**4) * (-11 / 6) * E * C_V
        b0.append(float(abs(II + V) / abs(II)))
out['beta0_RowII_plus_RowV2b_per_ordering'] = max(b0)
dl, x0 = sp.symbols('delta xi0', positive=True)
B_II_nonlog = sp.Rational(11, 3) - 8 * dl + 2 * dl**2 - sp.Rational(4, 3) * dl**3      # G5a exact bracket, non-log part
B_V_nonlog = -sp.Rational(11, 6) + 4 * dl - dl**2 + sp.Rational(2, 3) * dl**3          # 2b UV exact bracket, non-log part
out['beta0_and_powers_exact_in_Lambda'] = str(sp.simplify(B_II_nonlog + 2 * B_V_nonlog))

json.dump(out, open(os.path.join(HERE, 'check_v2.json'), 'w'), indent=1)
for k_, v_ in out.items():
    print(k_, ':', v_)
