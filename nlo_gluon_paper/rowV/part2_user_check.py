"""Check of the user's Row V part 2 (2a and 2b) after the plus prescription, at fixed random transverse points.

Common normalization: everything is multiplied by (2pi)^3 (8 pi^4/g^4) k+ (2 pi), i.e. the p+ integrand is
   phase(p+)/D(p+) * sum_j c_j(p+) b_j * O   (+ counterterm), integrated over dp+,
and an endpoint term (g^4/16 pi^5)(1/k+) [..] becomes [..] (no k+ factors).
Colour operators are p+-independent numbers at fixed points: O_a (2a), O_x = N_c[U(y)-U(x)], O_z = N_c[U(y)-U(z)].
-> part2_user_check.json"""
import json
import os

import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(2024)


def cq(f, a, b):
    return quad(lambda t: f(t).real, a, b, limit=400, epsabs=1e-13, epsrel=1e-12)[0] + \
        1j * quad(lambda t: f(t).imag, a, b, limit=400, epsabs=1e-13, epsrel=1e-12)[0]


def sq(v):
    return float(v @ v)


def setup():
    x, y, z, xp, wp = (rng.normal(size=2) for _ in range(5))
    k = rng.normal(size=2)
    return x, y, z, xp, wp, k


def contractions(x, y, z, xp, wp, dperp=2.0):
    """b_j = K^i(x'-w') (z-x)^m/(z-x)^2 B_j^{im}, j = 1..6 (scalars)"""
    K = (xp - wp) / sq(xp - wp)
    w = (z - x) / sq(z - x)
    yx, yz, xz = y - x, y - z, x - z
    I = np.eye(2)
    B1 = I * dperp / (2 * sq(xz)) * (sq(yx) - sq(yz))
    B2 = np.outer(yx, xz) / sq(xz) - np.outer(yx, yz) / (2 * sq(yz))
    B3 = I * (yz @ xz / sq(xz) + yx @ yz / (2 * sq(yx)) - (sq(yx) - sq(yz)) / (2 * sq(xz)))
    B4 = np.outer(yz, xz) / sq(xz) + np.outer(yz, yx) / (2 * sq(yx))
    B5 = I * (yx @ xz / sq(xz) - yx @ yz / (2 * sq(yz)) - (sq(yx) - sq(yz)) / (2 * sq(xz)))
    B6 = (np.outer(xz, yx) / sq(xz) - np.outer(yz, yx) / (2 * sq(yz)) + np.outer(xz, yz) / sq(xz)
          + np.outer(yx, yz) / (2 * sq(yx)))
    return [K @ B @ w for B in (B1, B2, B3, B4, B5, B6)]


def run(Lam_over_kp, kp=1.7, nrep=6):
    out = []
    for _ in range(nrep):
        x, y, z, xp, wp, k = setup()
        Lam = Lam_over_kp * kp
        L = np.log((kp - Lam) / Lam)
        lam = Lam / kp
        b = contractions(x, y, z, xp, wp)
        Oa, Ox, Oz = (complex(*rng.normal(size=2)) for _ in range(3))
        ph_wz = np.exp(-1j * k @ (wp - z))
        ph_wx = np.exp(-1j * k @ (wp - x))
        ph = lambda p: ph_wz * np.exp(-1j * k @ (z - x) * p / kp)                     # e^{-ik(w'-z)} e^{-ik(z-x)p+/k+}
        phq = lambda q: ph_wx * np.exp(-1j * k @ (x - z) * q / kp)                    # same, with p+ = k+ - q+
        D = lambda p: p * sq(y - x) + (kp - p) * sq(y - z)
        Dq = lambda q: (kp - q) * sq(y - x) + q * sq(y - z)
        c = lambda p: [p * (kp - p) / kp**2, -p / (kp - p), (kp - p) / kp, -(kp - p) / p, p / kp, -1.0]
        Kyz = ((xp - wp) @ (y - z)) / (sq(xp - wp) * sq(y - z))
        S = ((y - x) @ (y - z)) / (sq(y - x) * sq(x - z))
        Phi = lambda p: p * (kp - p) / kp**2 + p / (kp - p) + (kp - p) / p

        # ---------------- 2a: original vs. the user's plus form (user's limits 0..k+)
        orig_a = cq(lambda p: ph(p) / D(p) * sum(ci * bi for ci, bi in zip(c(p), b)) * Oa, Lam, kp - Lam)
        f4 = lambda p: ph(p) / D(p) * b[3]
        g2 = lambda q: phq(q) / Dq(q) * b[1]
        reg = lambda p: ph(p) / D(p) * (p * (kp - p) / kp**2 * b[0] + (kp - p) / kp * b[2] + b[3] + p / kp * b[4] - b[5])
        blk1 = cq(lambda p: reg(p) - kp / p * (f4(p) - f4(0)), 0, kp)
        blk2 = -cq(lambda q: kp / q * (g2(q) - g2(0)) - g2(q), 0, kp)                  # int_{k+}^{0} {[k+/q+]_+ - 1} g2
        E4a = -L * ph_wz * b[3] / sq(y - z)
        E2a = -L * ph_wx * b[1] / sq(y - x)
        user_a = (blk1 + blk2) * Oa + (E4a + E2a) * Oa

        # ---------------- 2b finite: original (2b - CT) on [Lambda, k+-Lambda]
        orig_b = cq(lambda p: ph(p) / D(p) * sum(ci * bi for ci, bi in zip(c(p), b)) * Ox
                    - (1 / kp) * Phi(p) * ph_wz * Kyz * S * Oz, Lam, kp - Lam)

        # correct plus form (CT integrated over p+): P1 + P2 + P3 + E4 + E2 + CT
        regb = lambda p: ph(p) / D(p) * (p * (kp - p) / kp**2 * b[0] + (kp - p) / kp * b[2] + b[3] + p / kp * b[4] - b[5])
        P1 = cq(lambda p: (regb(p) - kp / p * (f4(p) - f4(0))) * Ox, Lam, kp - Lam)
        P23 = -cq(lambda q: (kp / q * (g2(q) - g2(0)) - g2(q)) * Ox, Lam, kp - Lam)
        E4 = -L * ph_wz * b[3] / sq(y - z) * Ox
        E2 = -L * ph_wx * b[1] / sq(y - x) * Ox
        CTint = -(2 * L - 11 / 6 + 4 * lam - lam**2 + 2 * lam**3 / 3) * ph_wz * Kyz * S * Oz
        correct_b = P1 + P23 + E4 + E2 + CTint

        # the user's final form, as written
        P23_user = -cq(lambda q: (kp / q * (g2(q) - g2(0))) * Oz, 0, kp) + cq(lambda q: g2(q) * Oz, Lam, kp - Lam)
        E2_user = -L * ph_wx * b[1] / sq(y - x) * Oz
        E4_user = -L * ph_wx * b[3] / sq(y - z) * Oz
        CT_user = -(11 / 6 - L / 3) * ph_wz * Kyz * S * Oz
        user_b = P1 + P23_user + E2_user + E4_user + CT_user

        # the user's form with each fix applied in turn
        fix1 = P1 + (-cq(lambda q: (kp / q * (g2(q) - g2(0))) * Ox, 0, kp) + cq(lambda q: g2(q) * Ox, Lam, kp - Lam)) \
            + (-L * ph_wx * b[1] / sq(y - x) * Ox) + E4_user + CT_user                          # colour O_x in the B2 pieces
        fix2 = fix1 - E4_user + E4                                                                # B4 endpoint: phase w'-z, colour O_x
        fix3 = fix2 - CT_user + CTint                                                             # counterterm coefficient
        out.append(dict(
            rel_2a_user=abs(user_a - orig_a) / abs(orig_a),
            rel_2b_correct=abs(correct_b - orig_b) / abs(orig_b),
            rel_2b_user=abs(user_b - orig_b) / abs(orig_b),
            rel_2b_fix_colourB2=abs(fix1 - orig_b) / abs(orig_b),
            rel_2b_fix_colourB2_E4=abs(fix2 - orig_b) / abs(orig_b),
            rel_2b_fix_all=abs(fix3 - orig_b) / abs(orig_b)))
    return out


R = {}
R['Lambda/k+=1e-3'] = run(1e-3)
R['Lambda/k+=1e-5'] = run(1e-5)
R['Lambda/k+=0.05'] = run(0.05)
json.dump(R, open(os.path.join(HERE, 'part2_user_check.json'), 'w'), indent=1)
for k, rows in R.items():
    print(k)
    for key in rows[0]:
        print('   %-26s max %.2e  min %.2e' % (key, max(r[key] for r in rows), min(r[key] for r in rows)))
