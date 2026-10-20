"""Check of the user's revised 2b finite (after E1-E3): same setup as part2_user_check.py. -> part2_user_check_v2.json"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'part2_user_check.py')).read()
exec(src[:src.index('def run(')])            # cq, sq, setup, contractions, rng


def run(Lam_over_kp, kp=1.7, nrep=6):
    out = []
    for _ in range(nrep):
        x, y, z, xp, wp, k = setup()
        Lam = Lam_over_kp * kp
        L = np.log((kp - Lam) / Lam)
        lam = Lam / kp
        b = contractions(x, y, z, xp, wp)
        Ox, Oz = (complex(*rng.normal(size=2)) for _ in range(2))
        ph_wz = np.exp(-1j * k @ (wp - z)); ph_wx = np.exp(-1j * k @ (wp - x))
        ph = lambda p: ph_wz * np.exp(-1j * k @ (z - x) * p / kp)
        phq = lambda q: ph_wx * np.exp(-1j * k @ (x - z) * q / kp)
        D = lambda p: p * sq(y - x) + (kp - p) * sq(y - z)
        Dq = lambda q: (kp - q) * sq(y - x) + q * sq(y - z)
        c = lambda p: [p * (kp - p) / kp**2, -p / (kp - p), (kp - p) / kp, -(kp - p) / p, p / kp, -1.0]
        Kyz = ((xp - wp) @ (y - z)) / (sq(xp - wp) * sq(y - z))
        S = ((y - x) @ (y - z)) / (sq(y - x) * sq(x - z))
        Phi = lambda p: p * (kp - p) / kp**2 + p / (kp - p) + (kp - p) / p
        f4 = lambda p: ph(p) / D(p) * b[3]
        g2 = lambda q: phq(q) / Dq(q) * b[1]
        exact = cq(lambda p: ph(p) / D(p) * sum(ci * bi for ci, bi in zip(c(p), b)) * Ox
                   - (1 / kp) * Phi(p) * ph_wz * Kyz * S * Oz, Lam, kp - Lam)
        regb = lambda p: ph(p) / D(p) * (p * (kp - p) / kp**2 * b[0] + (kp - p) / kp * b[2] + b[3] + p / kp * b[4] - b[5])
        P1 = cq(lambda p: (regb(p) - kp / p * (f4(p) - f4(0))) * Ox, Lam, kp - Lam)
        P2 = -cq(lambda q: kp / q * (g2(q) - g2(0)) * Ox, 0, kp)            # int_{k+}^{0} [k+/q+]_+
        P3 = cq(lambda q: g2(q) * Ox, Lam, kp - Lam)                          # - int_{k+-Lam}^{Lam} (...)
        E2 = -L * ph_wx * b[1] / sq(y - x) * Ox
        E4 = -L * ph_wz * b[3] / sq(y - z) * Ox
        CT_user = -(2 * L - 11 / 6) * ph_wz * Kyz * S * Ox                   # as written: colour U(x)
        CT_fix = -(2 * L - 11 / 6) * ph_wz * Kyz * S * Oz                    # colour U(z)
        CT_exact = -(2 * L - 11 / 6 + 4 * lam - lam**2 + 2 * lam**3 / 3) * ph_wz * Kyz * S * Oz
        user = P1 + P2 + P3 + E2 + E4 + CT_user
        fixed = P1 + P2 + P3 + E2 + E4 + CT_fix
        out.append(dict(rel_user=abs(user - exact) / abs(exact), rel_fixed=abs(fixed - exact) / abs(exact),
                        rel_fixed_with_lambda_terms=abs(fixed - CT_fix + CT_exact - exact) / abs(exact)))
    return out


R = {lab: run(v) for lab, v in (('Lambda/k+=0.05', 0.05), ('Lambda/k+=1e-3', 1e-3), ('Lambda/k+=1e-5', 1e-5))}
json.dump(R, open(os.path.join(HERE, 'part2_user_check_v2.json'), 'w'), indent=1)
for lab, rows in R.items():
    print(lab, {k: '%.1e..%.1e' % (min(r[k] for r in rows), max(r[k] for r in rows)) for k in rows[0]})
