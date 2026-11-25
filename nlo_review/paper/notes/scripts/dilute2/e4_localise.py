#!/usr/bin/env python3
"""dilute2 / e4: localisation of the mismatch.  Residual of the test  E = c*x + e0  (x = ln(qa^2/qb^2)) and
one-parameter extensions  E + kappa*T_j = c*x + e0  with candidate structures T_j (see list).  Reads the table printed by e3.
usage: e4_localise.py <dir new outputs> <dir ref3_dilute>"""
import sys, os, io, contextlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e3_assemble as e3
from e0_points import PTS
S, S0 = sys.argv[1], sys.argv[2]
abar = e3.abar
names = []
for n in PTS:
    try:
        with contextlib.redirect_stdout(io.StringIO()): e3.groupI(n)
        names.append(n)
    except Exception: pass
rows = []
for n in names:
    with contextlib.redirect_stdout(io.StringIO()): r = e3.groupI(n)
    P = PTS[n]; qa, qb = P['qa'], P['qb']; a2, b2, k2 = qa @ qa, qb @ qb, P['kn'] ** 2
    a, b, lk = np.log(a2), np.log(b2), np.log(k2); x = a - b
    X = abar / 4 * e3.X_bsv(qa, qb)
    E = r['C'] - X + abar / 4 * x * (2 * lk - a - b)
    den = a2 - b2
    T = {
        'conv: (abar/2)[ln k^2 (a-b) - (a^2-b^2)/2]': abar / 2 * (lk * x - (a * a - b * b) / 2),
        'proj only: (abar/2)[ln k^2 a - a^2/2]': abar / 2 * (lk * a - a * a / 2),
        'targ only: (abar/2)[ln k^2 b - b^2/2]': abar / 2 * (lk * b - b * b / 2),
        '(abar/4)(a-b)^2': abar / 4 * x * x,
        '(abar/4) k^2/(qa^2-qb^2) ln(qa^2/qb^2)': abar / 4 * k2 / den * x if abs(den) > 1e-9 else abar / 4 * k2 / a2,
        '(abar/4)*[rational bracket of X]': abar / 4 * (e3.X_bsv(qa, qb, False) - (3 * k2 / den * x - x * x)) if abs(den) > 1e-9 else np.nan,
        '(abar/4) ln^2(k^2/qa^2)': abar / 4 * (lk - a) ** 2,
        '(abar/4) ln^2(k^2/qb^2)': abar / 4 * (lk - b) ** 2,
        '(abar/4) ln(k^2/qa^2) ln(k^2/qb^2)': abar / 4 * (lk - a) * (lk - b),
        '(abar/4) [ln^2(k^2/qa^2)+ln^2(k^2/qb^2)]': abar / 4 * ((lk - a) ** 2 + (lk - b) ** 2),
        '(abar/4) [ln^2(k^2/qa^2)-ln^2(k^2/qb^2)]': abar / 4 * ((lk - a) ** 2 - (lk - b) ** 2),
    }
    rows.append((n, x, E, r['dC'], T, r['C'], X))
xs = np.array([r[1] for r in rows]); Es = np.array([r[2] for r in rows]); dE = np.array([r[3] for r in rows])
print("points:", names)
A = np.vstack([xs, np.ones_like(xs)]).T
(c, e0), *_ = np.linalg.lstsq(A, Es, rcond=None); res = Es - A @ [c, e0]
print(f"baseline E = c x + e0: c = {c:.5f}, e0 = {e0:.5f}, residuals {np.round(res, 5)}, rms {np.sqrt(np.mean(res ** 2)):.5f}  (errors of E: {np.round(dE, 5)})")
if len(rows) >= 4:
    for key in rows[0][4]:
        t = np.array([r[4][key] for r in rows])
        ok = np.isfinite(t)
        if ok.sum() < 4: continue
        B = np.vstack([t[ok], xs[ok], np.ones(ok.sum())]).T
        (kap, c, e0), *_ = np.linalg.lstsq(B, Es[ok], rcond=None); res = Es[ok] - B @ [kap, c, e0]
        print(f"  E - kappa*T, T = {key:45s}: kappa = {kap: .4f}  c = {c: .5f}  e0 = {e0: .5f}  rms residual {np.sqrt(np.mean(res ** 2)):.5f}  ({ok.sum()} points, 3 parameters)")
