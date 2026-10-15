"""Row V Wilson-line simplification (line 113) and Group II Row IV (lines 98, 104, tensor z integral)."""
import numpy as np
from scipy import integrate
import fp
import user23 as us
import fixed23 as fx
from fp import F, NC, G, PI, K, sq

k = 1.7
for sym in (True, False):
    print('==== U symmetric' if sym else '==== general U')
    for seed in (1, 2):
        c = fp.Cfg(seed, sym)
        U, r = c.U, c.rho
        # ---- Row V: line 109 colour vs line 113 colour (as written, and with corrected placements)
        col109 = np.einsum('CDa,CDb,ab->', F, us.Lc_written(c),
                           np.einsum('bc,ce,ad,e,d->ab', U['z'], U['x'], U['y'], r['x'], r['y'])
                           - np.einsum('ad,ce,cb,e,d->ab', U['w'], U['x'], U['z'], r['x'], r['y']))
        L113 = np.einsum('EC,E->C', U['yp'], r['xp']) - np.einsum('CE,E->C', U['xp'], r['xp'])
        R113 = np.einsum('ce,ad,e,d->ca', U['x'], U['y'], r['x'], r['y']) - np.einsum('ad,ce,e,d->ca', U['w'], U['x'], r['x'], r['y'])
        col113 = np.einsum('Cca,C,ca->', F, L113, R113)
        col109f = np.einsum('CDa,CDb,ab->', F, fx.Lc_fix(c), fx.X1amp_fix(c))
        L113f = np.einsum('CE,E->C', U['yp'], r['xp']) - np.einsum('CE,E->C', U['xp'], r['xp'])
        col113f = np.einsum('Cca,C,ca->', F, L113f, R113)
        print('  seed %d Row V line 113/line 109 colour: as written %s, corrected placements %s'
              % (seed, np.round(col113 / col109, 10), np.round(col113f / col109f, 10)))
        # ---- Row IV: line 98 vs line 92 at w*, w'* (Jacobian (P/k)^4)
        p = 0.63; P = p + k
        ws = (P * c.pt['y'] - p * c.pt['z']) / k
        wps = (P * c.pt['yp'] - p * c.pt['z']) / k
        c.pt['w'] = ws; c.pt['wp'] = wps
        l92 = fx.g2_4_line92(c, k, p) * (P / k)**4
        km, kk, kmp, kkp = K(c.pt['y'] - c.pt['z']), K(c.pt['x'] - c.pt['y']), K(c.pt['yp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp'])
        I2 = np.eye(2)
        d = 2
        br = (d * p / P**2 * np.einsum('km,KM->kmKM', I2, I2) + 2 / k * (np.einsum('mK,Mk->kmKM', I2, I2) - np.einsum('MK,km->kmKM', I2, I2))
              + (p / k**2 + 1 / p) * np.einsum('mM,kK->kmKM', I2, I2))
        kern98 = np.einsum('m,k,M,K,kmKM->', km, kk, kmp, kkp, br)
        Lp = np.einsum('CDa,CDb->ab', F, fx.Lc_fix(c))
        Ra = (np.einsum('cda,ce,db,e->ab', F, U['y'], U['z'], r['x']) - np.einsum('cda,db,ce,e->ab', F, U['z'], U['x'], r['x']))
        col92 = np.sum(Lp * Ra)
        l98 = 1 / us.pi3 * G**4 / (4 * PI**4) / k * kern98 * col92
        phase_ok = np.exp(-1j * c.kv @ (wps - ws)) / np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['y']) * P / k)
        col104f = NC * np.einsum('cE,E,ce,e->', U['yp'] - U['xp'], r['xp'], U['y'] - U['x'], r['x'])
        col104w = NC * np.einsum('Ec,E,ec,e->', U['yp'], r['xp'], U['y'], r['x']) \
            - NC * np.einsum('Ec,E,ce,e->', U['yp'], r['xp'], U['x'], r['x']) \
            - NC * np.einsum('Ec,E,ec,e->', U['xp'], r['xp'], U['y'], r['x']) \
            + NC * np.einsum('Ec,E,ce,e->', U['xp'], r['xp'], U['x'], r['x'])
        print('  seed %d Row IV line 98 / (line 92 x (P/k)^4 at w*, w\'*) = %s; phase %s' % (seed, np.round(l98 / l92, 10), np.round(phase_ok, 10)))
        print('  seed %d Row IV colour of line 104 / line 98: as written %s, corrected %s'
              % (seed, np.round(col104w / col92, 10), np.round(col104f / col92, 10)))

# ---- tensor z integral J^{ab}(y,y') = int d^2z (y-z)^a (y'-z)^b / ((y-z)^2 (y'-z)^2): traceless part
rng = np.random.default_rng(5)
y, yp = rng.normal(size=2), rng.normal(size=2)
s = y - yp
cen = (y + yp) / 2


def integrand(rr, th, a, b):
    z = cen + rr * np.array([np.cos(th), np.sin(th)])
    u, v = y - z, yp - z
    T = np.outer(u, v) / (sq(u) * sq(v))
    T = T - np.trace(T) / 2 * np.eye(2)
    return T[a, b] * rr


out = np.zeros((2, 2))
for a in range(2):
    for b in range(2):
        val = 0.0
        edges = [0, 0.5 * np.sqrt(sq(s)), np.sqrt(sq(s)), 3, 30, 300, 3000]
        for lo, hi in zip(edges[:-1], edges[1:]):
            val += integrate.dblquad(lambda th, rr: integrand(rr, th, a, b), lo, hi, 0, 2 * np.pi, epsabs=1e-10, epsrel=1e-10)[0]
        out[a, b] = val
pred = -np.pi * (np.outer(s, s) / sq(s) - np.eye(2) / 2)
print('traceless part of J^{ab}: numeric', np.round(out, 6).tolist(), ' predicted -pi(s^a s^b/s^2 - delta/2)', np.round(pred, 6).tolist())
