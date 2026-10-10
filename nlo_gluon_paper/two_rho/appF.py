"""Singlet projection: <rho^a(p) rho^b(q)> = delta^{ab} mu(p,q)/(Nc^2-1), mu = <rho^f rho^f>.  The bracket of
Appendix F multiplies P = alpha_s^2 dY/(8 pi^6 (Nc^2-1)) mu, i.e. it equals the operator target with rho^a rho^b -> delta^{ab}.
Compare at one point
(p1, p2, w, w', z) the operator-level JIMWLK target with Appendix F (as written / corrected), and
Eq. (LOBFKL) of Appendix E.  All in units of P (the prefactor of Appendix F) x mu."""
import numpy as np
from core import Config, facs, placements, sub_kern, kernel_value, ORD
import target as TG

NC = 3


def singlet_value(terms, cfg):
    """sum of terms with the two rho's contracted with delta^{ab}/(Nc^2-1); positions as in cfg"""
    tot = 0.0
    for t in terms:
        maps, wgt = placements(t)
        for mp in maps:
            g = lambda p, mp=mp: mp.get(p, p)
            fl = facs(t, mp)
            rh = [f for f in fl if f[0] == 'rho']
            rest = [f for f in fl if f[0] != 'rho']
            fl2 = rest + [('d', rh[0][1], rh[1][1])]
            col = cfg.colour(fl2)
            tot += t['c'] * wgt * col * kernel_value(sub_kern(t['kern'], g), cfg.P)
    return tot


def K(cfg, a, b, t='z'):
    P = cfg.P
    A, B = P[a] - P[t], P[b] - P[t]
    return (A @ B) / ((A @ A) * (B @ B))


def KK(cfg, a, b, c, d):
    P = cfg.P
    A, B = P[a] - P[b], P[c] - P[d]
    return (A @ B) / ((A @ A) * (B @ B))


def appF(cfg, x, y, w='w', wp="w'", t='z', corrected=True):
    """bracket of Appendix F at fixed x (tied to w) and y (tied to w'), without the prefactor P"""
    U = {k: v for k, v in cfg.U.items()}
    Ud = lambda p: U[p].T
    T = [-1j * cfg.f[e] for e in range(NC * NC - 1)]
    tr = np.trace
    TeAB = lambda A, B: sum(tr(T[e] @ A @ T[e] @ B) for e in range(len(T)))
    TT = lambda A, B: sum(tr(A @ T[d]) * tr(B @ T[d]) for d in range(len(T)))
    Mt = lambda u, v: (Ud(u) - Ud(t)) @ (U[v] - U[t])
    Mxy = (Ud(x) - Ud(t)) @ (U[y] - U[t])
    Myx = (Ud(y) - Ud(t)) @ (U[x] - U[t])
    d = np.zeros(5, complex)
    for v, sv in [(x, 1), (w, -1)]:
        for u, su in [(y, 1), (wp, -1)]:
            d[0] += sv * su * K(cfg, u, v, t) * TeAB(Ud(v) @ U[u], Mt(u, v) if corrected else Mxy)
    for u, su in [(x, 1), (w, -1)]:
        for v, sv in [(y, 1), (wp, -1)]:
            d[1] += su * sv * K(cfg, u, v, t) * TeAB(Ud(u) @ U[v], Mt(u, v).T if corrected else Myx)
    for v, sv in [(y, 1), (wp, -1)]:
        d[2] += -sv * K(cfg, v, v, t) * TeAB((Ud(x) - Ud(w)) @ U[v], Mt(v, v) if corrected else Mxy)
    for v, sv in [(x, 1), (w, -1)]:
        d[3] += -sv * K(cfg, v, v, t) * TeAB(Ud(v) @ (U[y] - U[wp]), Mt(v, v).T if corrected else Myx)
    s5 = 1 if corrected else -1
    for u, su in [(y, 1), (wp, -1)]:
        d[4] += s5 * su * K(cfg, u, u, t) * TT((Ud(u) - Ud(t)) @ U[u], (Ud(x) - Ud(w)) @ U[u])
    for u, su in [(x, 1), (w, -1)]:
        d[4] += -s5 * su * K(cfg, u, u, t) * TT((Ud(u) - Ud(t)) @ U[u], Ud(u) @ (U[y] - U[wp]))
    return d * KK(cfg, x, w, y, wp)


def lobfkl(cfg, x, y, z='z', w='w', wp="w'"):
    """bracket of Eq. (LOBFKL) at fixed x, y, z (prefactor -alpha_s^2 Nc dY/(8 pi^6 (Nc^2-1)))"""
    U = cfg.U
    Ud = lambda p: U[p].T
    tr = np.trace
    Kf = lambda a, b: K(cfg, a, b, z)
    KD = Kf(x, x) - 2 * Kf(x, y) + Kf(y, y)
    v = KK(cfg, x, w, y, wp) * KD * tr((Ud(x) - Ud(w)) @ (U[y] - U[wp]))
    v += -2 * KK(cfg, z, w, z, wp) * Kf(x, y) * tr((Ud(z) - Ud(w)) @ (U[z] - U[wp]))
    v += KK(cfg, z, w, y, wp) * (2 * Kf(x, y) - Kf(x, x)) * tr((Ud(z) - Ud(w)) @ (U[y] - U[wp]))
    v += KK(cfg, x, w, z, wp) * (2 * Kf(x, y) - Kf(y, y)) * tr((Ud(x) - Ud(w)) @ (U[z] - U[wp]))
    return v


if __name__ == '__main__':
    import json
    out = []
    for sym in (False, True):
        for seed in (1, 2, 3):
            cfg = Config(seed, sym)
            tgt = singlet_value(TG.jimwlk(), cfg)
            # Appendix F: x tied to w, y tied to w'; symmetrize over the charge positions p1, p2
            Fc = 0.5 * (appF(cfg, 'p2', 'p1') + appF(cfg, 'p1', 'p2'))
            Fw = 0.5 * (appF(cfg, 'p2', 'p1', corrected=False) + appF(cfg, 'p1', 'p2', corrected=False))
            B = 0.5 * (lobfkl(cfg, 'p2', 'p1') + lobfkl(cfg, 'p1', 'p2'))
            out.append(dict(sym=sym, seed=seed, target=tgt.real, corrected=Fc.sum().real, written=Fw.sum().real, bfkl=-NC * B.real,
                            pieces_corrected=[float(v.real) for v in Fc], pieces_written=[float(v.real) for v in Fw]))
            print('symmetric U' if sym else 'general U  ', 'seed', seed,
                  ' JIMWLK target %.6f | App F corrected %.6f | as written %.6f | -Nc*LOBFKL bracket %.6f'
                  % (tgt.real, Fc.sum().real, Fw.sum().real, -NC * B.real))
    json.dump(out, open('appF.json', 'w'), indent=1)
