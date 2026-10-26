"""Group I Row III of your finite terms (the row minus its log(vee/k+) terms, vee -> infinity, Lambda -> 0; closed form
../finite/simp_exact.py, checked against the file in ../finite_v2), near the two special places:
 1. the circle (x-z)^2 = (x-w)^2 of the edge kernel: R stays finite (the 1/[(x-z)^2-(x-w)^2]^n cancel among the terms);
 2. z -> w: the angular average <R> ~ log(1/r)/r^2 (r = |z-w|): not integrable.  Without the last term of d^3N_3b,
    -(i g^4 f/16 pi^5)(1/k+) log[(k+(x-w)-(k+-Lambda)(x-z))^2/(k+(x-w)-Lambda(x-z))^2] (x'-w').(y-z)/((x'-w')^2(y-z)^2(z-w)^2),
    it is integrable.  Output: edge_check.json"""
import copy
import json
import numpy as np
import env  # noqa: F401
import split as SP
import simp_exact as SE
import ft
from ft import PI3, Kv, vec
from fp import G, PI, sq

KP, V = 1.0, 1e10


def R(c):
    return SE.g1_iii(c, KP, V) - SP.C_g1_iii(c, KP) * np.log(V / KP)


def n3b_last(c):
    Xp, Yz = Kv(c, 'xp', 'wp'), Kv(c, 'y', 'z')
    zw2, B = sq(vec(c, 'z', 'w')), sq(vec(c, 'x', 'w'))
    return (-1j * G**4 / (16 * PI**5) / KP * np.log(zw2 / B) * (Xp @ Yz) / zw2
            * ft.phase(c, 'w', 'wp') / PI3 * ft.col_G1_III(c))


if __name__ == '__main__':
    out = {'circle': {}, 'z_to_w': {}}
    for seed in (3, 5):
        c0 = SP.cfg(seed)
        x, w = c0.pt['x'], c0.pt['w']
        R0 = np.linalg.norm(x - w)
        rows = []
        for d in (1e-1, 1e-2, 1e-3, 1e-4, -1e-4, -1e-3, -1e-2):
            c = copy.deepcopy(c0)
            c.pt['z'] = x + R0 * (1 + d) * np.array([np.cos(0.3), np.sin(0.3)])
            rows.append(dict(delta=d, absR=abs(R(c))))
        out['circle'][seed] = rows
        rz = []
        for r in (1e-1, 1e-2, 1e-3, 1e-4):
            v1, v2 = [], []
            for th in (np.arange(64) + 0.5) * 2 * np.pi / 64:
                c = copy.deepcopy(c0)
                c.pt['z'] = w + r * np.array([np.cos(th), np.sin(th)])
                v1.append(R(c))
                v2.append(R(c) - n3b_last(c))
            rz.append(dict(r=r, r2_avg=abs(np.mean(v1)) * r * r, r2_avg_over_log=abs(np.mean(v1)) * r * r / np.log(1 / r),
                           without_last_N3b_r2_avg=abs(np.mean(v2)) * r * r, without_last_N3b_r_avg=abs(np.mean(v2)) * r))
        out['z_to_w'][seed] = rz
        print('seed', seed)
        for q in rz:
            print('  r=%.0e  r^2<R>=%.3e  r^2<R>/log(1/r)=%.3e   without last N3b term: r^2<R>=%.2e r<R>=%.2e'
                  % (q['r'], q['r2_avg'], q['r2_avg_over_log'], q['without_last_N3b_r2_avg'], q['without_last_N3b_r_avg']))
        print('  circle |R|:', ', '.join('%.3e' % q['absR'] for q in rows))
    json.dump(out, open('edge_check.json', 'w'), indent=1)
