"""Behaviour of the three-rho rows near z -> w (classical charges = the symmetric three-rho part).
f(r) = sum over the three-rho rows at z = w + r (cos t, sin t), averaged over t.  r^2 f(r) -> const != 0 means a
non-integrable 1/(z-w)^2 singularity (log divergence of int_z)."""
import json
import copy
import numpy as np
import env  # noqa: F401
import split as SP
import user_lvk as UL

THREE_USER = ['G1-III', 'G1-IV', 'G2-I(1)', 'G3-I', 'G3-II']
THREE_EXACT = {'G1-III': SP.C_g1_iii, 'G1-IV': SP.C_g1_iv, 'G2-I(1)': SP.C_g2_i1, 'G3-I': SP.C_g3_i, 'G3-II': SP.C_g3_ii}


def at(c, r, t):
    d = copy.deepcopy(c)
    d.pt['z'] = c.pt['w'] + r * np.array([np.cos(t), np.sin(t)])
    return d


def scan(funcs, seed, rs=(1e-2, 1e-3, 1e-4), nt=64):
    c = SP.cfg(seed)
    out = []
    for r in rs:
        ts = (np.arange(nt) + 0.5) * 2 * np.pi / nt
        vals = [sum(f(at(c, r, t), 1.0) for f in funcs) for t in ts]
        out.append(dict(r=r, r2_avg=abs(r**2 * np.mean(vals)), r1_avg=abs(r * np.mean(vals))))
    return out


if __name__ == '__main__':
    res = {}
    for s in (1, 2, 3):
        u = scan([UL.ROWS[k] for k in THREE_USER], s)
        e = scan(list(THREE_EXACT.values()), s)
        res[s] = dict(user=u, exact=e)
        print('seed', s)
        for a, b in zip(u, e):
            print('   r=%.0e   your list: r^2<f> %.3e   exact: r^2<f> %.3e   r<f> %.3e' % (a['r'], a['r2_avg'], b['r2_avg'], b['r1_avg']))
    json.dump(res, open('zw_sing.json', 'w'), indent=1)
