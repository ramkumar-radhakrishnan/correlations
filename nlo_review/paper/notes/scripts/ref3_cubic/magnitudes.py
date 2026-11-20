# Referee 3, Major 7 (a): size of the symmetrised (k-odd) and contact (k-even) colour weights of the cubic terms,
# relative to the quadratic (LO-type) colour weight, for a colour-singlet dipole (exact charge space), for dense
# (Haar-random) and weak (V = exp(i eps alpha.t)) targets.  Colour weights only: the kernels are not included.
# usage: python3 -I magnitudes.py
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colour import SUN, structures, ORDER, POINTS
from quark_dipole import State, exact, split


def lo_weight(st, U, sites):
    # LO colour operator [U(y')-U(x')]^{ce'}[U(y)-U(x)]^{ce} rho^{e'}(x') rho^e(x), rho points (x', x)
    C = np.einsum('ca,cb->ab', U['yp'] - U['xp'], U['y'] - U['x'])
    s1, s2 = sites
    return np.einsum('ab,ab->', C, st.E2[s1, :, s2, :]).real


def sample(G, st, rng, eps, nsamp):
    Vsite = None
    out = {name: [[], []] for name in ORDER}; lo = []
    for _ in range(nsamp):
        Vs = [G.adj(G.randV(rng, eps)) for _ in range(st.P)]
        for name in ORDER:
            rho_pts = structures(G, {p: np.eye(G.n) for p in POINTS})[name][2]
            for sites in itertools.product(range(st.P), repeat=3):
                U = {p: G.adj(G.randV(rng, eps)) for p in POINTS}
                for p, s in zip(rho_pts, sites): U[p] = Vs[s]
                C, typ, _ = structures(G, U)[name]
                sy, co = split(st, C, typ, sites)
                out[name][0].append(sy); out[name][1].append(abs(co))
        for sites in itertools.product(range(st.P), repeat=2):
            U = {p: G.adj(G.randV(rng, eps)) for p in POINTS}
            U['xp'], U['x'] = Vs[sites[0]], Vs[sites[1]]
            lo.append(lo_weight(st, U, sites))
    rms = lambda a: np.sqrt(np.mean(np.abs(np.array(a)) ** 2))
    return {k: (rms(v[0]), rms(v[1])) for k, v in out.items()}, rms(lo)


def main():
    rng = np.random.default_rng(2024)
    print('Colour-singlet dipole: rms over target configurations and charge assignments of |sym| (k-odd, C-odd) and |contact| (k-even),')
    print('in units of the rms LO (quadratic) colour weight.  Dense target: Haar-random V at every point.')
    for N in (3, 4, 5, 6):
        G = SUN(N); st = State(G, 'dipole')
        r, lo = sample(G, st, rng, None, 6 if N < 5 else 3)
        ssym = np.mean([v[0] for v in r.values()]) / lo; scon = np.mean([v[1] for v in r.values()]) / lo
        print('  SU(%d): LO rms %.3f ; mean over the 13 structures: sym/LO %.3f, contact/LO %.3f ;  per structure sym/LO: %s' % (
            N, lo, ssym, scon, ' '.join('%.2f' % (v[0] / lo) for v in r.values())))
        print('         per structure contact/LO: %s' % ' '.join('%.2f' % (v[1] / lo) for v in r.values()))
    print('Weak target, SU(3): scaling with eps (V = exp(i eps alpha.t), alpha Gaussian)')
    G = SUN(3); st = State(G, 'dipole'); prev = None
    for eps in (0.4, 0.2, 0.1, 0.05):
        r, lo = sample(G, st, np.random.default_rng(77), eps, 4)
        sy = np.mean([v[0] for v in r.values()]); co = np.mean([v[1] for v in r.values()])
        line = '  eps=%.3f: LO %.3e  sym %.3e  contact %.3e   sym/LO %.3e  contact/LO %.3f' % (eps, lo, sy, co, sy / lo, co / lo)
        if prev: line += '   ratios to eps*2: LO %.2f sym %.2f contact %.2f' % (prev[0] / lo, prev[1] / sy, prev[2] / co)
        print(line); prev = (lo, sy, co)


if __name__ == '__main__':
    main()
