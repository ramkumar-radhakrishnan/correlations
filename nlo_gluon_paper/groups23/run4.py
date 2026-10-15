"""C rows: w' integration and plus prescription."""
import numpy as np
from scipy import integrate
import fp
import user23 as us
import user23c as uc
import fixed23 as fx

k, Lam, vee = 1.7, 0.11, 23.0
I = 1j


def quad(f, a, b):
    re = integrate.quad(lambda p: np.real(f(p)), a, b, limit=800, epsabs=0, epsrel=1e-11)[0]
    im = integrate.quad(lambda p: np.imag(f(p)), a, b, limit=800, epsabs=0, epsrel=1e-11)[0]
    return re + 1j * im


def wstar(c, p):
    P = p + k
    return (P * c.pt['yp'] - p * c.pt['z']) / k


for sym in (True,):
    for seed in (1, 2, 3):
        c = fp.Cfg(seed, sym)
        p0 = 0.63
        # ---- w' integration: written (lines 72, 85, 105, 124) vs pre-integration line at w' = w'*
        col68 = np.einsum('CDa,CDb,ab->', fp.F, fx.Lc_fix(c), fx.two_fix(c, 'w', 'y', 'x'))
        col80 = np.einsum('CDa,CDb,ab->', fp.F, fx.Lc_fix(c), fx.Rone_fix(c, 'w', 'z', 'x'))
        col101 = np.einsum('CDa,CDb,ab->', fp.F, fx.Lc_fix(c), fx.X1amp_fix(c))
        col120 = np.einsum('CDa,CDb,ab->', fp.F, fx.Lc_fix(c), fx.X2amp_fix(c))
        c.pt['wp'] = wstar(c, p0)
        pre = {'III(1)': (fx.g3_3_part1_line68(c, k, p0, sign=1), uc.ker72(c, k, p0, -I * fp.G**4 / (8 * np.pi**4)) * col68),
               'III(2)': (fx.g3_3_part2_line80(c, k, p0, sign=1), (uc.ker85_pieces(c, k, p0)[0] + uc.ker85_pieces(c, k, p0)[1] / p0) / uc.phase(c, k, p0) * col80),
               'V': (fx.g3_5_line101(c, k, p0, sign=1), uc.ker72(c, k, p0, -I * fp.G**4 / (4 * np.pi**4)) * col101),
               'VI': (fx.g3_6_line120(c, k, p0, sign=1), uc.ker72(c, k, p0, -I * fp.G**4 / (4 * np.pi**4)) * col120)}
        for key, (a, b) in pre.items():
            print('seed %d %-7s after w\' integration: written / (pre-integration line at w\'=w\'*) = %s' % (seed, key, np.round(b / a, 10)))
        # phase check
        ph = np.exp(-1j * c.kv @ (wstar(c, p0) - c.pt['w']))
        ph2 = np.exp(-1j * c.kv @ (c.pt['yp'] - c.pt['w'])) * uc.phase(c, k, p0)
        print('   phase e^{-ik(w\'*-w)} / [e^{-ik(y\'-w)} e^{-ik(y\'-z)p/k}] = %s' % np.round(ph / ph2, 12))
        # ---- plus prescription
        for key, pref, col in (('III(1)', -I * fp.G**4 / (8 * np.pi**4), col68), ('V', -I * fp.G**4 / (4 * np.pi**4), col101),
                               ('VI', -I * fp.G**4 / (4 * np.pi**4), col120)):
            full = lambda p: sum(uc.ker72_pieces(c, k, p, pref)[:1]) + uc.ker72_pieces(c, k, p, pref)[1] / p
            exact = quad(full, Lam, vee - k) / (2 * np.pi)
            g10 = uc.ker72_pieces(c, k, 0.0, pref)[1]
            sub = lambda p: uc.ker72_pieces(c, k, p, pref)[0] + (uc.ker72_pieces(c, k, p, pref)[1] - g10) / p
            preflog = {'III(1)': I * fp.G**4 / (16 * np.pi**5), 'V': I * fp.G**4 / (8 * np.pi**5), 'VI': I * fp.G**4 / (8 * np.pi**5)}[key]
            logt = uc.logterm72(c, k, Lam, vee, preflog)
            written = quad(sub, 1e-14, vee - k) / (2 * np.pi) + logt
            fixed = quad(sub, Lam, vee - k) / (2 * np.pi) + logt
            print('seed %d %-7s plus prescription: written(int from 0)/exact = %s   with lower limit Lambda: %s' % (
                seed, key, np.round(written / exact, 8), np.round(fixed / exact, 12)))
        reg = lambda p: uc.ker85_pieces(c, k, p)[0]
        g1 = lambda p: uc.ker85_pieces(c, k, p)[1]
        exact = quad(lambda p: reg(p) + g1(p) / p, Lam, vee - k) / (2 * np.pi)
        g10 = g1(0.0)
        sub = lambda p: reg(p) + (g1(p) - g10) / p
        logt = uc.logterm90(c, k, Lam, vee)
        written = quad(sub, 1e-14, vee - k) / (2 * np.pi) + logt
        fixed = quad(sub, Lam, vee - k) / (2 * np.pi) + logt
        print('seed %d III(2)  plus prescription: written(int from 0)/exact = %s   with lower limit Lambda: %s' % (
            seed, np.round(written / exact, 8), np.round(fixed / exact, 12)))
