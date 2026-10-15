"""Integrated lines: written vs numerical p+ integral of the (verified) integrands, Lambda < p+ < vee - k+."""
import numpy as np
from scipy import integrate
import fp
import user23 as us
import user23b as ub
import fixed23 as fx

k, Lam, vee = 1.7, 0.11, 23.0


def pint(f):
    re = integrate.quad(lambda p: np.real(f(p)), Lam, vee - k, limit=400, epsabs=0, epsrel=1e-12)[0]
    im = integrate.quad(lambda p: np.imag(f(p)), Lam, vee - k, limit=400, epsabs=0, epsrel=1e-12)[0]
    return (re + 1j * im) / (2 * np.pi)


for sym in (True, False):
    print('==== U symmetric' if sym else '==== general U')
    for seed in (1, 2):
        c = fp.Cfg(seed, sym)
        # G2-I part 1
        A_ex = pint(lambda p: fx.g2_1_termA(c, k, p)); B_ex = pint(lambda p: fx.g2_1_termB(c, k, p))
        A_w, B_w = ub.g2_1_line50(c, k, Lam, vee)
        print('  seed %d G2-I(1) line 50: block A written/exact = %s, block B = %s' % (seed, np.round(A_w / A_ex, 10), np.round(B_w / B_ex, 10)))
        # G2-I part 2: groups
        col = fx.g2_2_colour58(c)
        w = {'2a': lambda p: (p + k)**2 / (p * k), '2b': lambda p: (p + k)**2 / p**2, '2c': lambda p: (p + k)**2 / k**2,
             '2d': lambda p: (p + k) / p, '2e': lambda p: (p + k) / k, '2f': lambda p: 1.0}
        gname = {'2a': 'pk', '2b': 'pp', '2c': 'kk', '2d': 'p', '2e': 'k', '2f': 'd'}
        written = ub.g2_2_integrated(c, k, Lam, vee)
        for key in w:
            ex = pint(lambda p: us.g2_2_pref(c, k, p) * w[key](p) * us.g2_2_groups62(c, k, p)[gname[key]])
            print('  seed %d G2-I(2) %s: written/exact = %s' % (seed, key, np.round(written[key] / ex, 10)))
        # line 62 vs 54 at one p
        p = 0.63
        print('  seed %d G2-I(2) line 62/line 54 (fixed colour) = %s' % (seed, np.round(
            us.g2_2_pref(c, k, p) * sum(v for v in [1]) * 0 + us.g2_2_line62(c, k, p) / fx.g2_2_colour58(c) * fx.g2_2_colour58(c) / fx.g2_2_line54(c, k, p), 10)))
        print('  seed %d G2-I(2) colour58/colour54 (fixed placements) = %s' % (seed, np.round(fx.g2_2_colour58(c) / fx.g2_2_colour(c), 10)))
        # G3 rows I, II
        ex1 = pint(lambda p: fx.g3_1_line44(c, k, p) / us.g3_pref(c, k, p) * us.g3_pref(c, k, p))
        L1 = np.sum((np.einsum('aB,CD,Cb,B,D->ab', c.U['xp'], c.U['y'], c.U['z'], c.rho['xp'], c.rho['y'])
                     - np.einsum('Db,DE,aC,C,E->ab', c.U['z'], c.U['y'], c.U['wp'], c.rho['xp'], c.rho['y'])) * fx.Rone_fix(c, 'w', 'z', 'x'))
        print('  seed %d G3-I line 48: written/exact = %s' % (seed, np.round(ub.g3_integrated(c, k, Lam, vee) * L1 / ex1, 10)))
