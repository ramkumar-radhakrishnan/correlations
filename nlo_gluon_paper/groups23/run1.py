"""Pointwise comparison: written equations (before p+ integration) vs first principles."""
import numpy as np
import fp
import user23 as us

k, Lam, vee = 1.7, 0.11, 23.0


def fp_rows(c, p):
    P = p + k
    J = (P / k)**2
    out = {}
    out['G2-I(1) A'] = fp.PREF * fp.dot(fp.frakB_two(c, k, p, 'wp', 'z', 'xp', 'y'), fp.frakB_one(c, k, p, 'w', 'z', 'x'))
    out['G2-I(1) B'] = fp.PREF * fp.dot(fp.frakB_one(c, k, p, 'wp', 'z', 'x'), fp.frakB_two(c, k, p, 'w', 'z', 'xp', 'y'))
    out['G2-I(2)'] = fp.PREF * fp.dot(fp.frakB_one(c, k, p, 'wp', 'z', 'xp'), fp.frakB_one(c, k, p, 'w', 'z', 'x'))
    out['G3-I'] = fp.PREF / 2 * fp.dot(fp.X1p(c, k, p, 'wp', 'z', 'xp', 'y'), fp.frakB_one(c, k, p, 'w', 'z', 'x'))
    out['G3-II'] = fp.PREF / 2 * fp.dot(fp.X2p(c, k, p, 'wp', 'z', 'xp', 'y'), fp.frakB_one(c, k, p, 'w', 'z', 'x'))
    Cc = fp.Cp_strip(c, k, p, 'wp', 'z', 'yp', 'xp') * J
    out['G3-III(1)'] = -fp.PREF * fp.dot(Cc, fp.frakB_two(c, k, p, 'w', 'z', 'y', 'x'))
    out['G3-III(2)'] = -fp.PREF * fp.dot(Cc, fp.frakB_one(c, k, p, 'w', 'z', 'x'))
    out['G3-V'] = -fp.PREF / 2 * fp.dot(Cc, fp.X1p(c, k, p, 'w', 'z', 'y', 'x'))
    out['G3-VI'] = -fp.PREF / 2 * fp.dot(Cc, fp.X2p(c, k, p, 'w', 'z', 'y', 'x'))
    out['G2-IV'] = fp.PREF * fp.dot(fp.Cp_strip(c, k, p, 'wp', 'z', 'yp', 'xp'), fp.Cp_strip(c, k, p, 'w', 'z', 'y', 'x'))
    return out


def user_rows(c, p):
    return {'G2-I(1) A': us.g2_1_termA(c, k, p), 'G2-I(1) B': us.g2_1_termB(c, k, p),
            'G2-I(2)': us.g2_2_line54(c, k, p), 'G3-I': us.g3_1_line44(c, k, p), 'G3-II': us.g3_2_line55(c, k, p),
            'G3-III(1)': us.g3_3_part1_line68(c, k, p), 'G3-III(2)': us.g3_3_part2_line80(c, k, p),
            'G3-V': us.g3_5_line101(c, k, p), 'G3-VI': us.g3_6_line120(c, k, p), 'G2-IV': us.g2_4_line92(c, k, p)}


if __name__ == '__main__':
    for sym in (False, True):
        print('==== U symmetric' if sym else '==== general U')
        for seed in (1, 2):
            c = fp.Cfg(seed, sym)
            p = 0.63
            a, b = fp_rows(c, p), user_rows(c, p)
            for key in a:
                ratio = b[key] / a[key]
                print('  seed %d  %-10s  written/first-principles = %s' % (seed, key, np.round(ratio, 6)))
