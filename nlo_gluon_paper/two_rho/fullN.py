"""Total of all two-rho terms vs JIMWLK x LO, z-free colour sector integrated (dim. reg.)."""
import numpy as np
from collections import defaultdict
from core import config_pair, accumulate_N, accumulate_dr, term
import reference as RF
import target as TG
from compare import report
import rowcmp


def g1_fixed_constant(ts):
    """Row II genuine with the bracket constant replaced by that of int_z K (dim. reg.):
    BG1 -> ELL/pi, i.e. -pi Nc [P - log] -> -Nc ell"""
    out = []
    for t in ts:
        kern = [('ELL',) + f[1:] if f[0] == 'BG1' else f for f in t['kern']]
        out.append(dict(t, kern=kern, c=t['c'] / np.pi))
    return out


def build(S, sets, fixes, seed=7, symmetric=True, target='JIMWLK'):
    cfg, cfg2 = config_pair(seed, symmetric)
    acc = defaultdict(lambda: defaultdict(complex))
    for k in sets:
        ts = S[k]
        if k == 'G:G3_III' and 'G3sign' in fixes:
            ts = [dict(t, c=-t['c']) for t in ts]
        if k in ('G:G1_II', 'G:G2_IV'):
            if k == 'G:G1_II' and 'G1const' in fixes:
                ts = g1_fixed_constant(ts)
            accumulate_dr(ts, cfg, acc, label=k)
        else:
            accumulate_N(ts, cfg, cfg2, acc, label=k)
    if target == 'JIMWLK':
        accumulate_N([dict(t, c=-t['c']) for t in TG.jimwlk()], cfg, cfg2, acc, label='-JIMWLK')
    return acc


if __name__ == '__main__':
    S = rowcmp.user_sets()
    R = RF.rows()
    # G2-IV alone vs reference G2-IV
    cfg, cfg2 = config_pair(7)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate_dr(S['G:G2_IV'], cfg, acc, label='notes')
    accumulate_N([dict(t, c=-t['c']) for t in R['G2-IV']], cfg, cfg2, acc, label='-ref')
    bad = report(acc, show=False)
    print('G2-IV (dim. reg.): keys', len(acc), 'nonzero', len(bad))
    for key, tot, sc, d in bad:
        print('   ', key, {k: round(v.real, 4) for k, v in d.items()})
    base = [k for k in S if not k.startswith(('A1', 'A2', 'B1', 'B2', 'D', '-'))]
    corr = base + ['D', '-RowIII2', '-RowIIaU']
    for name, sets, fixes in [('as written', base, ()),
                              ('3rho corrections', corr, ()),
                              ('3rho corrections + G3 III sign', corr, ('G3sign',)),
                              ('3rho corrections + G3 III sign + Row II constant', corr, ('G3sign', 'G1const'))]:
        acc = build(S, sets, fixes)
        bad = report(acc, show=False)
        print('=====', name, ': keys', len(acc), 'nonzero', len(bad))
        for key, tot, sc, d in bad[:12]:
            print('  %9.4f %s' % (abs(tot), key))
