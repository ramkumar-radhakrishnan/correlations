import numpy as np
import fp
import user23 as us
import fixed23 as fx
from run1 import fp_rows, user_rows, k

FIX = {'G2-I(1) A': fx.g2_1_termA, 'G2-I(1) B': fx.g2_1_termB, 'G2-I(2)': fx.g2_2_line54, 'G3-I': fx.g3_1_line44,
       'G3-II': fx.g3_2_line55, 'G3-III(1)': fx.g3_3_part1_line68, 'G3-III(2)': fx.g3_3_part2_line80,
       'G3-V': fx.g3_5_line101, 'G3-VI': fx.g3_6_line120, 'G2-IV': fx.g2_4_line92}
for sym in (False, True):
    print('==== U symmetric' if sym else '==== general U')
    for seed in (1, 2, 3):
        c = fp.Cfg(seed, sym)
        p = 0.63
        a, b = fp_rows(c, p), user_rows(c, p)
        for key in a:
            fixd = FIX[key](c, k, p)
            print('  seed %d %-10s written/FP = %-22s corrected/FP = %s' % (seed, key, np.round(b[key] / a[key], 5), np.round(fixd / a[key], 10)))
