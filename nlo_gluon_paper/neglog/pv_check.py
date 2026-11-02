r"""Each corrected logarithm = principal value of the p+ integral it comes from (minus the log(vee/k+) you split off),
in the three regions A < aB, aB < A < B, A > B  (A = (x-z)^2, B = (x-w)^2, a = k+/vee).   -> pv_check.json"""
import json
import numpy as np
from scipy import integrate

K, VEE = 1.0, 7.0
a = K / VEE


def pv(A, B, lo, hi):
    """PV int_lo^hi dp A/(pA - kB)"""
    p0 = K * B / A
    if lo < p0 < hi:
        return integrate.quad(lambda p: 1.0, lo, hi, weight='cauchy', wvar=p0)[0]
    return integrate.quad(lambda p: A / (p * A - K * B), lo, hi)[0]


out = []
for B in (1.0,):
    for A in (0.05, 0.3, 0.9, 1.7, 4.0):          # A<aB, aB<A<B (x2), A>B (x2)
        reg = 'A<aB' if A < a * B else ('aB<A<B' if A < B else 'A>B')
        rows = {
            'part1 line2 / part2 lines 3-8: |(A-aB)/(A-B)|':
                (pv(A, B, K, VEE) - np.log(VEE / K), np.log(abs((A - a * B) / (A - B)))),
            'part1 line3 / part3a: |(B-A)/B|':
                (pv(A, B, 1e-12, K), np.log(abs((B - A) / B))),
            'part1 line4: |(A-aB)/B|':
                (pv(A, B, 1e-12, VEE) - np.log(VEE / K), np.log(abs((A - a * B) / B))),
        }
        for name, (num, cf) in rows.items():
            out.append({'region': reg, 'A': A, 'term': name, 'PV integral': num, 'log|.|': cf, 'diff': abs(num - cf)})
            print('%-7s A=%.2f  %-48s PV=% .12f  log|.|=% .12f' % (reg, A, name, num, cf))
assert max(o['diff'] for o in out) < 1e-8
json.dump(out, open('pv_check.json', 'w'), indent=1)
print('max diff', max(o['diff'] for o in out))
