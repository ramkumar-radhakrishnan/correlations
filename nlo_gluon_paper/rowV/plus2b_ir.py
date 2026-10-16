"""Large-distance behaviour of the pieces of the plus prescription (x far from z, or z far from x).
Colour model: O(p) = a + b exp(-p^2) (U -> 1 far outside the target, so O does not vanish at large distance).
Prints r^2 x (angular average) at large r; a constant means the transverse integral is log divergent at large r."""
import json
import numpy as np
import plus2b_checks as C

P = C.Pt(1)
P.O = lambda p: complex(0.8 + 0.5 * np.exp(-C.dot(p, p)), 0.3)


def ang_far(fun, r, centre, n=2000):
    ts = (np.arange(n) + 0.5) * 2 * np.pi / n
    return r * r * np.mean([fun(centre + r * np.array([np.cos(t), np.sin(t)])) for t in ts])


xi0 = 0.37
rows = {
    'original integrand, xi = 0.37 (x far)': lambda x: C.I_full(C.with_x(P, x), xi0),
    'endpoint p+ -> 0, I0 (x far)': lambda x: C.I0_closed(P, x),
    'plus part, -(F4(xi)-F4(0))/xi at xi = 0.37 (x far)': lambda x: (lambda Q: -(Q_F4(Q, xi0) - Q_F4(Q, 0.0)) / xi0)(C.with_x(P, x)),
}


def Q_F4(Q, s):
    g = C.first(Q, s)
    return g['ph'] / g['D'] * g['B4'] * g['O']


OUT = {}
for name, f in rows.items():
    OUT[name] = [abs(ang_far(f, r, P.z)) for r in (10.0, 100.0, 1000.0)]
    print('%-55s r^2<.> at r = 10, 100, 1000: %s' % (name, ['%.3e' % v for v in OUT[name]]))


# endpoint p+ -> k+: z far from x
def I1_zfar(z):
    Q = C.with_x(P, P.x)
    Q.z = z
    return C.I1_closed(Q)


OUT['endpoint p+ -> k+, I1 (z far)'] = [abs(ang_far(I1_zfar, r, P.x)) for r in (10.0, 100.0, 1000.0)]
print('%-55s r^2<.> at r = 10, 100, 1000: %s' % ('endpoint p+ -> k+, I1 (z far)', ['%.3e' % v for v in OUT['endpoint p+ -> k+, I1 (z far)']]))
json.dump(OUT, open('plus2b_ir.json', 'w'), indent=1)
