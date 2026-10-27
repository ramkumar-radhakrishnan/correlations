"""Direct demonstration: the -(1/k+) term of Group III Rows III, V, VI contains log(vee/k+).
p+ integral first, then the y' integral near z: (2pi/i) int_0^inf dr/r [J0(|k| r) - J0(lam |k| r)] G(r), lam = vee/k+.
Output: direct_check.json"""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import j0

out = {}
for lam in (10, 100, 1000, 10000):
    f = lambda r: (j0(r) - j0(lam * r)) * np.exp(-r * r / 2) / r
    pts = np.concatenate([np.geomspace(1e-8, 1 / lam, 5), np.linspace(1 / lam, 40, 4000)])
    tot = sum(quad(f, a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))
    out[str(lam)] = dict(integral=tot, log=float(np.log(lam)), difference=tot - float(np.log(lam)))
    print(lam, out[str(lam)])
json.dump(out, open('direct_check.json', 'w'), indent=1)
