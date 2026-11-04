r"""The standard plus distribution on [a,1] is finite: check with an explicit cut 1-delta.
  int_a^{1-delta} dxi f/(1-xi) + f(1) log(delta)   ->   int_a^1 dxi (f-f(1))/(1-xi) + f(1) log(1-a)   as delta -> 0,
with f(xi) = F(xi)/xi.  The log(delta) of the cut integral and of the delta(1-xi) term cancel.
python3 plus_finite.py -> plus_finite.json"""
import json
import numpy as np
from scipy.integrate import quad

F = lambda x: np.exp(-0.7 / x) * np.cos(1.3 * x) + 0.2 * x**2       # noqa: E731
f = lambda x: F(x) / x                                               # noqa: E731
a = 1e-2
rhs = quad(lambda x: (f(x) - f(1.0)) / (1 - x), a, 1, limit=400)[0] + f(1.0) * np.log(1 - a)
out = []
for d in (1e-2, 1e-4, 1e-6, 1e-8):
    cut = quad(lambda x: f(x) / (1 - x), a, 1 - d, limit=800, points=[1 - 10 * d])[0]
    out.append({'delta': d, 'cut integral (diverges)': cut, 'cut + f(1) log delta': cut + f(1.0) * np.log(d), 'formula': rhs})
    print(out[-1])
json.dump(out, open('plus_finite.json', 'w'), indent=1)
