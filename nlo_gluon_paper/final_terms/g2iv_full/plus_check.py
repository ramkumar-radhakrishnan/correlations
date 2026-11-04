r"""Which plus prescription reproduces line 2 of the full Row IV, [log(k+/Lambda) + log(vee/k+ - 1)]?
p+ = k+(1-xi)/xi  (dxi/xi^2 * xi/(1-xi) = dp+/p+).   python3 plus_check.py -> plus_check.json"""
import json
import numpy as np
from scipy.integrate import quad

kp, Lam, vee = 1.0, 1e-6, 40.0
a = kp / vee
F = lambda x: np.cos(1.7 * x) * np.exp(-0.4 / x) + 0.3 * x                      # noqa: E731  smooth, F(1) = LO[B]
Fp = lambda p: F(kp / (p + kp))                                                  # noqa: E731  as a function of p+
out = {}
# the cut p+ >= Lambda (xi <= k+/(k+ + Lambda)) of  int dxi/xi^2 xi/(1-xi) F = int dp+/p+ F
cut = quad(lambda p: (Fp(p) - Fp(0.0)) / p, 0, vee - kp, limit=500)[0] + Fp(0.0) * np.log((vee - kp) / Lam)
pplus = quad(lambda p: (Fp(p) - Fp(0.0)) / p, 0, vee - kp, limit=500)[0]
out['cut - [p+ plus] = F(1)[log(k+/Lambda) + log(vee/k+ - 1)] (exact)'] = abs(cut - pplus - F(1.0) * (np.log(kp / Lam) + np.log(vee / kp - 1)))
# p+ plus vs the standard xi plus distribution of xi/(1-xi)_+ acting on h = F/xi^2, on [a,1]
std = quad(lambda x: (F(x) / x - F(1.0)) / (1 - x), a, 1, limit=500)[0] + F(1.0) * np.log(1 - a)
out['[p+ plus] = [standard xi plus] - F(1) log(vee/k+ - 1)'] = abs(pplus - (std - F(1.0) * np.log(vee / kp - 1)))
json.dump(out, open('plus_check.json', 'w'), indent=1)
print(out)
