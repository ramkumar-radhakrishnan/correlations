r"""The BFKL term of Row IV = the coefficient of log(vee), i.e. d/dlog(vee) of the full Row IV
(your three lines, your p+ prescription), with a = k+/vee.  F(xi) = LO_xi[B], F(1) = LO[B].
Claim 1 (exact identity): d(L1+L2)/dlog vee = c F(a) [1/(1-a) + (1-a) + (1-a)/a^2]
          -- the line-2 log and the plus subtraction cancel term by term.
Claim 2: with F(a) = a^2 G(a), G(0) = (4 pi/k^2) LO^y, this is c G(0) + O(a) = H4 + O(a).
Claim 3: Row IV - log(vee/k+) c G(0) is vee-independent (up to O(a)); Row IV - log(vee/k+) c F(1) is not.
python3 dlogvee.py -> dlogvee.json"""
import json
import numpy as np
from scipy.integrate import quad


def rowIV(F, a):                                                                # c = 1; line 3 omitted (O(a^2) derivative)
    L1a = quad(lambda x: (x * (1 - x) + (1 - x) / x) * F(x) / x**2, a, 1, limit=400)[0]
    L1b = quad(lambda x: (F(x) - F(1.0)) / (x * (1 - x)), a, 1, limit=400)[0]   # your p+ prescription
    L2 = np.log(1 / a - 1) * F(1.0)                                              # log(vee/k+ - 1) of line 2
    return L1a + L1b + L2


def dlogvee(f, a, h=1e-4):                                                       # d/dlog vee = -d/dlog a
    return (f(a * np.exp(-h)) - f(a * np.exp(h))) / (2 * h)


out = {'identity': [], 'bfkl': []}
# Claim 1: arbitrary smooth F
F1 = lambda x: np.exp(-0.7 / x) * np.cos(1.3 * x) + 0.2 * x**2                  # noqa: E731
for a in (1e-2, 1e-3):
    num = dlogvee(lambda b: rowIV(F1, b), a)
    pred = F1(a) * (1 / (1 - a) + (1 - a) + (1 - a) / a**2)
    out['identity'].append({'a': a, 'numerical': num, 'formula': pred, 'rel': abs(num - pred) / abs(pred)})
    print('identity', out['identity'][-1])
a = 1e-3
print('line-2 part alone:', F1(1.0) / (1 - a), '  plus-subtraction part alone:', -F1(1.0) / (1 - a))
out['cancel'] = {'line2': F1(1.0) / (1 - a), 'plus_subtraction': -F1(1.0) / (1 - a)}

# Claims 2-3: F behaving like the real LO_xi[B]: F(xi) = xi^2 G(xi), G(0) = 1.7 (plays (4pi/k^2) LO^y), F(1) = 0.9
F2 = lambda x: x**2 * (1.7 + (0.9 - 1.7) * x)                                    # noqa: E731  G(0)=1.7, F(1)=0.9
for a in (1e-2, 1e-3, 1e-4):
    tot = dlogvee(lambda b: rowIV(F2, b), a)
    rem_H4 = dlogvee(lambda b: rowIV(F2, b) - np.log(1 / b) * 1.7, a)          # subtract log(vee/k+) H4
    rem_LO = dlogvee(lambda b: rowIV(F2, b) - np.log(1 / b) * F2(1.0), a)      # subtract log(vee/k+) LO[B]
    out['bfkl'].append({'a': a, 'd RowIV/dlog vee': tot, 'H4 (=G(0))': 1.7, 'LO[B] (=F(1))': F2(1.0),
                        'd(Row IV - log*H4)/dlog vee': rem_H4, 'd(Row IV - log*LO[B])/dlog vee': rem_LO})
    print('bfkl', out['bfkl'][-1])
json.dump(out, open('dlogvee.json', 'w'), indent=1)
