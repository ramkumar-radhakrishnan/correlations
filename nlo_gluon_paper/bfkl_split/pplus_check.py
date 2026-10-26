"""Is the [1/p+]_+ term of Group III Rows III, V, VI finite for vee -> infinity?

    I(vee) = int_0^{vee-k+} dp+ [1/p+]_+ e^{-i kappa p+/k+} = int_0^{vee-k+} dp+ (e^{-i kappa p+/k+} - 1)/p+,
    kappa = k.(y'-z).

Numerical p+ integral (scipy quad with oscillatory weights, no closed form used) for growing vee/k+, compared with
-log(vee/k+).  Split at p+ = k+:  [0, k+] is vee-independent;  e^{-i kappa p/k}/p on [k+, vee-k+] converges;
-1/p on [k+, vee-k+] gives -log(vee/k+ - 1).  Same for part 2 with f(p) = e^{-i kappa p/k}/(p A + k B).
Group I Row V has [k+/p+]_+ on 0 < p+ < k+ (no vee): nothing grows there.  Output: pplus_check.json"""
import json
import numpy as np
from scipy.integrate import quad


def I_num(kap, k, V):
    P = V - k
    w = kap / k
    # (e^{-iwp}-1)/p = (cos(wp)-1)/p - i sin(wp)/p
    lo_re = quad(lambda p: (np.cos(w * p) - 1) / p if p > 0 else 0.0, 0, k, limit=200)[0]
    lo_im = -quad(lambda p: np.sin(w * p) / p if p > 0 else w, 0, k, limit=200)[0]
    hi_re = quad(lambda p: 1 / p, k, P, weight='cos', wvar=w, limit=2000)[0] - np.log(P / k)
    hi_im = -quad(lambda p: 1 / p, k, P, weight='sin', wvar=w, limit=2000)[0]
    return complex(lo_re + hi_re, lo_im + hi_im)


def I2_num(kap, k, V, A, B):
    """int_0^{vee-k} dp [f(p)-f(0)]/p,  f = e^{-i kap p/k}/(p A + k B);  times k B"""
    P = V - k
    w = kap / k
    g = lambda p: 1 / (p * A + k * B)
    re_ = quad(lambda p: (np.cos(w * p) * g(p) - g(0)) / p if p > 0 else -A * g(0)**2, 0, k, limit=200)[0]
    im_ = -quad(lambda p: np.sin(w * p) * g(p) / p if p > 0 else w * g(0), 0, k, limit=200)[0]
    re_ += quad(lambda p: g(p) / p, k, P, weight='cos', wvar=w, limit=2000)[0] - g(0) * np.log(P / k)
    im_ += -quad(lambda p: g(p) / p, k, P, weight='sin', wvar=w, limit=2000)[0]
    return complex(re_, im_) * k * B


if __name__ == '__main__':
    out = {}
    k = 1.0
    for kap in (0.7, -2.3):
        rows = []
        for r in (1e1, 1e2, 1e3, 1e4, 1e5, 1e6):
            I = I_num(kap, k, r * k)
            rows.append(dict(vee_over_k=r, I=[I.real, I.imag], I_plus_log=[(I + np.log(r)).real, (I + np.log(r)).imag]))
        out['part1 kappa=%g' % kap] = dict(rows=rows, limit_minus_gamma_minus_log_ikappa=[-np.euler_gamma - np.log(abs(kap)), -np.pi / 2 * np.sign(kap)])
    rows = []
    for r in (1e1, 1e2, 1e3, 1e4, 1e5, 1e6):
        I = I2_num(0.7, k, r * k, 1.3, 0.6)
        rows.append(dict(vee_over_k=r, kB_I=[I.real, I.imag], kB_I_plus_log=[(I + np.log(r)).real, (I + np.log(r)).imag]))
    out['part2 kappa=0.7 A=1.3 B=0.6'] = dict(rows=rows)
    for name, d in out.items():
        print(name)
        for row in d['rows']:
            a = row.get('I', row.get('kB_I'))
            b = row.get('I_plus_log', row.get('kB_I_plus_log'))
            print('  vee/k+ = %8.0e   I = %9.4f %+9.4fi    I + log(vee/k+) = %8.5f %+8.5fi' % (row['vee_over_k'], a[0], a[1], b[0], b[1]))
        if 'limit_minus_gamma_minus_log_ikappa' in d:
            print('  -gamma_E - log(i kappa) = %8.5f %+8.5fi' % tuple(d['limit_minus_gamma_minus_log_ikappa']))
    json.dump(out, open('pplus_check.json', 'w'), indent=1)


# ---------------------------------------------------------------- keeping the term unevaluated
# [f(p)-f(0)]/p = [f(p) - f(0) k/(p+k)]/p - f(0)/(p+k),   int_0^{vee-k} dp/(p+k) = log(vee/k)  exactly.
# So  int_0^{vee-k} dp [1/p]_+ f = int_0^{vee-k} dp [1/p]_+^{(k)} f - f(0) log(vee/k),
#     int dp [1/p]_+^{(k)} f  ==  int dp [f(p) - f(0) k/(p+k)]/p   (no log(vee/k): finite for vee -> infinity).
def Itil_num(kap, k, V, g=lambda p: 1.0):
    """int_0^{V-k} dp [f(p) - f(0) k/(p+k)]/p,  f = e^{-i kap p/k} g(p), numerically"""
    P = V - k
    w = kap / k
    lo_re = quad(lambda p: (np.cos(w * p) * g(p) - g(0) * k / (p + k)) / p if p > 0 else 0.0, 0, k, limit=400)[0]
    lo_im = -quad(lambda p: np.sin(w * p) * g(p) / p if p > 0 else w * g(0), 0, k, limit=400)[0]
    hi_re = (quad(lambda p: g(p) / p, k, P, weight='cos', wvar=w, limit=4000)[0]
             - g(0) * np.log(2 * P / (P + k)))          # int_k^P dp k/(p(p+k)), elementary
    hi_im = -quad(lambda p: g(p) / p, k, P, weight='sin', wvar=w, limit=4000)[0]
    return complex(lo_re + hi_re, lo_im + hi_im)


if __name__ == '__main__':
    res = {}
    k = 1.0
    A, B = 1.3, 0.6
    g2 = lambda p: 1 / (p * A + k * B)
    for name, kap, g, I_old in (('part1 kappa=0.7', 0.7, lambda p: 1.0, lambda V: I_num(0.7, k, V)),
                                ('part2 kappa=0.7 A=1.3 B=0.6', 0.7, g2, lambda V: I2_num(0.7, k, V, A, B) / (k * B))):
        rows = []
        for r in (1e1, 1e2, 1e3, 1e4, 1e5, 1e6):
            It = Itil_num(kap, k, r * k, g)
            Io = I_old(r * k)
            ident = abs(Io - (It - g(0) * np.log(r))) / abs(Io)
            rows.append(dict(vee_over_k=r, Itilde=[It.real, It.imag], identity_rel=ident))
            print('%-28s vee/k+ = %7.0e   I~ = %9.5f %+9.5fi   |I - (I~ - f(0) log(vee/k+))|/|I| = %.1e'
                  % (name, r, It.real, It.imag, ident))
        res[name] = rows
    d = json.load(open('pplus_check.json'))
    d['kept unevaluated: I~ = int [1/p]_+^(k) f'] = res
    json.dump(d, open('pplus_check.json', 'w'), indent=1)
