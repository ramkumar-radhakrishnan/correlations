r"""Every logarithm of Finite_terms_1.tex that contains vee, split exactly as

        log(...) = c * log(vee/k+) + R,

with R finite for vee -> infinity (vee enters R only through k+/vee).  c and R are given as LaTeX; this script
evaluates the LaTeX strings themselves (tex2py) and checks the identity at random points, for several vee, Lambda.
It also checks the two p+ integrals of Group III and the xi identity of Group II Row IV that the remaining terms use.
Output: split_logs.json"""
import json
import os
import re
import numpy as np
from scipy.integrate import quad
from scipy.special import exp1

HERE = os.path.dirname(os.path.abspath(__file__))
LVEE = r"\log{\left[\frac{\vee}{k^{+}}\right]}"


def lam(P, Q, order='PQ'):
    """log[((vee-k+) P + k+ Q)/(Lambda P + k+ Q)] -> remainder, keeping the order of the terms of the original"""
    if order == 'PQ':      # [(vee-k+) P + k+ Q] / [Lambda P + k+ Q]
        num = r"\left(1-\frac{k^{+}}{\vee}\right)%s+\frac{k^{+}}{\vee}\,%s" % (P, Q)
        den = r"\frac{\Lambda}{k^{+}}\,%s+%s" % (P, Q)
    elif order == 'QP':    # [(vee-k+) P + k+ Q] / [k+ Q + Lambda P]
        num = r"\left(1-\frac{k^{+}}{\vee}\right)%s+\frac{k^{+}}{\vee}\,%s" % (P, Q)
        den = r"%s+\frac{\Lambda}{k^{+}}\,%s" % (Q, P)
    else:                  # [k+ Q + (vee-k+) P] / [k+ Q + Lambda P]
        num = r"\frac{k^{+}}{\vee}\,%s+\left(1-\frac{k^{+}}{\vee}\right)%s" % (Q, P)
        den = r"%s+\frac{\Lambda}{k^{+}}\,%s" % (Q, P)
    return r"\log{\left[\frac{%s}{%s}\right]}" % (num, den)


A, B = r"(x-z)^{2}", r"(x-w)^{2}"
Ap, Bp = r"(x'-z)^{2}", r"(x'-w')^{2}"
Bw = r"(x-w')^{2}"
L1P = r"\left(-\log{\left[1+\frac{\Lambda}{k^{+}}\right]}\right)"
L1M = r"\log{\left[1-\frac{\Lambda}{k^{+}}\right]}"

# original (exact string of the file) -> (c, R)
SPLIT = {
    r"\log{\left[\frac{\left[\vee\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}{\left[(\Lambda + k^{+})\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}\right]}":
        (1, r"\log{\left[\frac{\left[(x-z)^{2} - \frac{k^{+}}{\vee}\, (x-w)^{2}\right]}{\left[\left(1+\frac{\Lambda}{k^{+}}\right) (x-z)^{2} - (x-w)^{2}\right]}\right]}"),
    r"\log{\left[\frac{\vee}{k^{+}-\Lambda}\right]}": (1, r"\left(-\log{\left[1-\frac{\Lambda}{k^{+}}\right]}\right)"),
    r"\log{\left[\frac{\left[k^{+}\,(x-w)^{2}-\vee\, (x-z)^{2}\right]}{\left[k^{+}\,(x-w)^{2}-\Lambda\, (x-z)^{2}\right]}\right]}":
        (1, r"\log{\left[\frac{\left[\frac{k^{+}}{\vee}\,(x-w)^{2}-(x-z)^{2}\right]}{\left[(x-w)^{2}-\frac{\Lambda}{k^{+}}\, (x-z)^{2}\right]}\right]}"),
    r"\log{\left[\frac{\vee}{\Lambda + k^{+}}\right]}": (1, L1P),
    r"\log{\left[\frac{\vee}{k^{+}+\Lambda}\right]}": (1, L1P),
    r"\log{\left[\frac{(x-w)^{2}\,k^{+}-(x-z)^{2}\,\vee}{(x-w)^{2}\,k^{+} - (x-z)^{2}\,(k^{+}+\Lambda)}\right]}":
        (1, r"\log{\left[\frac{(x-w)^{2}\,\frac{k^{+}}{\vee}-(x-z)^{2}}{(x-w)^{2} - (x-z)^{2}\,\left(1+\frac{\Lambda}{k^{+}}\right)}\right]}"),
    r"\log{\left[\frac{k^{+}-\Lambda}{\vee}\right]}": (-1, L1M),
    r"\log{\left(\frac{(k^{+}-\Lambda)}{\vee}\right)}": (-1, L1M),
    r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\,(x-w)^{2}}\right]}": (1, lam(A, B)),
    r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{\Lambda\, (x-z)^{2}+k^{+}\,(x-w)^{2}}\right]}": (1, lam(A, B)),
    r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w')^{2}}{\Lambda\, (x-z)^{2} + k^{+}\,(x-w')^{2}}\right]}": (1, lam(A, Bw)),
    r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2} + k^{+}\, (x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\, (x-w)^{2}}\right]}": (1, lam(A, B)),
    r"\log{\left[\frac{(\vee-k^{+})\,(x'-z)^{2} + k^{+}\, (x'-w')^{2}}{\Lambda\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}\right]}": (1, lam(Ap, Bp)),
    r"\log{\left[\frac{(\vee-k^{+})\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}{\Lambda\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}\right]}": (1, lam(Ap, Bp)),
    r"\log{\left[\frac{(\vee-k^{+})\, (x-z)^{2} + k^{+}\, (x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\, (x-w)^{2}}\right]}": (1, lam(A, B)),
    r"\log{\left[\frac{(\vee - k^{+})\, (x-z)^{2} + k^{+}\, (x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\, (x-w)^{2}}\right]}": (1, lam(A, B)),
    r"\log{\left[\frac{(\vee - k^{+})\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}{\Lambda\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}\right]}": (1, lam(Ap, Bp)),
    r"\log{\left[\frac{k^{+}\, (x'-w')^{2} + (\vee - k^{+})\, (x'-z)^{2}}{k^{+}\, (x'-w')^{2} + \Lambda\, (x'-z)^{2}}\right]}": (1, lam(Ap, Bp, 'QPn')),
    r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{k^{+}\,(x-w)^{2}+\Lambda\, (x-z)^{2}}\right]}": (1, lam(A, B, 'QP')),
    # Group II Row I 2a: vee cancels between numerator and denominator
    r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{k^{+}\, (x'-w')^{2} + (\vee-k^{+})\, (x'-z)^{2}}\right]}":
        (0, r"\log{\left[\frac{\left(1-\frac{k^{+}}{\vee}\right)(x-z)^{2}+\frac{k^{+}}{\vee}\,(x-w)^{2}}{\frac{k^{+}}{\vee}\, (x'-w')^{2} + \left(1-\frac{k^{+}}{\vee}\right) (x'-z)^{2}}\right]}"),
}
# context replacements (exact strings of the file)
CONTEXT = {
    # Group I Row II part 1
    r"\left[\frac{11}{3}+\log{\left[\frac{\vee}{k^{+}}\right]}+\log{\left[\frac{\vee}{k^{+}}-1\right]}\right]":
        (r"2\,\log{\left[\frac{\vee}{k^{+}}\right]}", r"\left[\frac{11}{3}+\log{\left[1-\frac{k^{+}}{\vee}\right]}\right]"),
    # Group I Row V 2b UV
    r"\left[-2\log\Big(\frac{\vee}{k^{+}}\Big)+2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]":
        (r"-2\log\Big(\frac{\vee}{k^{+}}\Big)", r"\left[2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]"),
}
# vee-dependent logs that stay as they are (no log(vee/k+): finite or vanishing for vee -> infinity)
KEEP = [
    r"\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}",
    r"\log{\left[1-\frac{k^{+}}{\vee}\right]}",
    r"\log\left(\frac{p^{2}}{k^{2}+\frac{k^{+}}{\vee}\,p^{2}}\right)",
    r"\log\left[\frac{k^{2}}{(k-p)^{2}}+\frac{k^{+}}{\vee}\left(1-\frac{k^{2}}{(k-p)^{2}}\right)\right]",
    r"\log\left(1-\frac{k^{+}}{\vee}\right)",
    r"\log{\left(1-\frac{k^{+}}{\vee}\right)}",
]


# ------------------------------------------------------------------ LaTeX -> python (for the identity check)
SQ = {('x', 'z'): 'A', ('x', 'w'): 'B', ("x'", 'z'): 'Ap', ("x'", "w'"): 'Bp', ('x', "w'"): 'Bw'}


def frac(t):
    while r'\frac{' in t:
        i = t.index(r'\frac{')
        j = i + 6
        d, k = 1, j
        while d:
            d += {'{': 1, '}': -1}.get(t[k], 0)
            k += 1
        num = t[j:k - 1]
        assert t[k] == '{'
        d, m = 1, k + 1
        while d:
            d += {'{': 1, '}': -1}.get(t[m], 0)
            m += 1
        den = t[k + 1:m - 1]
        t = t[:i] + '((' + num + ')/(' + den + '))' + t[m:]
    return t


def tex2py(t):
    t = t.replace(r'\,', ' ').replace(r'\!', '')
    t = re.sub(r"\((x'?)-(\w'?)\)\^\{2\}", lambda m: ' ' + SQ[(m.group(1), m.group(2))] + ' ', t)
    t = t.replace('k^{+}', ' k ').replace(r'\vee', ' V ').replace(r'\Lambda', ' L ')
    for a, b in ((r'\left[', '('), (r'\right]', ')'), (r'\left(', '('), (r'\right)', ')'), (r'\Big(', '('), (r'\Big)', ')'),
                 (r'\log^{2}', ' LOG2 '), (r'\log', ' LOG ')):
        t = t.replace(a, b)
    t = frac(t)
    t = t.replace('{', '(').replace('}', ')').replace('[', '(').replace(']', ')')
    # implicit products
    toks = re.findall(r"LOG2|LOG|[A-Za-z]+|\d+|[()+\-*/]", t)
    out = []
    for tk in toks:
        if out and (re.match(r'[A-Za-z0-9)]', out[-1][-1]) and not out[-1] in ('LOG', 'LOG2')) and \
                (re.match(r'[A-Za-z0-9(]', tk)):
            out.append('*')
        out.append(tk)
    return ''.join(out)


def ev(t, env):
    e = dict(env, LOG=lambda v: np.log(abs(v)), LOG2=lambda v: np.log(abs(v))**2)
    return eval(tex2py(t), {}, e)


def check_identities(n=200, seed=1):
    rng = np.random.default_rng(seed)
    worst = {}
    for orig, (c, R) in SPLIT.items():
        w = 0.0
        for _ in range(n):
            env = dict(A=rng.uniform(.1, 3), B=rng.uniform(.1, 3), Ap=rng.uniform(.1, 3), Bp=rng.uniform(.1, 3),
                       Bw=rng.uniform(.1, 3), k=rng.uniform(.5, 2))
            env['V'] = env['k'] * 10**rng.uniform(.1, 8)
            env['L'] = env['k'] * 10**rng.uniform(-6, -0.5)
            lhs = ev(orig, env)
            rhs = c * np.log(env['V'] / env['k']) + ev(R, env)
            w = max(w, abs(lhs - rhs) / max(1, abs(lhs)))
        worst[orig[:70]] = w
    for orig, (bf, R) in CONTEXT.items():
        w = 0.0
        for _ in range(n):
            env = dict(k=rng.uniform(.5, 2))
            env['V'] = env['k'] * 10**rng.uniform(.1, 8)
            env['L'] = env['k'] * 10**rng.uniform(-6, -0.5)
            lhs, b, r = ev(orig, env), ev(bf, env), ev(R, env)
            w = max(w, abs(lhs - b - r) / max(1, abs(lhs)))
        worst[orig[:70]] = w
    return worst


def check_remainders_finite(seed=2):
    """R has a limit for vee -> infinity, Lambda -> 0 (no hidden log(vee/k+) left)"""
    rng = np.random.default_rng(seed)
    out = {}
    env0 = dict(A=1.3, B=0.7, Ap=0.4, Bp=2.1, Bw=1.7, k=1.0, L=1e-12)
    for orig, (c, R) in SPLIT.items():
        vals = [ev(R, dict(env0, V=10.0**p)) for p in (6, 9, 12)]
        out[orig[:70]] = float(abs(vals[2] - vals[0]))
    return out


# ------------------------------------------------------------------ the p+ integrals (Group III) and xi identity (Group II Row IV)
def ein(z):
    return np.euler_gamma + np.log(z) + exp1(z)


def cquad(f, a, b):
    re_ = quad(lambda x: f(x).real, a, b, limit=2000)[0]
    im_ = quad(lambda x: f(x).imag, a, b, limit=2000)[0]
    return re_ + 1j * im_


def check_pplus():
    out = {}
    worst1 = worst2 = 0.0
    for kap, k, V, Aa, Bb in ((0.7, 1.0, 9.0, 1.3, 0.6), (-1.9, 0.8, 30.0, 0.5, 2.2), (2.4, 1.2, 6.0, 0.9, 0.9)):
        P = V - k
        # -int_0^P dp [1/p]_+ e^{-i kap p/k}  =  gamma + log(i kap P/k) + E1(i kap P/k)
        num = -cquad(lambda p: (np.exp(-1j * kap * p / k) - 1) / p if p > 0 else -1j * kap / k, 0, P)
        z = 1j * kap * P / k
        worst1 = max(worst1, abs(num - ein(z)) / abs(num))
        # int_0^P dp [f(p)-f(0)]/p,  f = e^{-i kap p/k}/(p A + k B)
        #   = -(1/(k B)) [gamma + log(i kap P/k) + E1(i kap P/k) + J(b)],  b = B/A,  J(b) = e^{i kap b}[E1(i kap b) - E1(i kap (b+P/k))]
        f = lambda p: np.exp(-1j * kap * p / k) / (p * Aa + k * Bb)
        num2 = cquad(lambda p: (f(p) - f(0)) / p if p > 0 else (-1j * kap / k / (k * Bb) - Aa / (k * Bb)**2), 0, P)
        b = Bb / Aa
        J = np.exp(1j * kap * b) * (exp1(1j * kap * b) - exp1(1j * kap * (b + P / k)))
        cl = -(1 / (k * Bb)) * (ein(z) + J)
        worst2 = max(worst2, abs(num2 - cl) / abs(num2))
    out['-int [1/p]_+ e  = gamma+log(i kappa (vee-k)/k)+E1'] = worst1
    out['int [1/p]_+ e/(pA+kB) = -(1/kB)[...+J(b)]'] = worst2
    # Group II Row IV: (e^{-i al/xi}-e^{-i al})/(xi(1-xi)) = (e^{-i al/xi}/xi - e^{-i al})/(1-xi) - e^{-i al}/xi
    w = 0.0
    for al in (0.3, 1.7, -2.2):
        for xi in (0.01, 0.3, 0.77, 0.999):
            l = (np.exp(-1j * al / xi) - np.exp(-1j * al)) / (xi * (1 - xi))
            r = (np.exp(-1j * al / xi) / xi - np.exp(-1j * al)) / (1 - xi) - np.exp(-1j * al) / xi
            w = max(w, abs(l - r) / abs(l))
    out['xi identity of Group II Row IV'] = w
    return out


if __name__ == '__main__':
    res = dict(identities=check_identities(), remainders_vee_to_inf=check_remainders_finite(), pplus=check_pplus())
    for k, v in res.items():
        print(k)
        for a, b in v.items():
            print('   %.1e  %s' % (b, a))
    json.dump(res, open(os.path.join(HERE, 'split_logs.json'), 'w'), indent=1)
