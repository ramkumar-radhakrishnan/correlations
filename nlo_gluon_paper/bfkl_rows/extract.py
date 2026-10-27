r"""The BFKL terms hidden in Other_terms.tex, and the remaining terms.

For every term of Other_terms.tex the logarithm is split as in ../llsplit/logsplit.py:
    log X = cV log(vee/k+) + R,
cV log(vee/k+) goes to the BFKL list, R (Lambda -> 0) stays.  (cL = 0 everywhere: the log(k+/Lambda) terms are gone.)
Exceptions (hidden logarithms): in Group III Rows III, V, VI and Group II Row IV the log(vee/k+ - 1) of line 2 is kept,
because it cancels the log(vee/k+) of the subtraction in [1/p+]_+ (resp. xi/(1-xi)_+) of line 1.  Line 1 of these
rows contains the BFKL term (B) (resp. H4); it goes to the BFKL list and is subtracted from the remaining terms.

EX[name] = list of dict(B=..., F=..., Fx=...) for the terms of UO[name]."""
import re
import rows as R
import gen as GN
import logsplit as LS

LKL, LVK, LVK1 = LS.LKL, LS.LVK, LS.LVK1
HIDDEN_LINE = {('G2 IV', 1), ('G3 III(1)', 1), ('G3 III(2)', 1), ('G3 V', 1), ('G3 VI', 1)}
LAM0 = [(r'\int_{\Lambda}^{k^{+}-\Lambda}', r'\int_{\new{0}}^{\new{k^{+}}}'),
        (r'\int_{k^{+}-\Lambda}^{\Lambda}', r'\int_{\new{k^{+}}}^{\new{0}}'),
        (r'\left[\frac{1}{\vee} - \frac{1}{k^{+} + \Lambda}\right]', r'\left[\frac{1}{\vee} - \new{\frac{1}{k^{+}}}\right]'),
        (r'\int_{\frac{k^{+}}{\vee}}^{\frac{k^{+}}{k^{+}+\Lambda}}', r'\int_{\frac{k^{+}}{\vee}}^{\new{1}}'),
        (r'\left[2\log\Big(1-\frac{\Lambda}{k^{+}}\Big)-\frac{11}{6}\right]', r'\left[\new{-\frac{11}{6}}\right]'),
        (r'\mathcal{A}', r'\new{\mathcal{A}_{0}}')]


def lam0(t):
    for a, b in LAM0:
        t = t.replace(a, b)
    # small logarithms (no large part) with Lambda: their Lambda -> 0 form
    for lg in GN.BL.logs_in(t):
        k = LS.norm(lg)
        if k in LS.SPLIT and LS.SPLIT[k][:2] == (0, 0) and LS.SPLIT[k][2] is not None and r'\Lambda' in lg:
            t = t.replace(lg, r'\new{%s}' % LS.SPLIT[k][2])
    assert r'\Lambda' not in t.replace(r'\new{', ''), t[:300]
    return t


def split_other(name, n, t):
    t = GN.strip_user_colours(t)
    if (name, n) in HIDDEN_LINE:
        return dict(B=None, F=lam0(t), Fx=t, hidden=True)
    if name == 'G1 II(1)':
        # log(1 - k+/vee) log(vee/k+) -> 0 for vee -> infinity: no BFKL term
        return dict(B=None, F=lam0(t), Fx=t)
    u = GN.find_unit(t)
    if u is None:
        return dict(B=None, F=lam0(t), Fx=t)
    i, j, (cL, cV, Rm, Rx), txt = u
    assert cL == 0, (name, n, txt)
    if cV == 0:
        return dict(B=None, F=lam0(t), Fx=t)
    B = t[:i] + r'\new{%s}' % GN.coef(cV, LVK) + t[j:]
    F = lam0(t[:i] + r'\new{%s}' % Rm + t[j:]) if Rm is not None else None
    Fx = t[:i] + Rx + t[j:] if Rx is not None else None
    return dict(B=B, F=F, Fx=Fx, log=txt)


def extract():
    UB, UO = R.user_blocks()
    EX = {name: [split_other(name, n, t) for n, t in enumerate(B.terms)] for name, B in UO.items()}
    return UB, UO, EX


if __name__ == '__main__':
    UB, UO, EX = extract()
    for name, L in EX.items():
        print('==', name)
        for n, s in enumerate(L):
            print('   %2d B:%-5s F:%-5s %s' % (n, s['B'] is not None, s['F'] is not None,
                                           re.sub(r'\s+', '', s.get('log', ''))[:70] + (' (hidden)' if s.get('hidden') else '')))
