r"""Your new 4rho.tex (input/), checked.

 1. Every equation against the previous, checked version (../llsplit/input/4rho.tex), numerically, pointwise
    (classical charges, random points and adjoint Wilson lines; log(k+/Lambda) and log(vee/k+) at generic values).
    Your new bracket [log(k+/Lambda) + log(vee/k+)] equals the old log(vee/Lambda).
 2. The four-rho terms of Eq. (4rho) cancel after symmetrization over the four charges (transcription
    ../three_rho/fourrho.py, checked against Eq. (4rho) in step 1 through the old file).
 3. The two-rho lists ('Two rho contribution', first and second part) against the checked lists of
    Final_evolution_check.pdf, Eq. (fourlists) (the earlier lists with F1, F2): F2 is always needed (an index occurs
    four times); the first part has lost its log factor in the new file.
Output: four_file.json"""
import json
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ('llsplit', 'bfkl_rows', 'three_rho', 'groups23'):
    p = os.path.normpath(os.path.join(HERE, '..', sub))
    if p not in sys.path:
        sys.path.append(p)
import texeval as TE  # noqa: E402
import fp  # noqa: E402

NEW = open(os.path.join(HERE, 'input', '4rho.tex')).read()
OLD = open(os.path.join(HERE, '..', 'llsplit', 'input', '4rho.tex')).read()
NAMES = TE.NAMES


def equations(src):
    """the align environments outside comments: (first line, text after the first '& =' or '&='"""
    lines = src.split('\n')
    out, i = [], 0
    while i < len(lines):
        if lines[i].lstrip().startswith(r'\begin{align}') and not lines[i].lstrip().startswith('%'):
            j = i
            while r'\end{align}' not in lines[j]:
                j += 1
            t = '\n'.join(lines[i:j + 1])
            t = t[t.index(r'\begin{align}') + len(r'\begin{align}'):t.rindex(r'\end{align}')]
            t = re.sub(r'\\label\{[^}]*\}', '', t)
            out.append((i + 1, t))
            i = j + 1
        else:
            i += 1
    return out


def rhs(t):
    m = re.search(r'&\s*=', t)
    return t[m.end():] if m else t


def normalize(t):
    """your notation -> the evaluator's: x^{\prime} -> x', ( - a + b) -> (b-a); the point u -> y' (y' does not occur)"""
    t = re.sub(r'\^\{\\prime\}', "'", t)
    t = re.sub(r'\+\s*\\(rm|mathrm|text)\{\s*c\.\s*c\.?\s*\}', '', t)       # "+ c.c." markers: same in both versions
    t = re.sub(r"\(\s*-\s*([a-z]'?)\s*\+\s*([a-z]'?)\s*\)", r"(\2-\1)", t)
    assert "y'" not in t
    return re.sub(r"(?<![A-Za-z\\])u(?![A-Za-z])", "y'", t)


def value(t, c, V=37.0, L=1e-3):
    t = TE.unwrap_colours(normalize(t), ('blue', 'red', 'green', 'orange', 'black'))
    return TE.term_value(t, TE.Env(c, KP=1.0, V=V, L=L))


def cfg(seed):
    return fp.Cfg(seed, sym=True, names=NAMES)


LOGB = r"{\color{blue}\left[\log{\left[\frac{k^{+}}{\Lambda}\right]}+\log{\left[\frac{\vee}{k^{+}}\right]}\right]}"
OLDLOG = r'\frac{\log{\left(\frac{\vee}{\Lambda}\right)}}{2\pi^{2}}'
OLDLOG2 = r'\log{\left(\frac{\vee}{\Lambda}\right)}'


def canon(t, new):
    """text without spaces and \, and without the log factor; and the number of log factors removed"""
    t = re.sub(r'\s+|\\,', '', t)
    if new:
        lb = re.sub(r'\s+|\\,', '', LOGB)
        n = t.count(lb)
        return t.replace(lb, ''), n
    o1, o2 = re.sub(r'\s+', '', OLDLOG), re.sub(r'\s+', '', OLDLOG2)
    n = t.count(o1)
    t = t.replace(o1, r'\frac{1}{2\pi^{2}}')
    n += t.count(o2)
    return t.replace(o2, ''), n


def compare_files():
    en, eo = equations(NEW), equations(OLD)
    assert len(en) == len(eo), (len(en), len(eo))
    out = []
    for (ln, tn), (lo, to) in zip(en, eo):
        r = []
        try:
            for seed in (1, 2):
                c = cfg(seed)
                a, b = value(rhs(tn), c), value(rhs(to), c)
                r.append(dict(new=complex(a), old=complex(b)))
        except TE.TexError as e:
            # e.g. "prefactor x rho rho rho rho [A] + [B]" meaning prefactor x (...)([A] + [B]): compare the text
            cn, nn = canon(rhs(tn), True)
            co, no = canon(rhs(to), False)
            out.append(dict(line_new=ln, line_old=lo, text_identical=(cn == co), logs_new=nn, logs_old=no,
                            not_evaluated=str(e)[:120]))
            print('new line %3d  old line %3d  text identical (without the log factor): %s  logs new/old %d/%d'
                  % (ln, lo, cn == co, nn, no))
            continue
        zero = all(max(abs(x['new']), abs(x['old'])) < 1e-12 for x in r)
        rel = 0.0 if zero else max(abs(x['new'] - x['old']) / max(abs(x['new']), abs(x['old']), 1e-300) for x in r)
        ratio = [x['new'] / x['old'] if abs(x['old']) > 1e-12 else 0.0 for x in r]
        out.append(dict(line_new=ln, line_old=lo, rel=float(rel), both_zero=zero,
                        size=float(max(abs(x['new']) for x in r)),
                        ratio=[[float(np.real(q)), float(np.imag(q))] for q in ratio],
                        log_vee_over_Lambda=float(np.log(37.0 / 1e-3))))
        print('new line %3d  old line %3d  rel %.1e  ratio %s  %s' % (ln, lo, rel, np.round(ratio, 6), 'both zero' if zero else ''))
    return out


def four_cancel():
    import four_check as FC
    return [dict(seed=s, residual=float(a), scale=float(b)) for s in (1, 2, 3) for a, b in [FC.check_four(s)]]


# ---------------------------------------------------------------- the two-rho lists
F1_START = r"&- 8 \frac{( - w' + x') \cdot ( - w + x') (y - z) \cdot (x' - z)}"
F2_OLD = r"f^{a b c} f^{d e b'} \rho^{b'}(x') \rho^{c}(y) U^{a b'}(y) U^{d b'}(w') U^{e b}(z)"
F2_NEW = r"f^{a b c} f^{d e g} \rho^{b'}(x') \rho^{c}(y) U^{a g}(y) U^{d b'}(w') U^{e b}(z)"


def y_to_x(t):
    """your lists have the charges at x' and y (and no point x): rename y -> x for the comparison"""
    t = t.replace(r"\int_{x,y,z,x',w,w'}", '')
    assert not re.search(r'(?<![A-Za-z\\])x(?![A-Za-z\'])', t), 'point x present'
    return re.sub(r'(?<![A-Za-z\\])y(?![A-Za-z])', 'x', t)


LVK = r'\log{\left(\frac{\vee}{k^{+}}\right)}'
NOLOG = r'\frac{g^{4}}{8\pi^{3}}\frac{1}{2\pi^{2}}'


def two_lists(variant):
    """the log(vee/k+) coefficient lists (log(vee/k+) set by the evaluator; the bracket [log(k+/Lambda) + log(vee/k+)]
    -> log(vee/k+)).  variant: 'as written' | 'log restored' (first part x [log(k+/Lambda)+log(vee/k+)] as in the
    earlier file) | 'log restored, F1, F2'"""
    eqs = {ln: rhs(t) for ln, t in equations(NEW)}
    first, second1, second2 = eqs[134], eqs[143], eqs[166]
    assert first.count(NOLOG) == 1 and LOGB not in first
    assert second1.count(F2_OLD) == 1
    second1 = second1.replace(F2_OLD, F2_NEW)                # F2
    if variant == 'F2':
        first = None                                          # no log: no log(vee/k+) term
    else:
        first = first.replace(NOLOG, r'\frac{g^{4}}{8\pi^{3}}\frac{%s}{2\pi^{2}}' % LVK)
    if variant == 'log restored, F1, F2':
        i = first.index(F1_START)
        j = first.index(r'\Bigg]', i)
        first = first[:i] + '& ' + first[j:]                 # F1: drop the -8 term, keep the closing bracket
    out = []
    for t in (first, second1, second2):
        if t is None:
            continue
        assert t.count(LOGB) <= 1
        out.append(y_to_x(normalize(t.replace(LOGB, LVK))))
    return out


def ref_fourlists():
    """Final_evolution_check.tex Eq. (fourlists): the two-rho lists of the earlier 4rho.tex with F1, F2, prefactor
    N log(vee/Lambda); here the log is log(vee/k+) (the evaluator's value), like the lists compared"""
    import compare_4rho as C4
    S = open(os.path.join(HERE, '..', 'final_check', 'Final_evolution_check.tex')).read()
    j = S.index(r'\label{eq:fourlists}')
    i = S.rindex(r'\begin{align}', 0, j)
    t = S[i + len(r'\begin{align}'):j]
    t = t[t.index('={}&') + 4:]
    assert r'\NNv' not in t and t.count(r'\NN') == 1
    t = t.replace(r'\NN', C4.PRE).replace(r'\tfrac', r'\frac').replace('N_c', 'N_{c}')
    return TE.unwrap_colours(t, ('CC', 'JIM', 'BFK'))


def compare_two():
    import compare_4rho as C4
    ref = [ref_fourlists()]
    out = {}
    for variant in ('F2', 'log restored, F2', 'log restored, F1, F2'):
        texts = two_lists(variant)
        res = []
        for seed in (1, 2, 3):
            c = TE.cfg(seed)
            a = sum(C4.sym(t, c) for t in texts)
            b = sum(C4.sym(t, c) for t in ref)
            res.append(float(abs(a - b) / max(abs(a), abs(b))))
        out[variant] = res
        print('two-rho lists %-22s vs Eq. (fourlists): %s' % (variant, res), flush=True)
    return out


if __name__ == '__main__':
    res = dict(equations=compare_files(), four_rho_cancel=four_cancel(), two_rho_lists=compare_two())
    print(res['four_rho_cancel'])
    json.dump(res, open(os.path.join(HERE, 'four_file.json'), 'w'), indent=1)
