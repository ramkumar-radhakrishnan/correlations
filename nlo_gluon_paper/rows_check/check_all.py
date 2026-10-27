r"""Every step of your Group_I/II/III.tex (input/), numerically, pointwise in the transverse coordinates.

Random points, random adjoint SU(3) Wilson lines (U = U^T, your convention), random charges (../groups23/fp.py).
Kinds of step:
  int   both sides evaluated in full; remaining p+, q+ or xi integrals done numerically (principal value where needed)
  plus  as int, for a step that uses your plus prescription and Lambda -> 0: the difference must be O(Lambda)
        (shown at Lambda = 1e-6 and 1e-7)
  pt    both sides have the same longitudinal integral: the integrands are compared at fixed p+ (or xi); a delta^2 of
        the left side is done by moving the point it fixes (with its Jacobian)
Usage: python3 check_all.py Group_I|Group_II|Group_III   -> check_<group>.json"""
import json
import sys
import numpy as np
import pe
import texeval as TE

KP, V = 1.0, 20.0
# Group II Eq. after the z integral: the phase of its second term, as it should read
PHASE_FIX = (r"e^{-ik\,\left(y'-y\right)\,\frac{p^{+}}{(p^{+}+k^{+})}}", r"e^{-ik\,\left(y'-y\right)\,\frac{p^{+}+k^{+}}{k^{+}}}")


def P(c, a):
    return c.pt[a]


# ------------------------------------------------------------------ delta functions of the explicit forms
def d_g3(c, p, k):          # Group III Rows III, V, VI: delta^2[w' - ((p+k) y' - p z)/k]
    return [('wp', P(c, 'yp') + (p / k) * (P(c, 'yp') - P(c, 'z')))], 1.0


def d_g1v(c, p, k):         # Group I Row V: delta^2[w - z + (p/k)(z - x)]
    return [('w', P(c, 'z') - (p / k) * (P(c, 'z') - P(c, 'x')))], 1.0


def d_g2iv(c, p, k):        # Group II Row IV: delta^2[y - w + p/(p+k) (w - z)] delta^2[y' - w' + p/(p+k) (w' - z)]
    r = (p + k) / k
    return [('w', P(c, 'z') + r * (P(c, 'y') - P(c, 'z'))), ('wp', P(c, 'z') + r * (P(c, 'yp') - P(c, 'z')))], r**4


def pv_g1iii(ctx):          # Group I Row III: 1/[p (x-z)^2 - k (x-w)^2] -> principal value
    c = ctx.cfg
    s = ctx.KP * np.sum((P(c, 'x') - P(c, 'w'))**2) / np.sum((P(c, 'x') - P(c, 'z'))**2)
    return lambda var: [s] if var == 'PP' else []


STEPS = {
    'Group_III': [
        ('Row I: p+ integral', 'int', [(43, 45)], [(47, 49)], {}),
        ('Row II: p+ integral', 'int', [(54, 56)], [(58, 60)], {}),
        ("Row III part 1: w' integral", 'pt', [(67, 69)], [(71, 73)], dict(delta=d_g3)),
        ('Row III part 1: plus prescription', 'plus', [(71, 73)], [(75, 77)], {}),
        ("Row III part 2: w' integral", 'pt', [(79, 82)], [(84, 87)], dict(delta=d_g3)),
        ('Row III part 2: plus prescription', 'plus', [(84, 87)], [(89, 93)], {}),
        ("Row V: w' integral", 'pt', [(100, 102)], [(104, 106)], dict(delta=d_g3)),
        ('Row V: plus prescription', 'plus', [(104, 106)], [(108, 110)], {}),
        ('Row V: Wilson lines', 'int', [(108, 110)], [(112, 114)], {}),
        ("Row VI: w' integral", 'pt', [(119, 121)], [(123, 125)], dict(delta=d_g3)),
        ('Row VI: plus prescription', 'plus', [(123, 125)], [(127, 129)], {}),
    ],
    'Group_II': [
        ('Row I part 1: p+ integral', 'int', [(45, 47)], [(49, 51)], {}),
        ('Row I part 2: Wilson lines', 'pt', [(53, 55)], [(57, 59)], {}),
        ('Row I part 2: simplifying', 'pt', [(57, 59)], [(61, 63)], {}),
        ('Row I part 2: p+ integral = 2a+...+2f', 'int', [(61, 63)], [(65, 67), (68, 70), (71, 73), (74, 76), (77, 79), (80, 82)], {}),
        ("Row IV: w and w' integrals", 'pt', [(91, 95)], [(97, 101)], dict(delta=d_g2iv)),
        ('Row IV: Wilson lines', 'pt', [(97, 101)], [(103, 105)], {}),
        ('Row IV: after the z integral -> simplifying', 'pt', [(111, 113)], [(115, 117)], {}),
        ('Row IV: the same, phase of the 2nd term of Eq. (111) corrected', 'pt', [(111, 113)], [(115, 117)], dict(fix=[PHASE_FIX])),
        ('Row IV: p+ -> xi', 'int', [(115, 117)], [(119, 121)], {}),
        ('Row IV: rewriting with P_gg', 'int', [(119, 121)], [(123, 125)], {}),
        ('Row IV: plus prescription in xi', 'plus', [(123, 125)], [(131, 133)], {}),
    ],
    'Group_I': [
        ('Row II part 1: rewriting', 'int', [(45, 47)], [(49, 51)], {}),
        ('Row III: p+ integrals = N1+N2+N3a+N3b', 'int', [(91, 93)], [(95, 97), (98, 100), (101, 120), (121, 129)], dict(pv=pv_g1iii, L=1e-2)),
        ('Row IV: p+ integral', 'int', [(138, 140)], [(142, 145)], {}),
        ('Row V part 1: w integral', 'pt', [(152, 154)], [(156, 158)], dict(delta=d_g1v, Ps=(0.13, 0.5, 0.81))),
        ('Row V part 1: plus prescription', 'plus', [(156, 158)], [(160, 162)], {}),
        ('Row V part 1: final form', 'int', [(160, 162)], [(244, 246)], {}),
        ('Row V part 2: w integral', 'pt', [(164, 166)], [(168, 170)], dict(delta=d_g1v, Ps=(0.13, 0.5, 0.81))),
        ('Row V part 2: rewriting', 'pt', [(168, 170)], [(172, 174)], dict(Ps=(0.13, 0.5, 0.81))),
        ('Row V part 2: identity f f U U = N_c U', 'pt', [(172, 174)], [(177, 179)], dict(Ps=(0.13, 0.5, 0.81))),
        ('Row V part 2: = 2a + 2b', 'pt', [(177, 179)], [(181, 183), (184, 186)], dict(Ps=(0.13, 0.5, 0.81))),
        ('Row V 2a: plus prescription', 'plus', [(181, 183)], [(188, 190)], {}),
        ('Row V 2a: final form', 'int', [(188, 190)], [(248, 250)], {}),
        ('Row V 2b finite: simplifying (q+ = k+ - p+)', 'int', [(205, 218)], [(220, 232)], dict(L=1e-2)),
        ('Row V 2b finite: plus prescription', 'plus', [(220, 232)], [(234, 241)], {}),
        ('Row V 2b finite: final form', 'int', [(234, 241)], [(256, 263)], {}),
    ],
}

# Row V 2b = 2b finite + counterterm (the term subtracted in 2b finite, added back as 2b UV after the x integral)
CT = r"""\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\int_{\Lambda}^{k^{+}-\Lambda}\frac{dp^{+}}{2\pi}\int_{x,y,z,x',w'}e^{-ik\cdot(w'-z)}\,
\frac{1}{k^{+}}\left[\frac{p^{+}(k^{+}-p^{+})}{k^{+2}}+\frac{p^{+}}{k^{+}-p^{+}}+\frac{k^{+}-p^{+}}{p^{+}}\right]
\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}\,(y-z)^{2}}\,\frac{(y-x)\cdot(y-z)}{(y-x)^{2}\,(x-z)^{2}}
\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[N_{c}\,U^{af}(y)-N_{c}\,U^{af}(z)\Big]\rho^{b'}(x')\,\rho^{f}(y)"""


def rel(a, b):
    a = complex(np.ravel(a)[0]) if np.ndim(a) else complex(a)
    b = complex(np.ravel(b)[0]) if np.ndim(b) else complex(b)
    return float(abs(a - b) / max(abs(a), abs(b), 1e-300))


# typos that stop the evaluation; corrected here and listed in the report
TYPOS = {'Group_II': [(r'\frac{(k^{+}+p^{+})}{k{+}}', r'\frac{(k^{+}+p^{+})}{k^{+}}'),
                      # your definition of the splitting function (text after Eq. "The pole is present ...")
                      (r'\mathcal{P}_{gg}(\xi)', r'2\,N_{c}\,\left[\xi\,(1-\xi)+\frac{(1-\xi)}{\xi}+\frac{\xi}{(1-\xi)}\right]')]}


def sides(group, eqs, fix=()):
    out = []
    for a, b in eqs:
        t = pe.align_rhs(group, a, b)
        for x, y in TYPOS.get(group, []) + list(fix):
            t = t.replace(x, y)
        t, n = pe.balance(t)
        if n:
            UNMATCHED.add((group, a, b))
        out.append(t)
    return out


UNMATCHED = set()          # equations with a closing \Bigg\} that has no opening one (display only)


def run_step(group, kind, e1, e2, opt, seed):
    c = TE.cfg(seed)
    L = opt.get('L', 1e-3)
    # d_perp = 2 - eps: the numerical index sums are 2-dimensional, so eps = 0 unless the step carries a 1/eps pole
    txt = ''.join(sides(group, e1, opt.get('fix', ())) + sides(group, e2))
    eps = opt.get('eps', 0.37 if r'\epsilon' in txt else 0.0)
    if kind in ('int', 'plus'):
        Ls = (1e-6, 1e-7) if kind == 'plus' else (L,)
        out = {}
        for Lv in Ls:
            ctx = pe.Ctx(c, KP, V, Lv, eps)
            pv = opt['pv'](ctx) if 'pv' in opt else None
            a = sum(pe.side_value(t, ctx, pv) for t in sides(group, e1, opt.get('fix', ())))
            b = sum(pe.side_value(t, ctx, pv) for t in sides(group, e2))
            out['rel_L%g' % Lv] = rel(a, b)
        return out
    ctx = pe.Ctx(c, KP, V, L, eps)
    Ps = np.array(opt.get('Ps', (0.3, 1.7, 6.0)))
    a = sum(pe.side_integrand(t, ctx, 'PP', Ps, opt.get('delta')) for t in sides(group, e1, opt.get('fix', ())))
    b = sum(pe.side_integrand(t, ctx, 'PP', Ps, opt.get('delta')) for t in sides(group, e2))
    return dict(rel=max(rel(x, y) for x, y in zip(np.atleast_1d(a), np.atleast_1d(b))))


def extra_group_I(seed):
    """Row V 2b: (2b) = (2b finite) + CT pointwise; CT's p+ coefficient vs the brackets of 2b UV"""
    c = TE.cfg(seed)
    ctx = pe.Ctx(c, KP, V, 1e-3)
    Ps = np.array((0.13, 0.5, 0.81))
    a = pe.side_integrand(pe.align_rhs('Group_I', 184, 186), ctx, 'PP', Ps)
    b = pe.side_integrand(pe.align_rhs('Group_I', 205, 218), ctx, 'PP', Ps) + pe.side_integrand(CT, ctx, 'PP', Ps)
    out = dict(two_b_equals_finite_plus_CT=max(rel(x, y) for x, y in zip(a, b)))
    # coefficient of the UV term: (1/k+) int_Lambda^{k+-Lambda} dp+ [...] vs your brackets
    for L in (1e-2, 1e-3):
        eps = L / KP
        xs, ws = pe.nodes(eps, 1 - eps)
        exact = float(np.sum(ws * (xs * (1 - xs) + xs / (1 - xs) + (1 - xs) / xs)))
        orange = 2 * np.log(V / L) - 2 * np.log(V / KP) + 2 * np.log(1 - L / KP) - 11 / 6 + 4 * eps - eps**2 + 2 * eps**3 / 3
        yours = np.log(KP / L) - np.log(V / KP) + 2 * np.log(1 - L / KP) - 11 / 6
        e1 = 2 * np.log(KP / L) + 2 * np.log(1 - L / KP) - 11 / 6
        out['UV_coefficient_L%g' % L] = dict(p_integral=exact, your_Eq_1_31=float(orange), your_final_bracket=float(yours),
                                             with_2log=float(e1), yours_minus_correct=float(yours - e1),
                                             minus_log_vee_over_Lambda=float(-np.log(V / L)))
    return out


if __name__ == '__main__':
    group = sys.argv[1]
    out = {}
    for name, kind, e1, e2, opt in STEPS[group]:
        res = [run_step(group, kind, e1, e2, opt, seed) for seed in (1, 2)]
        out[name] = res
        print('%-48s %s' % (name, res), flush=True)
    if group == 'Group_I':
        out['Row V 2b: UV term'] = [extra_group_I(s) for s in (1, 2)]
        print(out['Row V 2b: UV term'], flush=True)
    json.dump(out, open('check_%s.json' % group, 'w'), indent=1)
