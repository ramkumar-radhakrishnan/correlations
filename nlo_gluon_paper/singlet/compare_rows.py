"""Row by row: your two-rho evolution terms vs the region-T reference rows (LL wave function, index placement for a
general U), in the colour-singlet projection.

your row R = genuine terms (Sec. 1) + "from three rho's" (Sec. 2) + the two-rho part of the four-rho line of row R
             (Eq. (4rho), Weyl-symmetrized line by line; their sum must equal your "From 4 rho's" + the 4rho.tex lists).
Output: compare_rows.json"""
import json
import os
import sys
from collections import OrderedDict
import env  # noqa: F401
from core import Config, term, symmetrize
import sets as S2
from parse_user import sec2_terms
from orient import singlet_sum

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'two_rho'))
sys.path.insert(0, os.path.join(HERE, '..', 'three_rho'))
import reference as REF  # noqa: E402
import fourrho  # noqa: E402

EVO = os.path.join(HERE, '..', 'final_check', 'input', 'Evolution_terms.tex')
LINE_ROW = {'L1': 'G1-I', 'L2': 'G1-I', 'L3': 'G1-IV', 'L4': 'G1-III', 'L5': 'G2-I', 'L6': 'G2-II', 'L7': 'G2-III',
            'L8': 'G3-I', 'L9': 'G3-II', 'L10': 'G3-IV'}


def four_two_rho():
    """two-rho part of every four-rho word, per row"""
    out = OrderedDict()
    for c, row, lab, facs, word, kern, cc in fourrho.words():
        Us = [(f[1], f[2], f[3]) for f in facs if f[0] == 'U']
        k = [('KK', a, b, cc_, d) for (a, b), (cc_, d) in kern]
        t = term(c, Us, [], word, k, 'L')
        out.setdefault(LINE_ROW[lab.split('.')[0]], []).extend(symmetrize(t, (2,)))
    return out


def user_rows():
    A = S2.A_sets()
    sec2 = sec2_terms(EVO)
    by = {}
    for lab, t in sec2:
        by.setdefault(lab, []).append(t)
    g = lambda *labs: [t for l in labs for t in by.get(l, [])]
    F = four_two_rho()
    R = OrderedDict()
    R['G1-I'] = F['G1-I']
    R['G1-II'] = A['g:I.II'] + g('Group I (with c.c.)|Row II|')
    R['G1-III'] = g('Group I (with c.c.)|Row III|N1', 'Group I (with c.c.)|Row III|N2') + F['G1-III']
    R['G1-IV'] = g('Group I (with c.c.)|Row IV|N1', 'Group I (with c.c.)|Row IV|N2') + F['G1-IV']
    R['G1-V'] = A['g:I.V'] + g('Group I (with c.c.)|Row V|')
    R['G2-I'] = A['g:II.I'] + g('Group II|Row I|') + F['G2-I']
    R['G2-II'] = F['G2-II']
    R['G2-III'] = F['G2-III']
    R['G2-IV'] = A['g:II.IV']
    R['G3-I'] = g('Group III\\, (with c.c.)|Row I|N1', 'Group III\\, (with c.c.)|Row I|N2') + F['G3-I']
    R['G3-II'] = g('Group III\\, (with c.c.)|Row II|N1', 'Group III\\, (with c.c.)|Row II|N2') + F['G3-II']
    R['G3-III'] = A['g:III.III'] + g('Group III\\, (with c.c.)|Row III|')
    R['G3-IV'] = F['G3-IV']
    R['G3-V'] = g('Group III\\, (with c.c.)|Row V|')
    R['G3-VI'] = g('Group III\\, (with c.c.)|Row VI|')
    four_user = A['4rho'] + g("From $4\\, \\rho's$||")
    return R, F, four_user


def ref_rows():
    R = REF.rows()
    out = OrderedDict()
    out['G1-I'] = R['G1-I_Nafter'] + R['G1-I_Nbarbefore']
    out['G1-II'] = R['G1-II_A3'] + R['G1-II_A3bar']
    for k in ('G1-III', 'G1-IV', 'G1-V'):
        out[k] = R[k]
    out["G1-I'"] = R['G1-Ia']
    for k in ('G2-I', 'G2-II', 'G2-III', 'G2-IV', 'G3-I', 'G3-II', 'G3-III', 'G3-IV', 'G3-V', 'G3-VI'):
        out[k] = R[k]
    return out


def rel(a_terms, b_terms, sym, seeds=(3, 7, 11)):
    r = []
    for s in seeds:
        cfg = Config(s, sym)
        a, b = singlet_sum(a_terms, cfg), singlet_sum(b_terms, cfg)
        r.append(abs(a - b) / max(abs(b), abs(a), 1e-12))
    return max(r)


if __name__ == '__main__':
    U, F, four_user = user_rows()
    Rf = ref_rows()
    out = {}
    allF = [t for ts in F.values() for t in ts]
    out['four_rho_lines_vs_your_lists'] = dict(sym=rel(allF, four_user, True), general=rel(allF, four_user, False))
    print('two-rho part of Eq. (4rho), line by line, vs your "From 4 rho" + 4rho.tex lists:', out['four_rho_lines_vs_your_lists'])
    rows = OrderedDict()
    for k in Rf:
        u = U.get(k, [])
        rows[k] = dict(sym=rel(u, Rf[k], True), general=rel(u, Rf[k], False), n_user=len(u), n_ref=len(Rf[k]))
        print('%-8s your terms %4d  ref terms %6d   U=U^T: %.1e   general U: %.1e' % (k, len(u), len(Rf[k]), rows[k]['sym'], rows[k]['general']))
    out['rows'] = rows
    # G1-I + G1-I' together
    u = U['G1-I']
    r = Rf['G1-I'] + Rf["G1-I'"]
    out["G1-I + G1-I'"] = dict(sym=rel(u, r, True), general=rel(u, r, False))
    print("G1-I vs ref I(a)+I(b)+I':", out["G1-I + G1-I'"])
    json.dump(out, open(os.path.join(HERE, 'compare_rows.json'), 'w'), indent=1)
