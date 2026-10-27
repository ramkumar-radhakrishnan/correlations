r"""The three-rho terms of the BFKL list cancel?  The test of ../llsplit/three_check.py (pointwise, charges symmetrized,
c.c. for Groups I and III, log(vee/k+) = 1), applied to
  final:     your BFKL.tex (minus the two removed two-rho terms) + the log(vee/k+) parts extracted from Other_terms.tex
             + the terms (B), + your four-rho remainders d^3N_a, d^3N_b;
  as given:  your BFKL.tex alone (with d^3N_a, d^3N_b), i.e. without the extracted parts;
  final without (B): the final list without the terms (B) (their three-rho parts cancel among themselves).
Output: three_bfkl.json"""
import json
import os
import re
import numpy as np
import rows as R
import extract as X
import check_rows as CR
import three_check as T3
import texeval as TE

HERE = R.HERE
UB, UO, EX = CR.UB, CR.UO, CR.EX


def rhs(name):
    a, b = R.B_RANGES[name]
    t = R.align_text('BFKL', a, b)
    t = re.sub(r"\{\\color\{(red|orange)\}\s*\\alpha'?\}", '', t)     # your labels alpha, alpha'
    return t[re.search(r'&\s*=', t).end():]


def full(B, t):
    return B.prefix + ' (' + t + ') ' + B.suffix


def term_lists():
    final, given = [], []
    for name, B in UB.items():
        if name.startswith('4rho'):
            continue
        for n, t in enumerate(B.terms):
            if (name, n) in CR.REMOVED:
                continue
            tt = T3.z_form(full(B, t)) if name == 'G1 II(2a)' else full(B, t)
            item = ('BFKL.tex %s #%d' % (name, n + 1), tt, not name.startswith('G2'))
            final.append(item)
            given.append(item)
    for name, L in EX.items():
        for n, s in enumerate(L):
            if s.get('B') is not None:
                final.append(('extracted %s #%d' % (name, n + 1), full(UO[name], s['B']), not name.startswith('G2')))
    for name, t in CR.HB.items():
        final.append(('(B) %s' % name, t.replace("y'", "w'"), True))
    four = [('BFKL.tex d3N_a', rhs('4rho a'), False), ('BFKL.tex d3N_b', rhs('4rho b'), False)]
    return final + four, given + four


def run(items, seeds=(11, 12, 13)):
    c0 = T3.configs(seeds[0])
    three = [(lab, t, cc) for lab, t, cc in items if T3.degree(t, c0) == 3]
    out = []
    table = {}
    for seed in seeds:
        c = T3.configs(seed)
        unit = TE.fp.G**4 / ((2 * np.pi)**3 * 16 * np.pi**5) * np.exp(-1j * (c.kv @ (c.pt['wp'] - c.pt['w'])))
        tot, big = 0.0, 0.0
        for lab, t, cc in three:
            v = T3.sym_value(t, c, cc) / unit
            tot += v
            big = max(big, abs(v))
            if seed == seeds[0]:
                table[lab] = float(v.imag)
        out.append(dict(seed=seed, residual=float(abs(tot)), largest=float(big), relative=float(abs(tot) / big)))
    return dict(terms=[lab for lab, _, _ in three], per_seed=out, worst=max(r['relative'] for r in out), seed11=table)


if __name__ == '__main__':
    final, given = term_lists()
    res = {}
    no_b = [it for it in final if not it[0].startswith('(B)')]     # the three-rho parts of (B) cancel among themselves
    for key, items in (('final', final), ('as given', given), ('final without (B)', no_b)):
        res[key] = run(items)
        print('%-9s %d three-rho terms, worst relative residual %.2e' % (key, len(res[key]['terms']), res[key]['worst']))
    json.dump(res, open(os.path.join(HERE, 'three_bfkl.json'), 'w'), indent=1)
