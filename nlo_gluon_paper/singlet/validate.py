"""Check the trace form of every term: eval_words(singlet(t)) == c * N_c^ncp * colour(rho rho -> delta),
for general (non-symmetric) random adjoint Wilson lines.  Output: validate.json"""
import json
import os
import numpy as np
import env  # noqa: F401
from core import Config, facs, NC, placements
import sets as S2
import target as TG
from parse_user import sec2_terms
import trace as TR

HERE = os.path.dirname(os.path.abspath(__file__))
EVO = os.path.join(HERE, '..', 'final_check', 'input', 'Evolution_terms.tex')
GEN = ['g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III']


def direct(t, cfg, mp):
    fl = facs(t, mp)
    rh = [f for f in fl if f[0] == 'rho']
    fl2 = [f for f in fl if f[0] != 'rho'] + [('d', rh[0][1], rh[1][1])]
    return t['c'] * NC**t.get('ncp', 0) * cfg.colour(fl2)


def all_terms():
    A = S2.A_sets()
    out = [('sec2:' + l, t) for l, t in sec2_terms(EVO)]
    for k in GEN:
        out += [(k, t) for t in A[k]]
    out += [('4rho', t) for t in A['4rho']]
    return out


def check(terms, seeds=(1, 2)):
    worst = 0.0
    nzero = 0
    for seed in seeds:
        cfg = Config(seed, False)
        for lab, t in terms:
            mp = placements(t)[0][0]
            c, ws = TR.singlet(t)
            ws = [[(l[0], mp.get(l[1], l[1])) if l[0] != 'T' else l for l in w] for w in ws]
            v = TR.eval_words(c, ws, cfg.U, cfg.f) if c != 0 else 0.0
            d = direct(t, cfg, mp)
            if c == 0:
                nzero += 1
            worst = max(worst, abs(v - d) / (1 + abs(d)))
    return worst, nzero


if __name__ == '__main__':
    import refP
    L = all_terms()
    J = [('JIMWLK', t) for t in TG.jimwlk()]
    B = [('BFKL', t) for t in refP.bfkl()]
    out = {}
    for name, ts in (('lists', L), ('JIMWLK target', J), ('BFKL target', B)):
        w, nz = check(ts, seeds=(1, 2) if name != 'BFKL target' else (1,))
        out[name] = dict(terms=len(ts), max_abs_diff=w, zero_after_projection=nz)
        print('%-14s terms %5d   max |trace form - direct| = %.1e   (vanish after projection: %d)' % (name, len(ts), w, nz))
    json.dump(out, open(os.path.join(HERE, 'validate.json'), 'w'), indent=1)
