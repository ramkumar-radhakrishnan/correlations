"""M1: Secs. 1 + 2 of the new file miss the two-rho lists of 4rho.tex.
 * size of Secs. 1+2 - JIMWLK x LO (U = U^T), pointwise and in the singlet projection;
 * the kernel monomials that differ;
 * the explicit list (Eq. fourlists of Final_evolution_check.tex) parsed back and compared with the checked 4rho.tex lists.
Output: check_m1.json, fourlists_body.tex"""
import json
import os
import re
from collections import defaultdict
import env  # noqa: F401
from core import Config, pointwise_total, accumulate
import sets as S2
import target as TG
from parse_v4 import equations, parse_body
from check_all import both

HERE = os.path.dirname(os.path.abspath(__file__))
FC = os.path.join(HERE, '..', 'final_check', 'Final_evolution_check.tex')

if __name__ == '__main__':
    E = equations()
    A = S2.A_sets()
    S12 = [t for e in E if e['sec'].startswith(('Genuine', 'Symmetrizing')) for t in e['terms']]
    J = TG.jimwlk()
    out = {}
    r = []
    for s in (3, 7, 11, 19):
        cfg = Config(s, True)
        a, _ = pointwise_total(S12, cfg)
        j, _ = pointwise_total(J, cfg)
        r.append(abs(a - j) / abs(j))
    out['|Secs1+2 - J| / |J|, U=U^T, pointwise'] = r
    # differing kernel monomials
    cfg = Config(5, True)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(S12, cfg, acc, False, label='S12')
    accumulate([dict(t, c=-t['c']) for t in J], cfg, acc, False, label='-J')
    bad = []
    for k, d in acc.items():
        s = sum(d.values())
        if abs(s) > 1e-9 * (1 + sum(abs(v) for v in d.values())):
            bad.append(str(k))
    out['kernel monomials total'] = len(acc)
    out['kernel monomials differing'] = bad
    # the explicit list of the earlier note
    tex = open(FC).read()
    m = re.search(r"\\frac\{d\^3N\}\{d\^3k\}\\bigg\|_\{4\\rho\\text\{-lists\}\}=\{\}&\\NN\\int_\{x,x',z,w,w'\}e\^\{-ik\\cdot\(w'-w\)\}\\Bigg\[(.*?)\\Bigg\]\\label\{eq:fourlists\}", tex, re.S)
    body = m.group(1)
    clean = re.sub(r"\{\\color\{\w+\}", "{", body).replace(r'\tfrac', r'\frac').replace(r'\;', ' ').replace(r'\quad', '')
    L = parse_body(clean, 'fourlists')
    out['explicit list: terms'] = len(L)
    out['explicit list vs checked 4rho.tex lists'] = both(L, A['4rho'])
    out['Secs1+2 + explicit list vs J'] = both(S12 + L, J)
    open(os.path.join(HERE, 'fourlists_body.tex'), 'w').write(body)
    for k, v in out.items():
        print(k, ':', v)
    json.dump(out, open(os.path.join(HERE, 'check_m1.json'), 'w'), indent=1)
