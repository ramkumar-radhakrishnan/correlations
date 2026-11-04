r"""Your earlier BFKL.tex + Other_terms.tex (../bfkl_rows/input), whose exact BFKL bookkeeping was checked there, against
your NEW rows (input/Group_I/II/III.tex): the checks A, B, C of ../bfkl_rows/check_rows.py with the new rows.
 A: BFKL.tex + Other_terms.tex = new row - log(k+/Lambda) part   B: BFKL list = checked list   C: remainders equal.
So the BFKL list and the remaining terms written from BFKL.tex/Other_terms.tex are those of the new files.
Output: rows_vs_files.json"""
import json
import os
import ll  # noqa: F401  (new inputs; must come first)
import numpy as np
import texeval as TE
import check_rows as CR

HERE = ll.HERE

if __name__ == '__main__':
    out = {}
    for name in CR.ROWS:
        r = dict(A=0.0, A_given=0.0, B=0.0, C=0.0)
        for seed in (1, 2):
            c = TE.cfg(seed)
            CR.SCALE[0] = 0.0
            sA, sk = CR.sums(name, TE.Env(c, KP=1.0, V=37.0, L=1e-9))
            sB, _ = CR.sums(name, TE.Env(c, KP=1.0, V=np.e, L=1e-9))
            r['A'] = max(r['A'], CR.rel(*sA['A']))
            r['A_given'] = max(r['A_given'], CR.rel(*sA['A_given']))
            r['B'] = max(r['B'], CR.rel(*sB['B']))
            r['C'] = max(r['C'], CR.rel(*sA['C']))
        mine_txt = sorted(t for w, n, t in sk if w == 'mine')
        your_txt = sorted(t for w, n, t in sk if w == 'yours')
        r['not_evaluable_identical'] = all(a == b for a, b in zip(mine_txt, your_txt)) if len(mine_txt) == len(your_txt) else False
        r['n_not_evaluable'] = len(mine_txt)
        out[name] = r
        print('%-14s A %.1e (all of BFKL.tex %.1e)  B %.1e  C %.1e  non-evaluable identical: %s (%d)'
              % (name, r['A'], r['A_given'], r['B'], r['C'], r['not_evaluable_identical'], r['n_not_evaluable']), flush=True)
    json.dump(out, open(os.path.join(HERE, 'rows_vs_files.json'), 'w'), indent=1)
