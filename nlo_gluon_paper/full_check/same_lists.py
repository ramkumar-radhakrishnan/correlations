r"""The JIMWLK (log k+/Lambda) and BFKL (log vee/k+) parts of your new rows equal, row by row and pointwise, those of the
files checked before (../llsplit/input), whose two-rho sums are JIMWLK x LO and BFKL x LO (Final_evolution_check.pdf).
Also the remainder F.  Output: same_lists.json"""
import json
import os
import numpy as np
import ll
import blocks as BL
import gen as GN
import texeval as TE

HERE = ll.HERE


def splits(files, final):
    BL.FILES, BL.FINAL = files, final
    return GN.all_splits()


def row_value(B, S, part, env):
    tot = 0.0
    for n, s in enumerate(S):
        t = s.get(part)
        if t is None or s.get('power') or s.get('red'):
            continue
        try:
            tot += TE.term_value(B.prefix + ' (' + t + ') ' + B.suffix, env)
        except TE.TexError as e:
            if 'integral' in str(e):
                continue
            raise
    return tot


if __name__ == '__main__':
    new_files, new_final = dict(BL.FILES), dict(BL.FINAL)
    old_files = {g: open(os.path.join(ll.LL, 'input', '%s.tex' % g)).read().split('\n')
                 for g in ('Group_I', 'Group_II', 'Group_III', '4rho')}
    import importlib
    old_final = dict(importlib.reload(importlib.import_module('blocks')).FINAL)
    BLKn, SPLn = splits(new_files, new_final)
    BLKo, SPLo = splits(old_files, old_final)
    out = {}
    for name in BLKn:
        r = {}
        for part in ('J', 'B', 'F'):
            worst = 0.0
            for seed in (1, 2):
                env = TE.Env(TE.cfg(seed), KP=1.0, V=37.0, L=1e-9)
                a, b = row_value(BLKn[name], SPLn[name], part, env), row_value(BLKo[name], SPLo[name], part, env)
                if max(abs(a), abs(b)) > 0:
                    worst = max(worst, float(abs(a - b) / max(abs(a), abs(b))))
            r[part] = worst
        out[name] = r
        print('%-14s J %.1e  B %.1e  F %.1e' % (name, r['J'], r['B'], r['F']), flush=True)
    json.dump(out, open(os.path.join(HERE, 'same_lists.json'), 'w'), indent=1)
