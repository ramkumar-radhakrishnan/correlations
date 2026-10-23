"""Reading rules for the index order of the lists (U^{bc} = U^{cb} in the user's files) and their test for general U."""
import os
import numpy as np
import env  # noqa: F401
from core import Config, NC
import sets as S2
import target as TG
from parse_user import sec2_terms
from chains import chains  # noqa: F401
import trace as TR
from validate import direct
from core import placements

EVO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'final_check', 'input', 'Evolution_terms.tex')
GEN = ['g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III', '4rho']


def walk(t, start_idx, frm):
    """U's along the chain from index start_idx (coming from item frm): list of (k, forward?) and the end item"""
    occ = {}
    for k, (i, j, p) in enumerate(t['Us']):
        occ.setdefault(i, []).append(('U', k, 0))
        occ.setdefault(j, []).append(('U', k, 1))
    for m, f in enumerate(t['fs']):
        for s, x in enumerate(f):
            occ.setdefault(x, []).append(('f', m, s))
    for r, (i, p) in enumerate(t['rs']):
        occ.setdefault(i, []).append(('r', r, 0))
    out, cur = [], start_idx
    while True:
        nxt = [o for o in occ[cur] if o != frm]
        o = nxt[0]
        if o[0] in ('r', 'f'):
            return out, o
        k, sl = o[1], o[2]
        i, j, p = t['Us'][k]
        out.append((k, sl == 0))
        cur, frm = (j, ('U', k, 1)) if sl == 0 else (i, ('U', k, 0))


def reorient(t, rule):
    Us = [list(u) for u in t['Us']]
    flips = set()
    for r, (idx, p) in enumerate(t['rs']):
        path, end = walk(t, idx, ('r', r, 0))
        if not path:
            continue
        # path from the charge outwards; "acts on the charge" = traversed backwards from the charge
        if rule in ('R1',):
            k, fwd = path[0]
            if fwd:
                flips.add(k)
        elif rule == 'R2':           # every U from the charge to the next f or to the middle of a rho-rho chain
            n = len(path) if end[0] == 'f' else (len(path) + 1) // 2
            for k, fwd in path[:n]:
                if fwd:
                    flips.add(k)
    for k in flips:
        i, j, p = Us[k]
        Us[k] = [j, i, p]
    return dict(t, Us=[tuple(u) for u in Us])


def singlet_sum(ts, cfg):
    tot = 0.0
    for t in ts:
        maps, wgt = placements(t)
        for mp in maps:
            from core import kernel_value, sub_kern
            g = lambda p, mp=mp: mp.get(p, p)
            tot += wgt * direct(t, cfg, mp) * kernel_value(sub_kern(t['kern'], g), cfg.P)
    return tot


if __name__ == '__main__':
    A = S2.A_sets()
    L = [t for _, t in sec2_terms(EVO)] + [t for k in GEN for t in A[k]]
    J = TG.jimwlk()
    out = {}
    for rule in ('as written', 'R1', 'R2'):
        Lr = L if rule == 'as written' else [reorient(t, rule) for t in L]
        res = []
        for sym in (True, False):
            for seed in (3, 7, 11):
                cfg = Config(seed, sym)
                a, j = singlet_sum(Lr, cfg), singlet_sum(J, cfg)
                res.append((sym, abs(a - j) / abs(j)))
        print('%-10s' % rule, ' '.join('%s:%.1e' % ('sym' if s else 'gen', r) for s, r in res))
        out[rule] = dict(sym=[r for s, r in res if s], general=[r for s, r in res if not s])
    import json
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'orient.json'), 'w'), indent=1)
