r"""Row-by-row checks of your BFKL.tex and Other_terms.tex.  Exact identities only (vee > k+, p+ > Lambda).

 A. partition:  BFKL.tex + Other_terms.tex  =  row - (its log(k+/Lambda) part), pointwise          [Lambda = 1e-9]
    Without the removed line (G1 V 2b UV, E1) and without your two lines of Group II Row IV, which come out of line 1
    of that row by exact identities (verified in exact_check.py); 'A_given' with all of BFKL.tex.
 B. BFKL:       the BFKL list (BFKL.tex without the removed line + every explicit log(vee/k+) of Other_terms.tex
                + the subtraction terms -f(0) log(vee/k+) and (B) of Group III Rows III, V, VI)
                =  the previously checked list (../llsplit: hidden-log rule, (B), H4), whose two-rho sum is BFKL x LO.
 C. remaining:  Other_terms.tex minus the extracted parts (Lambda -> 0)  =  the remainder of the exact split of the row.

Terms with a remaining p+, xi or momentum integral cannot be evaluated pointwise; they are compared as text (identical
up to spaces, colours and \,).  Output: check_rows.json"""
import json
import os
import re
import sys
import numpy as np
import rows as R
import extract as X
import gen as GN
import texeval as TE
import bterms as BT

sys.path.insert(0, R.LL)
HERE = R.HERE
HIDDEN_PREV = set(GN.HIDDEN_ROWS)
GN.HIDDEN_ROWS = set()             # exact split of every explicit logarithm
BLK, SPL = GN.all_splits()
GN.HIDDEN_ROWS = HIDDEN_PREV       # the previously checked split (for check B)
_, SPL_PREV = GN.all_splits()
GN.HIDDEN_ROWS = set()
UB, UO, EX = X.extract()
E = np.e
REMOVED = {('G1 V(2b UV)', 0)}                       # E1: Row V has no vee
FROM_LINE1 = {('G2 IV', 0), ('G2 IV', 1)}            # your -f(0) log(vee/k+) and H4: exact identities on line 1
G3_LINE1 = ('G3 III(1)', 'G3 III(2)', 'G3 V', 'G3 VI')


HB = BT.hidden_B_text()
H4 = UB['G2 IV'].terms[1]
ROWS = [n for n in SPL]


def val(B, t, env):
    if t is None:
        return 0.0
    return TE.term_value(B.prefix + ' (' + t + ') ' + B.suffix, env)


def free_val(t, env):
    return TE.term_value(t, env)


SCALE = [0.0]


def safe(f, *a):
    try:
        v = f(*a)
        SCALE[0] = max(SCALE[0], abs(v))
        return v, None
    except TE.TexError as e:
        if 'integral' in str(e):
            return 0.0, 'integral'
        raise


def norm(t):
    t = GN.strip_user_colours(t)
    return re.sub(r'\s+|\\,|\\notag|\\\\|&', '', t)


def sums(name, env):
    """the three numbers of checks A, B, C for one row, and the list of terms that could not be evaluated"""
    skipped = []
    # mine: row - J  (B + exact remainder), B list, F list
    mine_rest, mine_B, mine_F = 0.0, 0.0, 0.0
    for n, s in enumerate(SPL[name]):
        if s.get('power') or s.get('red'):
            continue
        v, e = safe(val, BLK[name], s.get('B'), env)
        mine_B += v
        w, e2 = safe(val, BLK[name], s.get('Fx'), env)
        if e or e2:
            skipped.append(('mine', n, norm(BLK[name].terms[n])))
        mine_rest += v + w
        f, _ = safe(val, BLK[name], X.lam0(s['F']) if s.get('F') else None, env) if not e2 else (0.0, None)
        mine_F += f
    # the previously checked list
    prev_B = 0.0
    for n, sp in enumerate(SPL_PREV[name]):
        if sp.get('power') or sp.get('red'):
            continue
        prev_B += safe(val, BLK[name], sp.get('B'), env)[0]
    if name in HB:
        prev_B += free_val(HB[name], env)
    if name == 'G2 IV':
        prev_B += val(UB['G2 IV'], H4, env)
    # yours
    removed, yours_B_final, yours_rest, yours_F = 0.0, 0.0, 0.0, 0.0
    if name in UB:
        for n, t in enumerate(UB[name].terms):
            v, e = safe(val, UB[name], t, env)
            yours_rest += v
            if (name, n) not in REMOVED:
                yours_B_final += v
            if (name, n) in REMOVED or (name, n) in FROM_LINE1:
                removed += v
    if name in G3_LINE1:          # from line 1: -f(0) log(vee/k+) (= minus the log(vee/k+) of line 2) and (B)
        yours_B_final += -safe(val, UO[name], EX[name][1]['B'], env)[0] + free_val(HB[name], env)
    if name in UO:
        for n, t in enumerate(UO[name].terms):
            v, e = safe(val, UO[name], t, env)
            if e:
                skipped.append(('yours', n, norm(t)))
            yours_rest += v
            s = EX[name][n]
            b, _ = safe(val, UO[name], s.get('B'), env)
            yours_B_final += b
            f, _ = safe(val, UO[name], s.get('F'), env) if not e else (0.0, None)
            yours_F += f
    return dict(A=(yours_rest - removed, mine_rest), A_given=(yours_rest, mine_rest),
                B=(yours_B_final, prev_B), C=(yours_F, mine_F)), skipped


def rel(a, b):
    if max(abs(a), abs(b), SCALE[0]) < 1e-10:          # both sides vanish (only terms of order Lambda = 1e-9)
        return 0.0
    return float(abs(a - b) / max(abs(a), abs(b), SCALE[0], 1e-30))


if __name__ == '__main__':
    out = {}
    for name in ROWS:
        r = dict(A=0.0, A_given=0.0, B=0.0, C=0.0)
        for seed in (1, 2):
            c = TE.cfg(seed)
            SCALE[0] = 0.0
            sA, sk = sums(name, TE.Env(c, KP=1.0, V=37.0, L=1e-9))
            sB, _ = sums(name, TE.Env(c, KP=1.0, V=E, L=1e-9))
            sC, _ = sums(name, TE.Env(c, KP=1.0, V=37.0, L=1e-9))
            r['A'] = max(r['A'], rel(*sA['A']))
            r['A_given'] = max(r['A_given'], rel(*sA['A_given']))
            r['B'] = max(r['B'], rel(*sB['B']))
            r['C'] = max(r['C'], rel(*sC['C']))
        # text comparison of the terms that cannot be evaluated
        mine_txt = sorted(t for w, n, t in sk if w == 'mine')
        your_txt = sorted(t for w, n, t in sk if w == 'yours')
        r['not_evaluable_identical'] = [a == b for a, b in zip(mine_txt, your_txt)] if len(mine_txt) == len(your_txt) else 'count differs %d/%d' % (len(mine_txt), len(your_txt))
        r['n_not_evaluable'] = len(mine_txt)
        out[name] = r
        print('%-14s A %.1e (with all of BFKL.tex %.1e)  B %.1e  C %.1e  non-evaluable identical: %s'
              % (name, r['A'], r['A_given'], r['B'], r['C'], r['not_evaluable_identical']))
    json.dump(out, open(os.path.join(HERE, 'check_rows.json'), 'w'), indent=1)
