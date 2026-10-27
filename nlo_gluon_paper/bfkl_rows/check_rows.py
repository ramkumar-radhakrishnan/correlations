r"""Row-by-row checks of your BFKL.tex and Other_terms.tex against the exact split of your rows (../llsplit):

 A. partition:  BFKL.tex + Other_terms.tex  =  row - (its log(k+/Lambda) part)          [Lambda = 1e-9]
 B. BFKL:       BFKL.tex (minus the two removed terms) + the log(vee/k+) parts extracted from Other_terms.tex
                + (B) / H4  =  the log(vee/k+) list of LL_terms.pdf                         [log(vee/k+) = 1]
 C. remaining:  Other_terms.tex minus the extracted parts (Lambda -> 0) minus (B) / H4  =  Finite_terms_rows.pdf

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

sys.path.insert(0, R.LL)
HERE = R.HERE
BLK, SPL = GN.all_splits()
UB, UO, EX = X.extract()
E = np.e
REMOVED = {('G1 V(2b UV)', 0), ('G2 IV', 0)}          # BFKL.tex terms that are not BFKL terms (see the note)


def hidden_B():
    S = open(os.path.join(HERE, '..', 'bfkl_final', 'static.tex')).read()
    i = S.index(r'\paragraph{The new terms}')
    blk = S[S.index(r'\begin{align}', i) + len(r'\begin{align}'):S.index(r'\end{align}', i)]
    parts = re.split(r'\\text\{Row ([IV()12]+):\}\\quad&', blk)[1:]
    out = {}
    for k in range(0, len(parts), 2):
        t = parts[k + 1].replace(r'\Lv', GN.LVK).replace(r'\nc', '')
        out[parts[k]] = re.sub(r'[,.]\s*(\\\\)?\s*$', '', t.strip())
    return {'G3 III(1)': out['III(1)'], 'G3 III(2)': out['III(2)'], 'G3 V': out['V'], 'G3 VI': out['VI']}


HB = hidden_B()
H4 = UB['G2 IV'].terms[1]          # your H4 line (BFKL.tex, Group II Row IV, line 2)
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
    # yours
    removed, yours_B_final, yours_rest, yours_F = 0.0, 0.0, 0.0, 0.0
    if name in UB:
        for n, t in enumerate(UB[name].terms):
            v, e = safe(val, UB[name], t, env)
            yours_rest += v
            if (name, n) not in REMOVED:
                yours_B_final += v
            else:
                removed += v
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
    # (B) and H4
    if name in HB:
        hb = free_val(HB[name], env)
        mine_B += hb
        yours_B_final += hb
        mine_F -= hb
        yours_F -= hb
    if name == 'G2 IV':
        h4 = val(UB['G2 IV'], H4, env)
        mine_B += h4
        mine_F -= h4
        yours_F -= h4
    return dict(A=(yours_rest, mine_rest), B=(yours_B_final, mine_B, yours_B_final + removed), C=(yours_F, mine_F)), skipped


def rel(a, b):
    if max(abs(a), abs(b), SCALE[0]) < 1e-10:          # both sides vanish (only terms of order Lambda = 1e-9)
        return 0.0
    return float(abs(a - b) / max(abs(a), abs(b), SCALE[0], 1e-30))


if __name__ == '__main__':
    out = {}
    for name in ROWS:
        r = dict(A=0.0, B=0.0, C=0.0, B_as_given=0.0)
        for seed in (1, 2):
            c = TE.cfg(seed)
            SCALE[0] = 0.0
            sA, sk = sums(name, TE.Env(c, KP=1.0, V=37.0, L=1e-9))
            sB, _ = sums(name, TE.Env(c, KP=1.0, V=E, L=1e-9))
            sC, _ = sums(name, TE.Env(c, KP=1.0, V=37.0, L=1e-9))
            r['A'] = max(r['A'], rel(*sA['A']))
            r['B'] = max(r['B'], rel(*sB['B'][:2]))
            r['B_as_given'] = max(r['B_as_given'], rel(sB['B'][2], sB['B'][1]) if abs(sB['B'][1]) > 1e-30 or abs(sB['B'][2]) > 1e-30 else 0.0)
            r['C'] = max(r['C'], rel(*sC['C']))
        # text comparison of the terms that cannot be evaluated
        mine_txt = sorted(t for w, n, t in sk if w == 'mine')
        your_txt = sorted(t for w, n, t in sk if w == 'yours')
        r['not_evaluable_identical'] = [a == b for a, b in zip(mine_txt, your_txt)] if len(mine_txt) == len(your_txt) else 'count differs %d/%d' % (len(mine_txt), len(your_txt))
        r['n_not_evaluable'] = len(mine_txt)
        out[name] = r
        print('%-14s A %.1e  B %.1e  (BFKL.tex as given %.1e)  C %.1e  non-evaluable identical: %s'
              % (name, r['A'], r['B'], r['B_as_given'], r['C'], r['not_evaluable_identical']))
    json.dump(out, open(os.path.join(HERE, 'check_rows.json'), 'w'), indent=1)
