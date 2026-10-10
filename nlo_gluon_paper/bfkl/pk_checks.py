"""All numerical checks of the log(vee/k+) note -> checks.json"""
import env  # noqa
import json, sys
import sympy as sp
from collections import defaultdict
from core import Config, accumulate
from compare import report
import target as TG
import refP
import rowcmp
import userP as UP
import users as US
import lsets
import fourrho
from core import term, symmetrize
from acc3 import Config3, accumulate3, report3
from sources import SOURCES

out = {}


def dev2(us, rf, seed=5, symmetric=True):
    cfg = Config(seed, symmetric)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate(us, cfg, acc, False, label='a')
    accumulate([dict(t, c=-t['c']) for t in rf], cfg, acc, False, label='b')
    bad = report(acc, show=False)
    mx = max([abs(sum(d.values())) for d in acc.values()] + [0])
    return [len(acc), len(bad), mx]


def dev3(us, rf, seed=5):
    cfg = Config3(seed)
    acc = defaultdict(lambda: defaultdict(complex))
    accumulate3(us, cfg, acc, label='a')
    accumulate3([dict(t, c=-t['c']) for t in rf], cfg, acc, label='b')
    bad = report3(acc, show=False)
    mx = max([abs(sum(d.values())) for d in acc.values()] + [0])
    return [len(acc), len(bad), mx]


# ---------------- row by row, region P, two-rho and three-rho
GII = [x for x in SOURCES if x.key == 'GI.II'][0]
R2, R3 = refP.rows((2,)), refP.rows((3,))
S2 = rowcmp.user_sets()
W4 = {}
for c, row, lab, Us, word, kern, cc in fourrho.words():
    kk = [('KK', a, b, c_, d) for (a, b), (c_, d) in kern]
    t = term(c, [(u[1], u[2], u[3]) for u in Us], [], word, kk, '4rho')
    W4.setdefault(row, []).extend(symmetrize(t, (3,)))


def prow(keep, fixed):
    P = {k: UP.p_terms(k, keep) for k in UP.PSRC}
    four = (lambda r: sum((S2[k] for k in S2 if k.startswith('4:%s:' % r)), [])) if keep == (2,) else (lambda r: W4.get(r, []))
    g1 = P['GI.III.P'] + P['GI.IV.P'] + (sum((S2[k] for k in S2 if k.startswith('4:G1:')), []) if keep == (2,) else W4['G1'])
    if fixed:
        g1 = g1 + lsets.fixKh(keep) + US.neg(US.source_terms(GII, keep, monos=(1, 3), subs=(1, 3)))
    G = UP.genuineP() if keep == (2,) else defaultdict(list)
    R = R2 if keep == (2,) else R3
    return [
        ('Group I (all rows)', g1, sum((R[k] for k in R if k.startswith('G1-')), [])),
        ('Group II, Row I', P['GII.I1.P'] + P['GII.I2.P'] + four('G2-I') + G['G2_I'], R['G2-I']),
        ('Group II, Row II', four('G2-II'), R['G2-II']),
        ('Group II, Row III', four('G2-III'), R['G2-III']),
        ('Group II, Row IV', G['G2_IV'], R['G2-IV']),
        ('Group III, Row I', P['GIII.I.P'] + four('G3-I'), R['G3-I']),
        ('Group III, Row II', P['GIII.II.P'] + four('G3-II'), R['G3-II']),
        ('Group III, Row III', P['GIII.III.P'] + G['G3_III'], R['G3-III']),
        ('Group III, Row IV', four('G3-IV'), R['G3-IV']),
        ('Group III, Row V', P['GIII.V.P'], R['G3-V']),
        ('Group III, Row VI', P['GIII.VI.P'], R['G3-VI']),
    ]


rows = []
for keep, dev in (((2,), dev2), ((3,), dev3)):
    W, C = prow(keep, False), prow(keep, True)
    for (name, us, rf), (_, usc, _) in zip(W, C):
        if keep == (3,) and name == 'Group II, Row IV':
            continue
        a, b = dev(us, rf), dev(usc, rf)
        rows.append(dict(rho=keep[0], row=name, written=a, corrected=b))
        print(keep, name, a, b); sys.stdout.flush()
out['rowsP'] = rows

# ---------------- the log(vee/k+) total: states
B2, J2 = refP.bfkl(), TG.jimwlk()
Lw, Lc = lsets.sets((2,), False), lsets.sets((2,), True)


def total2(L, seed=11, symmetric=True, extra=()):
    cfg = Config(seed, symmetric)
    acc = defaultdict(lambda: defaultdict(complex))
    for k, ts in L.items():
        accumulate(ts, cfg, acc, False, label=k)
    for k, ts in extra:
        accumulate(ts, cfg, acc, False, label=k)
    accumulate([dict(t, c=-t['c']) for t in B2], cfg, acc, False, label='-B')
    accumulate(J2, cfg, acc, False, label='+J')
    bad = report(acc, show=False)
    return [len(acc), len(bad), max(abs(sum(d.values())) for d in acc.values())]


import cancel_data as CD
SRC = {s.key: s for s in SOURCES}
rowIII2 = US.source_terms(SRC['GI.III'], (2,), monos=(1,))
D2 = US.source_terms(CD.D_SOURCE, (2,))
G = US.genuine()
st = []
st.append(('as written', total2(Lw)))
st.append(('+ P1: Row III, sign of the K_h term', total2(Lw, extra=[('kh', lsets.fixKh((2,)))])))
st.append(('+ F3 carried over (Row III part 2 -> D)', total2(Lw, extra=[('kh', lsets.fixKh((2,))), ('r3', rowIII2), ('D', US.neg(D2))])))
st.append(('+ F5 carried over (Group III Row III part 2: sign)', total2(Lc)))
out['states2'] = st
for s in st:
    print(s)
# three-rho states
L3w, L3c = lsets.sets((3,), False), lsets.sets((3,), True)


def total3(L, seed=7, extra=()):
    cfg = Config3(seed)
    acc = defaultdict(lambda: defaultdict(complex))
    for k, ts in L.items():
        accumulate3(ts, cfg, acc, label=k)
    for k, ts in extra:
        accumulate3(ts, cfg, acc, label=k)
    accumulate3([dict(t, c=-t['c']) for t in refP.bfkl(keep=(3,))], cfg, acc, label='-B3')
    bad = report3(acc, show=False)
    return [len(acc), len(bad), max(abs(sum(d.values())) for d in acc.values())]


rowIII2_3 = US.source_terms(SRC['GI.III'], (3,), monos=(1,))
D3 = US.source_terms(CD.D_SOURCE, (3,))
st3 = [('as written', total3(L3w)),
       ('+ P1: Row III, sign of the K_h term', total3(L3w, extra=[('kh', lsets.fixKh((3,)))])),
       ('+ F3 carried over (Row III part 2 -> D)', total3(L3c))]
out['states3'] = st3
for s in st3:
    print(s)
# seeds and general U
out['seeds'] = [(sd, total2(Lc, seed=sd)) for sd in (21, 33)]
out['generalU'] = total2(Lc, seed=11, symmetric=False)
print(out['seeds'], out['generalU'])
json.dump(out, open('checks.json', 'w'), indent=1, default=float)
