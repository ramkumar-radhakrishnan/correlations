import pickle, time
import numpy as np
import sympy as sp
import users as US
import user4
import target as TG
import rowcmp
from fullpt import g1_pointwise, g2iv_pointwise
from tables import Book, tex_kkey
from engine import tex_colour

t0 = time.time()
S = rowcmp.user_sets()
U4 = user4.corrected_terms()
import re
REMOVED = re.compile(r'^(GI\.III\.\db|GI\.II\.[13][bd])\*?$')

GROUP = {'GI.II': 'G I, Row II', 'GI.III': 'G I, Row III', 'GI.IV': 'G I, Row IV', 'GI.V1': 'G I, Row V', 'GI.V2': 'G I, Row V',
         'GII.I1': 'G II, Row I', 'GII.I2': 'G II, Row I', 'GIII.I': 'G III, Row I', 'GIII.II': 'G III, Row II',
         'GIII.III': 'G III, Row III', 'GIII.V': 'G III, Row V', 'GIII.VI': 'G III, Row VI'}
B = Book()
for k, ts in S.items():
    if k.startswith(('A1', 'A2', 'B1', 'B2')):
        for t in ts:
            B.add(t, '4rho remainders A,B')
    elif k in GROUP:
        for t in ts:
            if not REMOVED.match(t['tag']):
                B.add(t, GROUP[k])
for t in U4:
    B.add(t, '4rho.tex 2rho')
for t in S['D']:
    B.add(t, 'D')
for t in g1_pointwise(S['G:G1_II']):
    B.add(t, 'G I, Row II')
for t in S['G:G1_V']:
    B.add(t, 'G I, Row V')
for t in S['G:G2_I']:
    B.add(t, 'G II, Row I')
for t in g2iv_pointwise(S['G:G2_IV'], 0.5):
    B.add(t, 'G II, Row IV')
for t in S['G:G3_III']:
    B.add(dict(t, c=-t['c']), 'G III, Row III')
for t in TG.jimwlk():
    B.add(dict(t, c=-t['c']), '-JIMWLK')
B.finalize()
print('built in %.0fs' % (time.time() - t0))
nk = len(B.data)
bad = 0
for key, nets in B.data.items():
    for net in nets:
        tot = sp.nsimplify(sum(net['groups'].values()))
        if tot != 0:
            bad += 1
    print(tex_kkey(key), ' networks', len(nets), ' JIMWLK networks', sum(1 for n in nets if '-JIMWLK' in n['groups']))
print('kernel structures', nk, 'networks', sum(len(v) for v in B.data.values()), 'nonzero totals', bad)
pickle.dump({k: [dict(col=n['col'], groups=dict(n['groups'])) for n in v] for k, v in B.data.items()}, open('book.pkl', 'wb'))
