# s16b_author2_sec7.py -- numerical check of the colour factors of Sec. 7 ("Average over the projectile") of the
# author's file Final_terms_explicit_2_author.tex (2026-11-23), as PRINTED (TeX parsed by s16_lib):
#   (3rho)  printed coefficient x printed trace bracket   ==  t * C^{abc} d^{abc}       (t = 2 for the anticommutator
#           types, 1 for ordered products; C^{abc} the verified colour tensor of the row, s3_colour.STRUCT)
#   (2rho, contact)  printed coefficient x bracket        ==  g_pair * C^{abc} f^{abc} |_{p_i = p_j}
#           (g_pair = i/2, i/2, -i/2 for 12, 23, 13 of ordered products; i, -i for the anticommutator types)
#   (2rho, genuine)  bracket                               ==  delta^{b'b}-contraction of the two-index colour factor
# random SU(3) and SU(4) adjoint Wilson lines.  Run: python3 -I s16b_author2_sec7.py -> s16b_author2_sec7.log
import os, re, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE); sys.path.insert(0, HERE)
from s16_lib import SUN, ev, evs, bracket_after, randUs
from s3_colour import STRUCT, E
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
SF = {s[0]: s for s in STRUCT}
GEN = {'ord': {'12': 0.5j, '23': 0.5j, '13': -0.5j}, 'r{rr}': {'12': 1j, '13': -1j}, '{rr}r': {'23': 1j, '13': -1j}}
TT = {'ord': 1, 'r{rr}': 2, '{rr}r': 2}
src = open(os.path.join(TOP, 'Final_terms_explicit_2_author.tex')).read()
sec7 = src[src.index(r'\section{Average over the projectile}'):src.index(r'\end{document}')]
PT = {"x'": 'xp', "w'": 'wp', "y'": 'yp', 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w'}
def items(block):
    """sequence of printed colour items in a block of Sec. 7: ('G', coef, args, bracket), ('C', coef, (p,q), nargs, bracket),
    ('Q', nargs, bracket)"""
    res = []
    pat = re.compile(r"\\Big\((?P<c3>[^()]*?)\\Big\)\\,o\((?P<oa>[^()]*)\)"
                     r"|\\Big\((?P<cc>[^()]*?)\\Big\)\\delta\^\{\(2\)\}\((?P<pq>[^()]*)\)\\,n\((?P<na>[^()]*)\)"
                     r"|\\frac\{n\((?P<qa>[^()]*)\)\}\{N_\{c\}\^\{2\}-1\}")
    for m in pat.finditer(block):
        j = m.end()
        while block[j] in ' \\,' and not block.startswith(r'\Big', j) and not block.startswith(r'\Bigg', j) and not block.startswith(r'{\rm', j): j += 1
        if block.startswith(r'{\rm Tr}', j):
            b0 = j + len(r'{\rm Tr}'); cont, end = bracket_after(block, b0); br = r'{\rm Tr}\{' + cont + r'\}'
        else:
            cont, end = bracket_after(block, j); br = cont
        if m.group('c3') is not None: res.append(('G', m.group('c3'), m.group('oa'), br))
        elif m.group('cc') is not None: res.append(('C', m.group('cc'), m.group('pq'), m.group('na'), br))
        else: res.append(('Q', m.group('qa'), br))
    return res
def eqblock(label):
    i = sec7.index(r'\label{%s}' % label); a = sec7.rfind(r'\begin{align}', 0, i); return sec7[a:i]
# genuine two-rho reference colour (delta-contracted): functions of (S, U)
def q_I2(S, U): return np.trace((U['xp'] - U['wp']).T @ (U['w'] - U['x']))
def q_I5(S, U):
    L = U['xp'] - U['wp']
    C = E('ab,af->bf', L, S.N * U['y']) - E('aB,adc,fbe,db,ce->Bf', L, S.f, S.f, U['x'], U['z']); return np.trace(C)
def q_I5b(S, U): return np.trace((U['xp'] - U['wp']).T @ (U['y'] - U['x']))
def q_II1(S, U):
    f = S.f
    C = (E('kbB,lbd,ak,al->Bd', f, f, U['wp'], U['w']) - E('kbB,acd,ak,cb,dE->BE', f, f, U['wp'], U['z'], U['x'])
         - E('akd,cbq,ac,kb,dB->Bq', f, f, U['w'], U['z'], U['xp']) + S.N * E('dB,dE->BE', U['xp'], U['x']))
    return np.trace(C)
def q_II4(S, U): return np.trace((U['yp'] - U['xp']).T @ (U['y'] - U['x']))
def q_II4r(S, U): return np.trace((U['z'] - U['xp']).T @ (U['z'] - U['x']))
def q_III3(v):
    def q(S, U):
        f = S.f; Lv = E('db,cB->cdbB', U['z'], U[v] - U['xp'])
        R2 = E('cbd,ac->abd', f, U['w']) - E('acd,cb,de->abe', f, U['z'], U['x'])
        return np.trace(E('cda,cdbB,abh->Bh', f, Lv, R2))
    return q
# expected sequences: (label, [(kind, struct-or-function)])
PLAN = [
    ('eq:E-I.II-3', [('G', 'I.II 2b')]),
    ('eq:E-I.II-2', [('C', 'I.II 2b'), ('C', 'I.II 2b')]),
    ('eq:E-I.II-2o', [('Q', q_I2)]),
    ('eq:E-I.III-3', [('G', 'I.III C_D (u=z)'), ('G', 'I.III C_A (u=w)')]),
    ('eq:E-I.III-2', [('C', 'I.III C_D (u=z)')] * 3 + [('C', 'I.III C_A (u=w)')] * 3),
    ('eq:E-I.V-3', [('G', 'I.V first')]),
    ('eq:E-I.V-2', [('C', 'I.V first')] * 2),
    ('eq:E-I.V-2o', [('Q', q_I5), ('Q', q_I5b)]),
    ('eq:E-II.I-3', [('G', 'II.I 1a'), ('G', 'II.I 1b')]),
    ('eq:E-II.I-2', [('C', 'II.I 1a')] * 2 + [('C', 'II.I 1b')] * 2),
    ('II.I genuine', [('Q', q_II1)] * 5),
    ('eq:E-II.IVrem', [('Q', q_II4), ('Q', q_II4r)]),
    ('eq:E-III.II-3', [('G', 'III.II')]),
    ('eq:E-III.II-2', [('C', 'III.II')] * 3),
    ('eq:E-III.III-3', [('G', "III.III L(y')"), ('G', 'III.III L(z)')]),
    ('eq:E-III.III-2', [('C', "III.III L(y')")] * 2 + [('C', 'III.III L(z)')] * 2),
    ('eq:E-III.III-2o', [('Q', q_III3('yp'))] * 3 + [('Q', q_III3('z'))] * 2),
    ('eq:E-III.V-3', [('G', "III.V L(y')"), ('G', 'III.V L(z)')]),
    ('eq:E-III.V-2', [('C', "III.V L(y')")] * 3 + [('C', 'III.V L(z)')] * 3),
    ('eq:E-III.VI-3', [('G', "III.VI L(y')"), ('G', 'III.VI L(z)')]),
    ('eq:E-III.VI-2', [('C', "III.VI L(y')")] * 3 + [('C', 'III.VI L(z)')] * 3),
]
def block_of(lab):
    if lab == 'II.I genuine':
        a = sec7.index(r'\label{eq:E-II.I-2}'); b = sec7.index(r'\subsubsection*{Row II}', a); return sec7[a:b]
    return eqblock(lab)
results = {}
nitems = 0
for N in (3, 4):
    S = SUN(N); rng = np.random.default_rng(1600 + N)
    for lab, plan in PLAN:
        its = items(block_of(lab))
        assert len(its) == len(plan), (lab, len(its), len(plan), [i[0] for i in its])
        for (kind, ref), it in zip(plan, its):
            assert it[0] == kind, (lab, kind, it[0])
            dev = 0; devb = 0; info = ''
            for _ in range(4):
                U = randUs(S, rng)
                if kind == 'G':
                    name, fn, typ, pts = SF[ref]
                    exp = TT[typ] * E('abc,abc->', fn(S, U), S.d)
                    got = evs(it[1], S, U) * evs(it[3], S, U)
                    info = 'coef %s, sigma args (%s)' % (it[1], it[2])
                elif kind == 'C':
                    name, fn, typ, pts = SF[ref]
                    p, q = [PT[s.strip()] for s in it[2].split('-')]
                    brpts = set(re.findall(r"\((x'|y'|w'|[xyzw])\)", it[4]))
                    # the eliminated point is the one of the delta that does not occur in the bracket
                    names = {PT[s]: s for s in brpts}
                    if p not in names: fr, to = p, q
                    elif q not in names: fr, to = q, p
                    else: fr, to = p, q  # both occur: ambiguous, test both below
                    i1, i2 = pts.index(p), pts.index(q)
                    pair = ''.join(sorted('%d%d' % (i1 + 1, i2 + 1)))
                    U2 = dict(U); U2[fr] = U[to]
                    Phi = E('abc,abc->', fn(S, U2), S.f)
                    exp = GEN[typ][pair] * Phi
                    got = evs(it[1], S, U) * evs(it[4], S, U2)
                    devb = max(devb, abs(evs(it[4], S, U2) - Phi) / max(1, abs(Phi)))
                    info = 'coef %s, delta(%s), n(%s), pair %s' % (it[1], it[2], it[3], pair)
                else:
                    exp = ref(S, U); got = evs(it[2], S, U); info = 'n(%s)' % it[1]
                dev = max(dev, abs(got - exp) / max(1, abs(exp)))
            key = (lab, len([k for k in results if k[0] == lab]) if N == 3 else None)
            results.setdefault((lab, its.index(it) if False else id(it)), None)
            pr('SU(%d) %-16s %-2s %-26s %-52s max rel dev %.1e' % (N, lab, kind, (ref if isinstance(ref, str) else ref.__name__), info, dev))
            nitems += 1
pr('items checked: %d (both N)' % nitems)
open(os.path.join(HERE, 's16b_author2_sec7.log'), 'w').write('\n'.join(out) + '\n')
