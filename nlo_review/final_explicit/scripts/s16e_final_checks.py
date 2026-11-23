# s16e_final_checks.py -- checks of the PRINTED file Final_terms_explicit_2_fixed.tex:
#  (1) Sec. 8.1: every colour block \Bigg\{...\Bigg\} (factorised C^{abc} or two-rho factor), TeX parsed, == the verified
#      colour tensor of the row (s3_colour.STRUCT / s16c.refq); the Row III.III second-part factors L(v) x R parsed;
#  (2) Sec. 8.2: every printed coefficient x trace bracket == t C.d (three-rho), g_pair C.f|_{p_i=p_j} (contact),
#      delta-contraction (genuine two-rho);
#  (3) operator level (SU(3) fundamental charges on sites, random C^{abc}): Eqs. (D-ord) as corrected, (D-r12), (D-r21);
#  (4) the missing brackets of Sec. 6 (author's file): as typeset vs intended, and with the red braces;
#  (5) no p+ left in Sec. 8; counts of xi integrals.
# Run: python3 -I s16e_final_checks.py -> s16e_final_checks.log
import os, re, sys, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE); sys.path.insert(0, HERE)
from s16_lib import SUN, ev, evs, bracket_after, randUs
from s3_colour import STRUCT, E
from s16c_sec8_colour import refq, GEN, TT
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
SF = {s[0]: s for s in STRUCT}
doc = open(os.path.join(TOP, 'Final_terms_explicit_2_fixed.tex')).read()
F1 = doc[doc.index(r'\label{sec:F1}'):doc.index(r'\label{sec:F2}')]
F2 = doc[doc.index(r'\label{sec:F2}'):doc.index(r'\label{sec:F3}')]
def eqblock(T, lab):
    i = T.index(r'\label{%s}' % lab); return T[T.rfind(r'\begin{align}', 0, i):i]
def bigg_blocks(blk):
    res = []; i = 0
    while True:
        j = blk.find(r'\Bigg\{', i)
        if j < 0: return res
        depth = 0; k = j
        while True:
            if blk.startswith(r'\Bigg\{', k): depth += 1; k += 6; continue
            if blk.startswith(r'\Bigg\}', k):
                depth -= 1; k += 6
                if depth == 0: break
                continue
            k += 1
        res.append(blk[j + 6:k - 6]); i = k
P1 = [('eq:H-I.II-3', ['I.II 2b']), ('eq:H-I.II-2', ['I.II 2b']), ('eq:H-I.III-3', ['I.III C_D (u=z)', 'I.III C_A (u=w)']),
      ('eq:H-I.III-2', ['I.III C_D (u=z)', 'I.III C_A (u=w)']), ('eq:H-I.V-3', ['I.V first']), ('eq:H-I.V-2', ['I.V first']),
      ('eq:H-I.V-2o', ['Q:I.V N2']), ('eq:H-II.I-3', ['II.I 1a', 'II.I 1b']), ('eq:H-II.I-2', ['II.I 1a', 'II.I 1b']),
      ('eq:H-III.II-3', ['III.II']), ('eq:H-III.II-2', ['III.II']), ('eq:H-III.III-3', ["III.III L(y')", 'III.III L(z)']),
      ('eq:H-III.III-2', ["III.III L(y')", 'III.III L(z)']), ('eq:H-III.V-3', ["III.V L(y')", 'III.V L(z)']),
      ('eq:H-III.V-2', ["III.V L(y')", 'III.V L(z)']), ('eq:H-III.VI-3', ["III.VI L(y')", 'III.VI L(z)']), ('eq:H-III.VI-2', ["III.VI L(y')", 'III.VI L(z)'])]
worst = {}
def upd(k, v): worst[k] = max(worst.get(k, 0), v)
n1 = 0
for N in (3, 4):
    S = SUN(N); rng = np.random.default_rng(1900 + N)
    for _ in range(3):
        U = randUs(S, rng)
        for lab, sts in P1:
            bl = [b for b in bigg_blocks(eqblock(F1, lab)) if 'U^{' in b and 'mathcal' not in b]; assert len(bl) == len(sts), (lab, len(bl))
            for b, st in zip(bl, sts):
                v, idx = ev(b, S, U)
                if st.startswith('Q:'):
                    ref = refq(st[2:], S, U); order = ["b'", 'f']
                else:
                    ref = SF[st][1](S, U); order = ['a', 'b', 'c']
                v = np.transpose(v, [idx.index(q) for q in order])
                upd('(1) Sec. 8.1 printed colour blocks == verified tensors', np.abs(v - ref).max() / max(1, np.abs(ref).max())); n1 += 1
        # II.I second part (5 equations, unlabelled) and III.III second part
        g = F1[F1.index(r'\label{eq:H-II.I-2}'):F1.index(r'\label{eq:H-II.IVrem}')]
        pieces = re.findall(r'(\\bigg\[\\Big\[U\^\{ac\'\}\(w\'\).*?)\\rho\^\{d\'\}\(x\'\)', g, re.S)
        assert len(pieces) == 5, len(pieces)
        for pc in pieces:
            v, idx = ev(pc, S, U); v = np.transpose(v, [idx.index(q) for q in ["d'", 'd']]); ref = refq('II.I second', S, U)
            upd('(1) Sec. 8.1 Group II Row I second part (5 eqs) colour == verified', np.abs(v - ref).max() / np.abs(ref).max())
        g = eqblock(F1, 'eq:H-III.III-2o')
        Ls = re.findall(r"\\,(f\^\{c'd'a\}\\,U\^\{d'b\}\(z\)\\Big\[U\^\{c'e'\}\((y'|z)\)-U\^\{c'e'\}\(x'\)\\Big\])\\rho\^\{e'\}\(x'\)", g)
        Rs = re.findall(r"(\\bigg\[\\Big\[U\^\{ac\}\(w\)-U\^\{ac\}\(z\)\\Big\]f\^\{cbd\}-f\^\{ace\}\\,U\^\{cb\}\(z\)\\Big\[U\^\{ed\}\(x\)-U\^\{ed\}\(z\)\\Big\]\\bigg\])\\rho\^\{d\}\(x\)", g)
        assert len(Ls) == 5 and len(Rs) == 4, (len(Ls), len(Rs))
        for Lt, v_ in Ls:
            v, idx = ev(Lt + Rs[0], S, U); v = np.transpose(v, [idx.index(q) for q in ["e'", 'd']])
            ref = refq('III.III second L' + ("(y')" if v_ == "y'" else '(z)'), S, U)
            upd('(1) Sec. 8.1 Group III Row III second part L(v) x R == verified', np.abs(v - ref).max() / np.abs(ref).max())
pr('(1) Sec. 8.1: %d colour blocks evaluated (17 equations x 3 configurations x SU(3), SU(4))' % n1)
# (2) Sec. 8.2
PT = {"x'": 'xp', "w'": 'wp', "y'": 'yp', 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w'}
PLAN = [('eq:S-I.II-3', ['G:I.II 2b']), ('eq:S-I.II-2', ['C:I.II 2b'] * 2), ('eq:S-I.II-2o', ['Q:I.II N1']),
        ('eq:S-I.III-3', ['G:I.III C_D (u=z)', 'G:I.III C_A (u=w)']), ('eq:S-I.III-2', ['C:I.III C_D (u=z)'] * 3 + ['C:I.III C_A (u=w)'] * 3),
        ('eq:S-I.V-3', ['G:I.V first']), ('eq:S-I.V-2', ['C:I.V first'] * 2), ('eq:S-I.V-2o', ['Q:I.V N2', 'Q:I.V 11/6']),
        ('eq:S-II.I-3', ['G:II.I 1a', 'G:II.I 1b']), ('eq:S-II.I-2', ['C:II.I 1a'] * 2 + ['C:II.I 1b'] * 2), ('II.I genuine', ['Q:II.I second'] * 5),
        ('eq:S-II.IVrem', ['Q:II.IV', 'Q:II.IV rem']), ('eq:S-III.II-3', ['G:III.II']), ('eq:S-III.II-2', ['C:III.II'] * 3),
        ('eq:S-III.III-3', ["G:III.III L(y')", 'G:III.III L(z)']), ('eq:S-III.III-2', ["C:III.III L(y')"] * 2 + ['C:III.III L(z)'] * 2),
        ('eq:S-III.III-2o', ["Q:III.III second L(y')"] * 3 + ['Q:III.III second L(z)'] * 2),
        ('eq:S-III.V-3', ["G:III.V L(y')", 'G:III.V L(z)']), ('eq:S-III.V-2', ["C:III.V L(y')"] * 3 + ['C:III.V L(z)'] * 3),
        ('eq:S-III.VI-3', ["G:III.VI L(y')", 'G:III.VI L(z)']), ('eq:S-III.VI-2', ["C:III.VI L(y')"] * 3 + ['C:III.VI L(z)'] * 3)]
pat = re.compile(r"\\Big\((?P<c3>[^()]*?)\\Big\)\\,\\sigma\((?P<oa>[^()]*)\)"
                 r"|\\Big\((?P<cc>[^()]*?)\\Big\)\\delta\^\{\(2\)\}\((?P<pq>[^()]*)\)\\,n\((?P<na>[^()]*)\)"
                 r"|\\frac\{n\((?P<qa>[^()]*)\)\}\{N_\{c\}\^\{2\}-1\}")
def blk2(lab):
    if lab == 'II.I genuine':
        a = F2.index(r'\label{eq:S-II.I-2}'); return F2[a:F2.index(r'\paragraph*{Row II}', a)]
    return eqblock(F2, lab)
n2 = 0
for N in (3, 4):
    S = SUN(N); rng = np.random.default_rng(2000 + N)
    for lab, plan in PLAN:
        blk = blk2(lab); ms = list(pat.finditer(blk)); assert len(ms) == len(plan), (lab, len(ms))
        for m, item in zip(ms, plan):
            kind, st = item.split(':', 1)
            j = m.end()
            while blk[j] in ' \\,' and not blk.startswith(r'\Big', j) and not blk.startswith(r'{\rm', j): j += 1
            if blk.startswith(r'{\rm Tr}', j):
                cont, end = bracket_after(blk, j + 8); br = r'{\rm Tr}\{' + cont + r'\}'
            else:
                br, end = bracket_after(blk, j)
            for _ in range(3):
                U = randUs(S, rng)
                if kind == 'G':
                    name, fn, typ, pts = SF[st]; exp = TT[typ] * E('abc,abc->', fn(S, U), S.d); got = evs(m.group('c3'), S, U) * evs(br, S, U)
                elif kind == 'C':
                    name, fn, typ, pts = SF[st]; p, q = [PT[x.strip()] for x in m.group('pq').split('-')]
                    names = {PT[x] for x in re.findall(r"\((x'|y'|w'|[xyzw])\)", br)}
                    fr, to = (p, q) if p not in names else (q, p)
                    pair = ''.join(sorted('%d%d' % (pts.index(p) + 1, pts.index(q) + 1)))
                    U2 = dict(U); U2[fr] = U[to]
                    exp = GEN[typ][pair] * E('abc,abc->', fn(S, U2), S.f); got = evs(m.group('cc'), S, U) * evs(br, S, U)
                else:
                    exp = np.trace(refq(st, S, U)); got = evs(br, S, U)
                upd('(2) Sec. 8.2 printed coefficient x bracket == average of the Sec. 8.1 colour factor', abs(got - exp) / max(1, abs(exp)))
            n2 += 1
pr('(2) Sec. 8.2: %d printed singlet items evaluated (x 3 configurations, SU(3) and SU(4))' % n2)
# (3) operator level
lam = np.zeros((8, 3, 3), complex)
lam[0][0, 1] = lam[0][1, 0] = 1; lam[1][0, 1] = -1j; lam[1][1, 0] = 1j; lam[2][0, 0] = 1; lam[2][1, 1] = -1
lam[3][0, 2] = lam[3][2, 0] = 1; lam[4][0, 2] = -1j; lam[4][2, 0] = 1j; lam[5][1, 2] = lam[5][2, 1] = 1
lam[6][1, 2] = -1j; lam[6][2, 1] = 1j; lam[7] = np.diag([1, 1, -2]) / np.sqrt(3)
t = lam / 2
fabc = np.real(-2j * np.einsum('aij,bjk,cki->abc', t, t, t) + 2j * np.einsum('bij,ajk,cki->abc', t, t, t))
def site_ops(ns):
    I = np.eye(3); o = []
    for s_ in range(ns):
        oo = []
        for a in range(8):
            m = np.array([[1.0]])
            for q in range(ns): m = np.kron(m, t[a] if q == s_ else I)
            oo.append(m)
        o.append(np.array(oo))
    return o
def sym(*Ms):
    perms = list(itertools.permutations(range(len(Ms))))
    return sum(np.linalg.multi_dot([Ms[i] for i in p]) if len(Ms) > 2 else Ms[p[0]] @ Ms[p[1]] for p in perms) / len(perms)
rng = np.random.default_rng(7); C = rng.normal(size=(8, 8, 8))
ff12 = np.einsum('abc,abr,rcs->s', C, fabc, fabc); ff23 = np.einsum('abc,bcr,ras->s', C, fabc, fabc); ff13 = np.einsum('abc,acr,rbs->s', C, fabc, fabc)
dev = {'ord': 0, 'r12': 0, 'r21': 0}
for assign in [(0, 1, 2), (0, 0, 1), (0, 1, 1), (0, 1, 0), (0, 0, 0)]:
    R = site_ops(max(assign) + 1); A, B, Cc = [R[s_] for s_ in assign]
    d12, d23, d13 = [float(assign[i] == assign[j]) for i, j in ((0, 1), (1, 2), (0, 2))]
    tri = lambda X, Y, Z: sum(C[a, b, c] * X[a] @ Y[b] @ Z[c] for a in range(8) for b in range(8) for c in range(8))
    S3 = sum(C[a, b, c] * sym(A[a], B[b], Cc[c]) for a in range(8) for b in range(8) for c in range(8))
    def ct(T2, X, Y): return sum(T2[r, c] * sym(X[r], Y[c]) for r in range(8) for c in range(8))
    c12 = ct(np.einsum('abc,abr->rc', C, fabc), A, Cc) * d12
    c23 = ct(np.einsum('abc,bcr->ar', C, fabc), A, B) * d23
    c13 = ct(np.einsum('abc,acr->rb', C, fabc), A, B) * d13
    one = lambda v: sum(v[s_] * A[s_] for s_ in range(8)) * d12 * d13
    lhs = tri(A, B, Cc); rhs = S3 + 0.5j * (c12 + c23 + c13) - one(ff12 - ff23) / 6
    dev['ord'] = max(dev['ord'], np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs))
    lhs = tri(A, B, Cc) + sum(C[a, b, c] * A[a] @ Cc[c] @ B[b] for a in range(8) for b in range(8) for c in range(8))
    rhs = 2 * S3 + 1j * (c12 + c13) - one(ff12 + ff13) / 6
    dev['r12'] = max(dev['r12'], np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs))
    lhs = tri(A, B, Cc) + sum(C[a, b, c] * B[b] @ A[a] @ Cc[c] for a in range(8) for b in range(8) for c in range(8))
    ff_r21 = np.einsum('abc,bcr,ras->s', C, fabc, fabc) + np.einsum('abc,acr,rbs->s', C, fabc, fabc)
    rhs = 2 * S3 + 1j * (c23 + c13) + one(ff_r21) / 6
    dev['r21'] = max(dev['r21'], np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs))
pr('(3) operator level, SU(3) charges on 1-3 sites, random C^{abc}, all coincidence patterns: (D-ord) as corrected %.1e, (D-r12) %.1e, (D-r21) %.1e' % (dev['ord'], dev['r12'], dev['r21']))
# (4) missing brackets of Sec. 6 (author's file)
au = open(os.path.join(TOP, 'Final_terms_explicit_2_author.tex')).read()
S = SUN(3); U = randUs(S, np.random.default_rng(11))
cases = [('D-II.I-3 (1a)', r"2\,U^{ga}(w')\notag \\ & \Big[f^{kbc}\,U^{gk}(w)-f^{gkl}\,U^{kb}(z)\,U^{lc}(x)\Big]-U^{ga}(x')\,U^{mh}(z)\,U^{mb}(y)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]", 'II.I 1a', 2),
         ('D-III.II-3', r"U^{ma}(y)\,U^{gb}(x')\notag \\ & U^{mh}(z)\,\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]-U^{gb}(x')\Big[f^{kac}\,U^{gk}(w)-f^{gkl}\,U^{ka}(z)\,U^{lc}(x)\Big]", 'III.II', 1)]
for nm, tx, st, c in cases:
    assert tx in au
    ref = c * SF[st][1](S, U)
    v, idx = ev(tx, S, U); v = np.transpose(v, [idx.index(q) for q in 'abc'])
    coef = tx.split(r'U^{', 1)[0]
    v2, idx2 = ev(coef + r'\Big[' + tx[len(coef):] + r'\Big]', S, U); v2 = np.transpose(v2, [idx2.index(q) for q in 'abc'])
    pr('(4) %s: as typeset (coefficient on the first term only) rel dev %.2f; with the red braces %.1e' % (nm, np.abs(v - ref).max() / np.abs(ref).max(), np.abs(v2 - ref).max() / np.abs(ref).max()))
# (5)
sec8 = doc[doc.index(r'\label{sec:F}'):doc.index(r'\end{document}')]
eqs8 = sec8[sec8.index(r'\label{sec:F1}'):]
pr('(5) p^{+} in the equations of Sec. 8.1/8.2: %d occurrences; xi integrals: int_0^1 %d, int_0^inf %d, int_{k+/vee}^1 %d' % (eqs8.count('p^{+}'), eqs8.count(r'\int_{0}^{1}'), eqs8.count(r'\int_{0}^{\infty}\frac{d\xi}'), eqs8.count(r'\int_{k^{+}/\vee}^{1}')))
for k, v in worst.items(): pr('%-100s max rel dev %.1e' % (k, v))
open(os.path.join(HERE, 's16e_final_checks.log'), 'w').write('\n'.join(out) + '\n')
