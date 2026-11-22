# s11_author1_checks.py -- checks of the author's edited file Final_terms_explicit_1_author.tex (2026-11-22)
# (1) TeX level: every labelled equation of the author's file vs the verified file (git 5ae9e4f version, the one the
#     author edited; = HEAD up to the xi -> p+ rewrite of Row I.V), after normalisation (\ph, colours, spacing removed)
#     and substitution of the colour tensors C^{abc}_{row} of the verified Sec. D by their definitions.
# (2) Operator level: the symmetrisation identity eq:D-ord with the bracket as typeset (second (i/2)C^{abc} inside
#     the bracket) vs the intended one, SU(3) fundamental charges on sites.
# (3) Row I.V second part: xi form (author's Sec. B, xi undefined) == p+ form for xi = p+/k+, != for xi = 1-p+/k+.
# (4) Row I.III last term: typeset e^{-ik.(z-w)} (p+/k+) vs e^{-ik.(z-w) p+/k+}: the p+ integrals differ.
# (5) the averaging formula of the author's Sec. 7: d^{abc} projection coefficient, f part zero.
# Run: python3 -I s11_author1_checks.py  ->  s11_author1_checks.log
import os, re, sys, subprocess, difflib, itertools
import numpy as np
from scipy import integrate
HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE)
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
# ---------------- (1)
def strip_cmd(s, cmd):
    o = []; i = 0; key = '\\' + cmd + '{'
    while True:
        j = s.find(key, i)
        if j < 0: o.append(s[i:]); break
        o.append(s[i:j]); k = j + len(key); d = 1; st = k
        while d:
            d += {'{': 1, '}': -1}.get(s[k], 0); k += 1
        o.append(strip_cmd(s[st:k - 1], cmd)); i = k
    return ''.join(o)
def norm(s):
    s = strip_cmd(s, 'ph'); s = re.sub(r'\\color\{[a-zA-Z]+\}', '', s)
    s = re.sub(r'\\(quad|qquad|notag|left|right|big|Big|bigg|Bigg)(?![a-zA-Z])', '', s); s = re.sub(r'\\[,;!]', '', s)
    s = s.replace('\\\\', '').replace('&', ''); s = re.sub(r'\s+', '', s)
    for _ in range(5): s = re.sub(r'\{\{([^{}]*)\}\}', r'{\1}', s)
    s = re.sub(r'\\cdot\{(\([^{}]*\))\}', r'\\cdot\1', s); s = s.replace("^{\\prime}", "'")
    s = s.replace('\\times', '')            # \times between factors is typographical
    s = s.replace('{', '').replace('}', '') # brace grouping (colour groups) is typographical
    return s
def eqs(t):
    res = {}
    for m in re.finditer(r'\\begin\{align\}(.*?)\\end\{align\}', t, re.S):
        body = m.group(1); labs = re.findall(r'\\label\{([^}]*)\}', body); parts = re.split(r'\\label\{[^}]*\}', body)
        for l, p in zip(labs, parts): res[l] = p
    return res
base = subprocess.run(['git', 'show', '5ae9e4f:nlo_review/final_explicit/Final_terms_explicit.tex'], cwd=TOP, capture_output=True, text=True).stdout
auth = open(os.path.join(TOP, 'Final_terms_explicit_1_author.tex')).read()
defs = {m.group(1): m.group(2).strip() for m in re.finditer(r'(C\^\{abc\}_\{\{\\rm [^}]*\}\})&\\equiv (.*?)\\notag', base)}
bsub = base
for k, v in defs.items(): bsub = bsub.replace(k, v)
B, A = eqs(bsub), eqs(auth)
def heads(s): return re.sub(r'^\\frac\{d\^\{3\}(\\Delta)?N\^\{\\rm[^}]*\}(_\{[^}]*\})?\}\{d\^\{3\}k\}(\\bigg\|[^=]*)?=', '', s)
pr('(1) equation-by-equation comparison (author vs verified 5ae9e4f; %d colour tensors substituted):' % len(defs))
for l in sorted(set(A) | set(B), key=lambda l: (auth.find('label{%s}' % l) if l in A else 10 ** 7)):
    if l not in A: pr('   ONLY in verified : %s' % l); continue
    if l not in B: pr('   ONLY in author   : %s' % l); continue
    a, b = norm(heads(A[l].strip())), norm(heads(B[l].strip()))
    if a == b: pr('   same            : %s' % l); continue
    sm = difflib.SequenceMatcher(None, b, a, autojunk=False)
    d = ['%s %r -> %r' % (op, b[i1:i2][:80], a[j1:j2][:80]) for op, i1, i2, j1, j2 in sm.get_opcodes() if op != 'equal']
    pr('   DIFF            : %s : %s' % (l, ' | '.join(d)[:600]))
# merged Group II Row IV: author = verified tensor term + remnant
for x, y in (('eq:B-II4', 'eq:B-II4rem'), ('eq:D-II.IV', 'eq:D-II.IVrem')):
    rem = norm(B[y]); rem = rem[rem.index('=') + 1:]
    pr('   %s (author) == verified %s + %s : %s' % (y, x, y, norm(heads(B[x].strip())) + rem == norm(heads(A[y].strip()))))
# ---------------- (2) symmetrisation identity, operator level
lam = np.zeros((8, 3, 3), complex)
lam[0][0, 1] = lam[0][1, 0] = 1; lam[1][0, 1] = -1j; lam[1][1, 0] = 1j; lam[2][0, 0] = 1; lam[2][1, 1] = -1
lam[3][0, 2] = lam[3][2, 0] = 1; lam[4][0, 2] = -1j; lam[4][2, 0] = 1j; lam[5][1, 2] = lam[5][2, 1] = 1
lam[6][1, 2] = -1j; lam[6][2, 1] = 1j; lam[7] = np.diag([1, 1, -2]) / np.sqrt(3)
t = lam / 2
fabc = np.real(-2j * np.einsum('aij,bjk,cki->abc', t, t, t) + 2j * np.einsum('bij,ajk,cki->abc', t, t, t))
def site_ops(nsite):
    I = np.eye(3); ops = []
    for s_ in range(nsite):
        o = []
        for a in range(8):
            m = np.array([[1.0]])
            for q in range(nsite): m = np.kron(m, t[a] if q == s_ else I)
            o.append(m)
        ops.append(np.array(o))
    return ops
def sym(*Ms):
    perms = list(itertools.permutations(range(len(Ms)))); return sum(np.linalg.multi_dot([Ms[i] for i in p]) if len(Ms) > 2 else Ms[p[0]] @ Ms[p[1]] for p in perms) / len(perms)
rng = np.random.default_rng(1); C = rng.normal(size=(8, 8, 8))
worst = {'typeset': 0, 'intended': 0}
for assign in [(0, 1, 2), (0, 0, 1), (0, 1, 1), (0, 1, 0), (0, 0, 0)]:
    ns_ = max(assign) + 1; R = site_ops(ns_); rho = [R[s_] for s_ in assign]
    d12, d23, d13 = [float(assign[i] == assign[j]) for i, j in ((0, 1), (1, 2), (0, 2))]
    lhs = sum(C[a, b, c] * rho[0][a] @ rho[1][b] @ rho[2][c] for a in range(8) for b in range(8) for c in range(8))
    S3 = sum(C[a, b, c] * sym(rho[0][a], rho[1][b], rho[2][c]) for a in range(8) for b in range(8) for c in range(8))
    def c2(F_, i, j, k, dd):  # sum C f {rho_r(p_i) rho(p_k)}
        tot = 0
        for a, b, c in itertools.product(range(8), repeat=3):
            if abs(C[a, b, c]) == 0: continue
        return tot
    A12 = sum(C[a, b, c] * fabc[a, b, r] * sym(rho[0][r], rho[2][c]) for a in range(8) for b in range(8) for c in range(8) for r in range(8) if fabc[a, b, r] != 0) * d12
    A23 = sum(C[a, b, c] * fabc[b, c, r] * sym(rho[0][a], rho[1][r]) for a in range(8) for b in range(8) for c in range(8) for r in range(8) if fabc[b, c, r] != 0) * d23
    A13 = sum(C[a, b, c] * fabc[a, c, r] * sym(rho[0][r], rho[1][b]) for a in range(8) for b in range(8) for c in range(8) for r in range(8) if fabc[a, c, r] != 0) * d13
    ff = np.einsum('abr,rcs->abcs', fabc, fabc) - np.einsum('bcr,ras->abcs', fabc, fabc)
    one = np.einsum('abc,abcs->s', C, ff); O1 = -sum(one[s_] * rho[0][s_] for s_ in range(8)) / 6 * d12 * d13
    intended = S3 + 0.5j * (A12 + A23 + A13) + O1
    sc = np.einsum('abc,abc->', C, C) ** 0  # typeset: (i/2) C [f^{abr}.. + f^{bcr}.. + (i/2) C f^{acr}..] -> extra factor (i/2)C^{abc} on the f^{acr} term
    # the typeset bracket multiplies the third term by an extra (i/2)C^{abc} (indices a,b,c already summed): the
    # term becomes (i/2)^2 C^{abc} C^{abc} f^{acr}{..}, which is not even an allowed index structure; we evaluate it
    # with the repeated tensor summed once more (C^{abc} C^{abc} -> sum over the inner a,b,c of C^2):
    A13t = sum(C[a, b, c] ** 2 * fabc[a, c, r] * sym(rho[0][r], rho[1][b]) for a in range(8) for b in range(8) for c in range(8) for r in range(8) if fabc[a, c, r] != 0) * d13
    typeset = S3 + 0.5j * (A12 + A23 + 0.5j * A13t) + O1
    nrm = np.linalg.norm(lhs)
    worst['intended'] = max(worst['intended'], np.linalg.norm(lhs - intended) / nrm)
    worst['typeset'] = max(worst['typeset'], np.linalg.norm(lhs - typeset) / nrm if d13 else 0)
pr('(2) eq:D-ord, operator level (3x3 fundamental charges, random C^{abc}, assignments distinct/12/23/13/123):')
pr('    intended  (i/2)C[f^{abr}d12{}+f^{bcr}d23{}+f^{acr}d13{}] : max rel. deviation %.1e' % worst['intended'])
pr('    as typeset (second (i/2)C^{abc} inside the bracket)        : max rel. deviation %.2f (at p1=p3)' % worst['typeset'])
# ---------------- (3) Row I.V second part, xi vs p+ (kernel ratio identity at random points)
dev_p, dev_q = 0, 0
for _ in range(200):
    kp = rng.uniform(0.5, 2); p = rng.uniform(0, kp); YX, YZ, DY = rng.uniform(0.1, 3, 3)
    pf = kp / (p * YX + (kp - p) * YZ); coeffs_p = (p * (kp - p) / kp ** 2, (kp - p) / kp, p / kp)
    for xi, store in ((p / kp, 'p'), (1 - p / kp, 'q')):
        xf = 1 / (xi * YX + (1 - xi) * YZ); coeffs_x = (xi * (1 - xi), 1 - xi, xi)
        d = abs(xf - pf) / pf + sum(abs(a - b) for a, b in zip(coeffs_x, coeffs_p))
        if store == 'p': dev_p = max(dev_p, d)
        else: dev_q = max(dev_q, d)
pr('(3) Row I.V 2nd part: xi form == p+ form: max deviation %.1e for xi = p+/k+, %.2f for xi = 1-p+/k+' % (dev_p, dev_q))
# ---------------- (4) exponent typo in Row I.III last term (Sec. B)
dmax = 0
for s in (0.3, -1.7, 4.0):
    good = integrate.quad(lambda u: (np.cos(s * u) - 1) / u, 0, 1)[0] - 1j * integrate.quad(lambda u: np.sin(s * u) / u, 0, 1)[0]
    typo = np.exp(-1j * s) * integrate.quad(lambda u: (u - 0) / u, 0, 1)[0]   # [1/p+]_+ acting on (p+/k+): (u - 0)/u
    dmax = max(dmax, abs(good - typo))
pr('(4) Row I.III last term: int_0^1 du [1/u]_+ e^{-isu} vs typeset e^{-is} int_0^1 du [1/u]_+ u: |difference| up to %.2f (s=0.3,-1.7,4.0): not the same term' % dmax)
hits = [auth[:m.start()].count('\n') + 1 for m in re.finditer(r'e\^\{-ik\\cdot\(?[^{}]*\)\}\\frac\{p\^\{\+\}\}\{k\^\{\+\}\}', auth)]
pr('    occurrences in the author file of a recoil factor p+/k+ typeset OUTSIDE the exponent (e^{-ik.(..)}\\frac{p+}{k+}): lines %s' % hits)
# ---------------- (5) averaging coefficient
d3 = np.real(2 * np.einsum('aij,bjk,cki->abc', t, t, t) + 2 * np.einsum('bij,ajk,cki->abc', t, t, t))
N = 3; pr('(5) d^{abc}d^{abc} = %.6f, (N^2-4)(N^2-1)/N = %.6f; f^{abc}d^{abc} = %.1e: <{rho rho rho}_sym> = d^{abc} N/((N^2-4)(N^2-1)) <d^{efg}{rho rho rho}_sym> is correct (author omits the opening angle bracket)' % (np.sum(d3 * d3), (N * N - 4) * (N * N - 1) / N, abs(np.sum(fabc * d3))))
open(os.path.join(HERE, 's11_author1_checks.log'), 'w').write('\n'.join(out) + '\n')
