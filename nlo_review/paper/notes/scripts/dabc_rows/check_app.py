# Operator-level check of the printed row-by-row decompositions (app:fin_sym) of the cubic colour operators:
#   ordered product (colour structure of the row, parsed from the TeX source of the finite terms, rows.py)
#   == printed tensor x [t {rho rho rho}_sym + printed commutator (two-rho) terms + printed one-rho terms],
# with rho^a(site) = SU(3) generators on 3 lattice sites (explicit matrices), delta = Kronecker in the sites,
# for every assignment of the charge points to sites, random Wilson lines (also merged at coincident sites).
# usage: python3 check_app.py <app_finite.tex>   (or a file containing the same equations)
import sys, os, re, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); DABC = os.path.dirname(HERE); SING = os.path.join(os.path.dirname(DABC), 'singlet')
sys.path.insert(0, DABC); sys.path.insert(0, SING)
import numpy as np
import fasteinsum
from rows import CUB
from engine import parse, eval_terms, Group
src = open(sys.argv[1]).read()
G = Group('su3'); n = G.n; P = 3; d = 3
def site_ops():
    M = np.zeros((P, n, d**P, d**P), complex)
    for p in range(P):
        for a in range(n):
            m = np.array([[1.0]], complex)
            for q in range(P): m = np.kron(m, G.t[a] if q == p else np.eye(d))
            M[p, a] = m
    return M
M = site_ops()
# sanity: [rho^a(p), rho^b(q)] = i f^{abc} rho^c(p) delta_pq
cm = np.einsum('aij,bjk->abik', M[0], M[0]) - np.einsum('bij,ajk->abik', M[0], M[0])
assert np.allclose(cm, 1j * np.einsum('abc,cij->abij', G.f, M[0]))
assert np.allclose(np.einsum('aij,bjk->abik', M[0], M[1]), np.einsum('bij,ajk->abik', M[1], M[0]))

def braces(t):  # TeX of the paper -> parser format (braced indices)
    t = re.sub(r"\\rho\^([a-z])", r"\\rho^{\1}", t)
    return t
def prod3(T, o1, o2, o3):
    return np.einsum('abc,aij,bjk,ckl->il', T, o1, o2, o3, optimize=True)
def sym3(T, ops):
    out = 0
    for perm in itertools.permutations(range(3)):
        # operator order perm: ops[perm[0]] ops[perm[1]] ops[perm[2]]; tensor indices follow ops
        letters = 'abc'
        sub = 'abc,' + ','.join(letters[q] + ij for q, ij in zip(perm, ('ij', 'jk', 'kl'))) + '->il'
        out = out + np.einsum(sub, T, *[ops[q] for q in perm], optimize=True)
    return out / 6
def sym2(K, o1, o2):
    return 0.5 * (np.einsum('ab,aij,bjk->ik', K, o1, o2, optimize=True) + np.einsum('ab,bij,ajk->ik', K, o2, o1, optimize=True))

# ---------- parse the printed blocks ----------
blocks = re.findall(r"\\paragraph\{(Row [^}]*?(?:\{[^}]*\}[^}]*?)*)\}(.*?)\\begin\{align\}\n(.*?)\\label\{eq:(app_sym_[A-Za-z0-9]+)\}", src, re.S)
NAMES = {'app_sym_IIb': ['I.II.N2b'], 'app_sym_IIII': ['I.III(z)', 'I.III(w)'], 'app_sym_IV1': ['I.V.N1'],
         'app_sym_III1a': ['II.I.N1a'], 'app_sym_III1b': ['II.I.N1b'], 'app_sym_IIIII': ['III.II'],
         'app_sym_IIIIII': ["III.III.1(y')", 'III.III.1(z)'], 'app_sym_IIIV': ["III.V(y')", 'III.V(z)'], 'app_sym_IIIVI': ["III.VI(y')", 'III.VI(z)']}
assert sorted(b[3] for b in blocks) == sorted(NAMES), [b[3] for b in blocks]
PT = r"([a-z]'?)"
worst = 0; ncases = 0
for head, sent, body, lab in blocks:
    lines = body.split('\n')
    # tensor definition: from '&\equiv ' up to the line ending with '\,,\notag\\'
    k = next(i for i, l in enumerate(lines) if l.rstrip().endswith(r'\,,\notag\\'))
    deftex = '\n'.join(lines[:k + 1])
    cname, cdef = deftex.split(r'&\equiv ', 1)
    cdef = cdef.rsplit(r'\,,\notag\\', 1)[0]
    rest = '\n'.join(lines[k + 1:])
    cn = re.escape(cname.strip())
    # LHS: type and points
    m = re.match(cn + r"\\,(\\rho\^a\(" + PT + r"\)\\rho\^b\(" + PT + r"\)\\rho\^c\(" + PT + r"\)|\\rho\^a\(" + PT + r"\)\\\{\\rho\^b\(" + PT + r"\),\\rho\^c\(" + PT + r"\)\\\}|\\\{\\rho\^a\(" + PT + r"\),\\rho\^b\(" + PT + r"\)\\\}\\rho\^c\(" + PT + r"\))\n", rest)
    assert m, (lab, rest[:200])
    if m.group(2): typ, pts = 'ord', m.group(2, 3, 4)
    elif m.group(5): typ, pts = 'A23', m.group(5, 6, 7)
    else: typ, pts = 'A12', m.group(8, 9, 10)
    rest = rest[m.end():]
    # the operator in the sentence must be C rho rho rho with the same points/type
    lhs_s = m.group(1)
    assert (cname.strip() + r'\,' + lhs_s) in sent, (lab, 'sentence')
    # symmetrised term
    ms = re.match(r"&=(2\\,)?" + cn + r"\\,\\\{\\rho\^a\(" + PT + r"\)\\rho\^b\(" + PT + r"\)\\rho\^c\(" + PT + r"\)\\\}_\{\\rm sym\}\\notag\\\\\n", rest)
    assert ms, (lab, rest[:120])
    tsym = 2 if ms.group(1) else 1
    assert ms.group(2, 3, 4) == pts
    rest = rest[ms.end():]
    # two-rho bracket
    mt = re.match(r"&\\quad\+(\\frac i2|i)\\," + cn + r"\\,\\Big\[(.*?)\\Big\]\\notag\\\\\n", rest, re.S)
    assert mt, (lab, rest[:150])
    c2 = 0.5j if mt.group(1) == r'\frac i2' else 1j
    two = re.findall(r"f\^\{([a-z])([a-z])([a-z])\}\\,\\delta\^\{\(2\)\}\(" + PT + "-" + PT + r"\)\\,\\\{\\rho\^\{([a-z])\}\(" + PT + r"\)\\rho\^\{([a-z])\}\(" + PT + r"\)\\\}_\{\\rm sym\}", mt.group(2))
    assert len(two) == (3 if typ == 'ord' else 2), (lab, two)
    rest = rest[mt.end():]
    # one-rho
    mo = re.match(r"&\\quad([+-])\\frac16\\," + cn + r"\\big\[f\^\{([a-z]{3})\}f\^\{([a-z]{3})\}([+-])f\^\{([a-z]{3})\}f\^\{([a-z]{3})\}\\big\]\\notag\\\\\n&\\qquad\\times\\delta\^\{\(2\)\}\(" + PT + "-" + PT + r"\)\\,\\delta\^\{\(2\)\}\(" + PT + "-" + PT + r"\)\\,\\rho\^s\(" + PT + r"\)\\,\.$", rest.strip())
    assert mo, (lab, rest)
    s1 = 1 if mo.group(1) == '+' else -1; s2 = 1 if mo.group(4) == '+' else -1
    ff = [(mo.group(2), mo.group(3), s1), (mo.group(5), mo.group(6), s1 * s2)]
    dd = [(mo.group(7), mo.group(8)), (mo.group(9), mo.group(10))]; ps = mo.group(11)
    for name in NAMES[lab]:
        tex, ctyp, cpts, _ = CUB[name]
        assert ctyp == {'ord': 'ord', 'A23': 'anti23', 'A12': 'anti12'}[typ] and tuple(cpts) == tuple(pts), (lab, name)
        sub = {}
        if '(u)' in cname: sub['u'] = 'z' if name.endswith('(z)') else 'w'
        if '(v)' in cname: sub['v'] = 'z' if name.endswith('(z)') else "y'"
        cd = cdef
        for k_, v_ in sub.items(): cd = re.sub(r"\(" + k_ + r"\)", "(" + v_ + ")", cd)
        for merge in (False, True):
            for assign in itertools.product(range(P), repeat=3):
                s_of = dict(zip(pts, assign))
                if len(set(pts)) < 3: continue
                rng = np.random.default_rng(hash((lab, name, assign, merge)) % 2**32)
                U = {}
                for q in ["x", "x'", "y", "y'", "z", "w", "w'"]: U[q] = G.randU(rng)
                if merge:
                    for p1_, p2_ in itertools.combinations(pts, 2):
                        if s_of[p1_] == s_of[p2_]: U[p2_] = U[p1_]
                # LHS: the row's ordered colour operator
                ev = eval_terms(G, parse(tex), U)
                L = 0
                for seq, T in ev.items():
                    L = L + prod3(T, M[s_of[seq[0]]], M[s_of[seq[1]]], M[s_of[seq[2]]])
                # printed tensor
                C = eval_terms(G, parse(braces(cd) + r"\rho^{a}(%s)\rho^{b}(%s)\rho^{c}(%s)" % pts), U)[tuple(pts)]
                ops = [M[s_of[q]] for q in pts]
                R = tsym * sym3(C, ops)
                for (i1, i2, r_, pa, pb, j1, q1, j2, q2) in two:
                    if s_of[pa] != s_of[pb]: continue
                    K = np.einsum('abc,%s%s%s->%s%s' % (i1, i2, r_, j1, j2), C, G.f)
                    R = R + c2 * sym2(K, M[s_of[q1]], M[s_of[q2]])
                if s_of[dd[0][0]] == s_of[dd[0][1]] and s_of[dd[1][0]] == s_of[dd[1][1]]:
                    Ks = sum(sg * np.einsum('abc,%s,%s->s' % (fa, fb), C, G.f, G.f) for fa, fb, sg in ff)
                    R = R + (1 / 6) * np.einsum('s,sij->ij', Ks, M[s_of[ps]])
                r = np.abs(L - R).max() / max(1, np.abs(L).max()); worst = max(worst, r); ncases += 1
    print('%-15s %-14s type %-3s points %s: checked' % (lab, ','.join(NAMES[lab]), typ, pts))
print('APPENDIX decompositions vs ordered operator products: %d cases, worst rel. deviation %.1e' % (ncases, worst))
