#!/usr/bin/env python3
"""Independent check of the one-gluon amplitude at O(g^3), Eq. (A3amp), and of the two-gluon amplitude
at O(g^2), Eq. (A2gamp), of Section 4.

Discrete model: M modes (a mode = colour index at one transverse point), [a_p,a^dag_q]=delta_pq,
non-commuting coefficients (formal words), eikonal S-matrix = real orthogonal matrix U on the mode index:
S^dag a_p S = U_pq a_q,  S^dag X S = Xbar for a rho-dependent coefficient X (independent symbols X~),
C_1 is rho independent (Cbar = C).  Omegabar^dag Omega|0> is normal ordered, the unitarity relations
(Section 3.4) are applied to the barred coefficients, and  sum_b U_ab <g_b|...>  and
sum_c U_ac <g_c g_b|...>  are compared with the expressions printed in the paper.
Run:  python3 -I check_amplitudes.py
"""
import numpy as np
from collections import defaultdict
from itertools import product
rng = np.random.default_rng(5)
M = (0, 1, 2)
Q, _ = np.linalg.qr(rng.normal(size=(3, 3))); U = Q            # real orthogonal
ORDER = {'A1': 1, 'A2': 1, 'a1': 3, 'a2': 3, 'N2': 2, 'B1': 2, 'B2': 2, 'B3': 2, 'C1': 1, 'C2': 1}
def sym(name, idx=(), dag=False):
    base = name.rstrip('~')
    if base in ('B1', 'B2'): idx = tuple(sorted(idx))
    if base in ('C1', 'C2'): idx = tuple(sorted(idx[:2]))+tuple(idx[2:])
    if base == 'N2': dag = False
    return (name, tuple(idx), dag)
def order(w): return sum(ORDER[s[0].rstrip('~')] for s in w)
def add(D, word, cre, ann, c=1.0):
    if order(word) > 3 or abs(c) < 1e-14: return
    D[(tuple(word), tuple(sorted(cre)), tuple(sorted(ann)))] += c
def contract(ann, cre):
    if not ann or not cre: return [(tuple(cre), tuple(ann), 1)]
    a, rest = ann[-1], ann[:-1]; out = []
    for (c2, a2, m) in contract(rest, cre): out.append((c2, a2+(a,), m))
    for i, q in enumerate(cre):
        if q == a:
            for (c2, a2, m) in contract(rest, cre[:i]+cre[i+1:]): out.append((c2, a2, m))
    return out
def mul(D1, D2):
    R = defaultdict(float)
    for (w1, c1, a1), x in D1.items():
        for (w2, c2, a2), y in D2.items():
            if order(w1)+order(w2) > 3: continue
            for (cl, al, m) in contract(a1, c2): add(R, w1+w2, c1+cl, al+a2, x*y*m)
    return R
def dagger(D):
    R = defaultdict(float)
    for (w, cre, ann), c in D.items():
        add(R, tuple((s[0], s[1], (not s[2]) if s[0].rstrip('~') != 'N2' else False) for s in reversed(w)), ann, cre, c)
    return R
def rot(modes):
    """all rotated mode tuples with weight prod U[p,q]"""
    out = [((), 1.0)]
    for p in modes: out = [(t+(q,), wgt*U[p, q]) for (t, wgt) in out for q in M]
    return out
def omega(bar):
    b = '~' if bar else ''
    D = defaultdict(float)
    def term(word, cre, ann):
        if not bar: add(D, word, cre, ann); return
        for (c2, w1) in rot(cre):
            for (a2, w2) in rot(ann): add(D, word, c2, a2, w1*w2)
    term((), (), ()); term((sym('N2'+b),), (), ())
    for p in M:
        for n1, n2 in (('A1', 'A2'), ('a1', 'a2')):
            term((sym(n1+b, (p,)),), (), (p,)); term((sym(n2+b, (p,)),), (p,), ())
    for p, q in product(M, M):
        term((sym('B1'+b, (p, q)),), (), (p, q)); term((sym('B2'+b, (p, q)),), (p, q), ()); term((sym('B3'+b, (p, q)),), (q,), (p,))
    for l, k, j in product(M, M, M):
        term((sym('C1', (l, k, j)),), (l, k), (j,)); term((sym('C2', (l, k, j)),), (j,), (k, l))     # C is rho independent
    return D
def substitute(poly, rules):
    changed = True
    while changed:
        changed = False; new = defaultdict(float)
        for w, c in poly.items():
            for i, s in enumerate(w):
                if s in rules:
                    for (c2, w2) in rules[s]:
                        nw = w[:i]+tuple(w2)+w[i+1:]
                        if order(nw) <= 3: new[nw] += c*c2
                    changed = True; break
            else: new[w] += c
        poly = defaultdict(float, {w: c for w, c in new.items() if abs(c) > 1e-13})
    return poly
b = '~'; rules = {}
for p in M:
    rules[sym('A1'+b, (p,), True)] = [(-1, (sym('A2'+b, (p,)),))]
    rules[sym('A1'+b, (p,), False)] = [(-1, (sym('A2'+b, (p,), True),))]
for l, k, j in product(M, M, M):
    rules[sym('C2', (l, k, j), True)] = [(-1, (sym('C1', (l, k, j)),))]
for p, q in product(M, M):
    if p <= q:
        r = [(-1, (sym('B2'+b, (p, q)),)), (.5, (sym('A2'+b, (p,)), sym('A2'+b, (q,)))), (.5, (sym('A2'+b, (q,)), sym('A2'+b, (p,))))]
        r += [(1, (sym('C1', (p, q, j)), sym('A2'+b, (j,)))) for j in M]
        rules[sym('B1'+b, (p, q), True)] = r
for p in M:
    r = [(-1, (sym('a2'+b, (p,)),)), (-1, (sym('N2'+b), sym('A2'+b, (p,)))), (1, (sym('A2'+b, (p,)), sym('N2'+b)))]
    for q in M:
        r += [(2, (sym('A1'+b, (q,)), sym('B2'+b, (p, q)))), (-1, (sym('B3'+b, (p, q), True), sym('A2'+b, (q,))))]
        r += [(-2, (sym('C1', (l, q, p), True), sym('B2'+b, (q, l)))) for l in M]
    rules[sym('a1'+b, (p,), True)] = r
# Omega|0>: only creation parts
Om0 = defaultdict(float, {k: v for k, v in omega(False).items() if not k[2]})
P = mul(dagger(omega(True)), Om0)
T1 = {m: defaultdict(float) for m in M}; T2 = {(c, d): defaultdict(float) for c in M for d in M}
for (w, cre, ann), c in P.items():
    if ann: continue
    if len(cre) == 1 and order(w) == 3: T1[cre[0]][w] += c
    if len(cre) == 2 and order(w) == 2:
        if cre[0] == cre[1]: T2[cre][w] += c
        else: T2[cre][w] += c/2; T2[cre[::-1]][w] += c/2      # symmetric coefficient of a^dag_c a^dag_d
for m in M: T1[m] = substitute(T1[m], rules)
for k in T2: T2[k] = substitute(T2[k], rules)
def diff(A, B):
    keys = set(A) | set(B); return max([abs(A.get(k, 0)-B.get(k, 0)) for k in keys] or [0])
A2 = lambda p, bar=False, dag=False: sym('A2'+('~' if bar else ''), (p,), dag)
# ---- one-gluon amplitude, Eq. (A3amp) ----------------------------------------------------------------
e1 = 0; nterms = 0
for a in M:
    lhs = defaultdict(float)
    for m in M:
        for w, c in T1[m].items(): lhs[w] += U[a, m]*c
    E = defaultdict(float)
    E[(A2(a, True), sym('N2'))] += -1                                  # A_N
    for bb in M: E[(sym('N2~'), A2(bb))] += U[a, bb]
    E[(sym('N2~'), A2(a, True))] += -1; E[(A2(a, True), sym('N2~'))] += 1
    for bb in M: E[(sym('a2', (bb,)),)] += U[a, bb]                    # A_{A(3)}
    E[(sym('a2~', (a,)),)] += -1
    for bb in M:                                                       # A_{B3}
        for c in M: E[(sym('B3~', (a, bb), True), A2(c))] += U[bb, c]
        E[(sym('B3~', (a, bb), True), A2(bb, True))] += -1
    for bb in M:                                                       # A_{AB2}
        for c in M:
            for d in M: E[(A2(bb, True, True), sym('B2', (c, d)))] += 2*U[a, d]*U[bb, c]
        E[(A2(bb, True, True), sym('B2~', (bb, a)))] += -2
    for d in M:                                                        # A_{CB2}
        for c in M:
            for bb in M:
                for e in M: E[(sym('C1', (d, c, a), True), sym('B2', (bb, e)))] += 2*U[d, bb]*U[c, e]
            E[(sym('C1', (d, c, a), True), sym('B2~', (d, c)))] += -2
    lhs = {k: v for k, v in lhs.items() if abs(v) > 1e-12}; nterms += len(lhs)
    e1 = max(e1, diff(lhs, E))
print(f"one-gluon amplitude O(g^3), Eq.(A3amp): {nterms} distinct operator words, max deviation {e1:.1e}")
# ---- two-gluon amplitude, Eq. (A2gamp) ---------------------------------------------------------------
e2 = 0; nterms = 0
for a in M:
    for bb in M:
        lhs = defaultdict(float)
        for c in M:
            for w, x in T2[(c, bb)].items(): lhs[w] += U[a, c]*x
        E = defaultdict(float)
        for c in M:
            E[(sym('B2', (bb, c)),)] += U[a, c]; E[(sym('B2~', (c, a)),)] += -U[c, bb]              # A_B
            E[(A2(c, True), A2(a, True))] += .5*U[c, bb]                                              # A_AA^(1)
            for d in M: E[(A2(d, True), A2(c))] += -.5*U[a, c]*U[d, bb]
            E[(A2(a, True), A2(c, True))] += .5*U[c, bb]                                              # A_AA^(2)
        E[(A2(a, True), A2(bb))] += -.5
        for d in M:                                                                                   # -A_CA
            for c in M:
                for e in M: E[(sym('C1', (d, a, c)), A2(e))] += -U[d, bb]*U[c, e]
                E[(sym('C1', (d, a, c)), A2(c, True))] += U[d, bb]
        lhs = {k: v for k, v in lhs.items() if abs(v) > 1e-12}; nterms += len(lhs)
        e2 = max(e2, diff(lhs, E))
print(f"two-gluon amplitude O(g^2), Eq.(A2gamp): {nterms} distinct operator words, max deviation {e2:.1e}")
# ---- U = 1: all amplitudes vanish ---------------------------------------------------------------------
print("ALL CHECKS PASSED" if max(e1, e2) < 1e-10 else "FAILED")
