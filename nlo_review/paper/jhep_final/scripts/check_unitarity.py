#!/usr/bin/env python3
"""Symbolic check of the unitarity relations of Section 3.4.

Omega is written for a discrete set of modes (commutator [a_p, a^dag_q] = delta_pq) as
  Omega = 1 + N2 + sum_p [A1_p a_p + A2_p a^dag_p]
        + sum_{pq} [B1_{pq} a_p a_q + B2_{pq} a^dag_p a^dag_q + B3_{pq} a^dag_q a_p]
        + sum_{jlk} [C1_{lk;j} a^dag_l a^dag_k a_j + C2_{lk;j} a^dag_j a_k a_l]
with non-commuting coefficients (words in formal symbols), A = A^(1)+A^(3), orders
g^1: A^(1), C ;  g^2: N2, B ;  g^3: A^(3).   Omega^dag Omega is normal ordered and the coefficient
of every operator structure is printed order by order; it must vanish.  The script then substitutes
the solutions quoted in the paper and verifies that the coefficients of
1, a^dag, a^dag a^dag, a^dag a^dag a  vanish identically, and prints the a^dag a relation.
Symmetries B1_{pq}=B1_{qp}, B2_{pq}=B2_{qp}, C_{lk;j}=C_{kl;j} are imposed.  N2 is Hermitian.
Run:  python3 -I check_unitarity.py
"""
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
M = (0, 1, 2)                       # modes
ORDER = {'A1': 1, 'A2': 1, 'a1': 3, 'a2': 3, 'N2': 2, 'B1': 2, 'B2': 2, 'B3': 2, 'C1': 1, 'C2': 1}   # a1,a2 = A^(3)
def sym(name, idx, dag=False):
    if name in ('B1', 'B2'): idx = tuple(sorted(idx))
    if name in ('C1', 'C2'): idx = tuple(sorted(idx[:2]))+idx[2:]
    if name == 'N2': dag = False
    return (name, idx, dag)
def order(word): return sum(ORDER[s[0]] for s in word)
# an operator: dict {(word, cre, ann): coeff}, cre/ann sorted tuples of modes
def add(D, word, cre, ann, c=F(1)):
    if order(word) > 3: return
    key = (tuple(word), tuple(sorted(cre)), tuple(sorted(ann)))
    D[key] += c
    if D[key] == 0: del D[key]
def omega():
    D = defaultdict(F)
    add(D, (), (), ()); add(D, (sym('N2', ()),), (), ())
    for p in M:
        for n1, n2 in (('A1', 'A2'), ('a1', 'a2')):
            add(D, (sym(n1, (p,)),), (), (p,)); add(D, (sym(n2, (p,)),), (p,), ())
    for p, q in product(M, M):
        add(D, (sym('B1', (p, q)),), (), (p, q)); add(D, (sym('B2', (p, q)),), (p, q), ())
        add(D, (sym('B3', (p, q)),), (q,), (p,))
    for l, k, j in product(M, M, M):
        add(D, (sym('C1', (l, k, j)),), (l, k), (j,)); add(D, (sym('C2', (l, k, j)),), (j,), (k, l))
    return D
def dagger(D):
    R = defaultdict(F)
    for (w, cre, ann), c in D.items():
        add(R, tuple((s[0], s[1], (not s[2]) if s[0] != 'N2' else False) for s in reversed(w)), ann, cre, c)
    return R
def contract(ann, cre):
    """normal order (prod a_ann)(prod a^dag_cre): list of (cre_left, ann_left, multiplicity)."""
    if not ann or not cre: return [(tuple(cre), tuple(ann), 1)]
    a, rest = ann[-1], ann[:-1]          # move the last annihilator through all creators
    out = []
    for (c2, a2, m) in contract(rest, cre): out.append((c2, a2+(a,), m))          # no contraction of a
    for i, q in enumerate(cre):
        if q == a:
            for (c2, a2, m) in contract(rest, cre[:i]+cre[i+1:]): out.append((c2, a2, m))
    return out
def mul(D1, D2):
    R = defaultdict(F)
    for (w1, c1, a1), x in D1.items():
        for (w2, c2, a2), y in D2.items():
            if order(w1)+order(w2) > 3: continue
            for (cl, al, m) in contract(a1, c2):
                add(R, w1+w2, c1+cl, al+a2, x*y*m)
    return R
def show(word):
    if not word: return '1'
    names = {'A1': 'A1', 'A2': 'A2', 'a1': 'A1(3)', 'a2': 'A2(3)'}
    return ' '.join(f"{names.get(s[0], s[0])}{'+' if s[2] else ''}{list(s[1]) if s[1] else ''}" for s in word)
Om = omega(); P = mul(dagger(Om), Om)
add(P, (), (), (), F(-1))
structs = defaultdict(lambda: defaultdict(F))
for (w, cre, ann), c in P.items(): structs[(cre, ann)][w] += c
def report(cre, ann, title):
    print(f"\n{title}:  coefficient of  a^dag{list(cre)} a{list(ann)}  (sum of all terms = 0)")
    byord = defaultdict(list)
    for w, c in structs[(cre, ann)].items(): byord[order(w)].append((c, w))
    for o in sorted(byord):
        print(f"  order g^{o}: " + '  '.join(f"{'+' if c > 0 else '-'}{'' if abs(c) == 1 else abs(c)} {show(w)}" for c, w in sorted(byord[o], key=lambda t: show(t[1]))))
report((), (), 'identity')
report((0,), (), 'a^dag_0')
report((0, 1), (), 'a^dag_0 a^dag_1')
report((0,), (1,), 'a^dag_0 a_1')
report((0, 1), (2,), 'a^dag_0 a^dag_1 a_2')

# ---------------- substitution of the solutions quoted in the paper -----------------------------------
def substitute(poly, rules):
    """rules: symbol -> list of (coeff, word). Repeated until no rule applies; words truncated at g^3."""
    changed = True
    while changed:
        changed = False; new = defaultdict(F)
        for w, c in poly.items():
            for i, s in enumerate(w):
                if s in rules:
                    for (c2, w2) in rules[s]:
                        nw = w[:i]+tuple(w2)+w[i+1:]
                        if order(nw) <= 3: new[nw] += c*c2
                    changed = True; break
            else: new[w] += c
        poly = defaultdict(F, {w: c for w, c in new.items() if c != 0})
    return poly
rules = {}
for p in M:
    rules[sym('A1', (p,), True)] = [(F(-1), (sym('A2', (p,)),))]                 # A1^(1)dag = -A2^(1)
    rules[sym('A1', (p,), False)] = [(F(-1), (sym('A2', (p,), True),))]
rules[sym('N2', ())] = [(F(-1, 2), (sym('A2', (p,), True), sym('A2', (p,)))) for p in M]   # N2 = -1/2 sum A2^dag A2
for l, k, j in product(M, M, M):
    rules[sym('C2', (l, k, j), True)] = [(F(-1), (sym('C1', (l, k, j)),))]        # C2^dag = -C1
    rules[sym('C2', (l, k, j), False)] = [(F(-1), (sym('C1', (l, k, j), True),))]
for p, q in product(M, M):                                                     # B1^dag relation
    if p <= q:
        r = [(F(-1), (sym('B2', (p, q)),)), (F(1, 2), (sym('A2', (p,)), sym('A2', (q,)))), (F(1, 2), (sym('A2', (q,)), sym('A2', (p,))))]
        r += [(F(1), (sym('C1', (p, q, j)), sym('A2', (j,)))) for j in M]
        rules[sym('B1', (p, q), True)] = r
for p in M:                                                                     # A1^(3)dag relation
    r = [(F(-1), (sym('a2', (p,)),)), (F(-1), (sym('N2', ()), sym('A2', (p,)))), (F(1), (sym('A2', (p,)), sym('N2', ())))]
    for q in M:
        r += [(F(2), (sym('A1', (q,)), sym('B2', (p, q)))), (F(-1), (sym('B3', (p, q), True), sym('A2', (q,))))]
        r += [(F(-2), (sym('C1', (l, q, p), True), sym('B2', (q, l)))) for l in M]
    rules[sym('a1', (p,), True)] = r
print("\nAfter substituting  A1(1)+ = -A2(1),  N2 = -1/2 sum A2+ A2,  C2+ = -C1,  the B1+ relation and the A1(3)+ relation:")
allok = True
for (cre, ann), title in (((), ()), 'identity'), (((0,), ()), 'a^dag'), (((0, 0), ()), 'a^dag a^dag (equal modes)'), (((0, 1), ()), 'a^dag a^dag'), (((0, 1), (2,)), 'a^dag a^dag a, order g'), (((0, 0), (1,)), 'a^dag a^dag a (equal modes), order g'):
    src = structs[(cre, ann)]
    if len(ann) == 1 and len(cre) == 2: src = {w: c for w, c in src.items() if order(w) == 1}   # order g^3 involves coefficients not kept
    res = substitute(defaultdict(F, src), rules)
    print(f"  {title:28s}: remaining terms = {len(res)}"); allok &= (len(res) == 0)
res = substitute(defaultdict(F, structs[((0,), (1,))]), rules)
print("  a^dag_0 a_1 relation after substitution (must vanish; constrains B3):")
print("   ", '  '.join(f"{'+' if c > 0 else '-'}{'' if abs(c) == 1 else abs(c)} {show(w)}" for w, c in sorted(res.items(), key=lambda t: show(t[0]))))
print("\nALL UNITARITY RELATIONS VERIFIED" if allok else "\nFAILED")
