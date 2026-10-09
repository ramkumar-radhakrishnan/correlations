"""Reordering of three-rho terms into symmetric part + two-rho terms.

Colour factors are lists of tuples
    ('U', i, j, pos)   adjoint Wilson line U^{ij}(pos)
    ('f', i, j, k)     structure constant f^{ijk}
    ('d', i, j)        Kronecker delta
    ('rho', i, pos)    valence charge rho^i(pos)
All colour indices are summed (each appears exactly twice).
"""
import itertools
import sympy as sp

Nc = sp.Symbol('N_c')
I = sp.I

CANON = ['a', 'b', 'c', 'd', 'e', 'h', 'l', 'm', 'n', 'r', 's', 't', 'o', 'q']
FRESH = ['h', 'l', 'm', 'n', 'r', 's', 't', 'o', 'q', 'g_1', 'g_2', 'g_3']


# ---------------------------------------------------------------- helpers
def indices_of(fac):
    if fac[0] == 'U':
        return [fac[1], fac[2]]
    if fac[0] == 'f':
        return [fac[1], fac[2], fac[3]]
    if fac[0] == 'd':
        return [fac[1], fac[2]]
    if fac[0] == 'rho':
        return [fac[1]]
    raise ValueError(fac)


def all_indices(facs):
    out = []
    for f in facs:
        out += indices_of(f)
    return out


def check_contracted(facs, where=''):
    from collections import Counter
    c = Counter(all_indices(facs))
    bad = {k: v for k, v in c.items() if v != 2}
    if bad:
        raise ValueError(f'index occurrence error {bad} in {where}: {facs}')


def fresh(used, k=1):
    out = []
    for n in FRESH:
        if n not in used and n not in out:
            out.append(n)
            if len(out) == k:
                return out
    raise RuntimeError('no fresh index')


def rename_index(fac, old, new):
    if fac[0] == 'U':
        return ('U', new if fac[1] == old else fac[1], new if fac[2] == old else fac[2], fac[3])
    if fac[0] == 'f':
        return ('f',) + tuple(new if x == old else x for x in fac[1:4])
    if fac[0] == 'd':
        return ('d', new if fac[1] == old else fac[1], new if fac[2] == old else fac[2])
    if fac[0] == 'rho':
        return ('rho', new if fac[1] == old else fac[1], fac[2])


def subst_pos(fac, old, new):
    if fac[0] == 'U':
        return fac[:3] + ((new if fac[3] == old else fac[3]),)
    if fac[0] == 'rho':
        return ('rho', fac[1], new if fac[2] == old else fac[2])
    return fac


def perm_sign(seq_from, seq_to):
    """sign of the permutation taking seq_from to seq_to (both length 3)."""
    p = [seq_from.index(x) for x in seq_to]
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def f_bring_front(fac, i, j):
    """write f^{...} as sign * f^{i j k}; returns (sign, k)."""
    idx = list(fac[1:4])
    k = [x for x in idx if x not in (i, j)]
    assert len(k) == 1, (fac, i, j)
    k = k[0]
    return perm_sign(idx, [i, j, k]), k


# ---------------------------------------------------------- simplification
def simplify(coeff, facs):
    """Apply colour identities.  Returns (coeff, facs, rules) where rules is a
    list of short tags describing the identities used, in order."""
    facs = list(facs)
    rules = []
    changed = True
    while changed:
        changed = False
        # (1) Kronecker deltas
        for n, fac in enumerate(facs):
            if fac[0] != 'd':
                continue
            i, j = fac[1], fac[2]
            rest = facs[:n] + facs[n + 1:]
            if i == j:
                coeff = coeff * (Nc**2 - 1)
                facs = rest
                rules.append('trd')
            else:
                if j in all_indices(rest):
                    facs = [rename_index(x, j, i) for x in rest]
                else:
                    facs = [rename_index(x, i, j) for x in rest]
                rules.append('delta')
            changed = True
            break
        if changed:
            continue
        # (2) f with a repeated index vanishes
        for fac in facs:
            if fac[0] == 'f' and len(set(fac[1:4])) < 3:
                return sp.Integer(0), [], rules + ['f_ii=0']
        # (3) two rho's at the same point contracted with one f vanish
        rhos = [x for x in facs if x[0] == 'rho']
        for fac in facs:
            if fac[0] != 'f':
                continue
            for r1, r2 in itertools.combinations(rhos, 2):
                if r1[2] == r2[2] and r1[1] in fac[1:4] and r2[1] in fac[1:4]:
                    return sp.Integer(0), [], rules + ['f_rhorho=0']
        # (4) f^{ijk} U^{ij}(p): zero only if U is symmetric
        for fac in facs:
            if fac[0] != 'f':
                continue
            for u in facs:
                if u[0] == 'U' and u[1] in fac[1:4] and u[2] in fac[1:4] and u[1] != u[2]:
                    return sp.Integer(0), [], rules + ['fU=0(sym)']
        # (5) U U at the same point sharing an index -> delta
        Us = [(n, x) for n, x in enumerate(facs) if x[0] == 'U']
        for (n1, u1), (n2, u2) in itertools.combinations(Us, 2):
            if u1[3] != u2[3]:
                continue
            shared = [s for s in (u1[1], u1[2]) if s in (u2[1], u2[2])]
            if not shared:
                continue
            s = shared[0]
            slot1 = 1 if u1[1] == s else 2
            slot2 = 1 if u2[1] == s else 2
            o1 = u1[2] if slot1 == 1 else u1[1]
            o2 = u2[2] if slot2 == 1 else u2[1]
            rest = [x for k, x in enumerate(facs) if k not in (n1, n2)]
            facs = rest + [('d', o1, o2)]
            rules.append('UU' if slot1 == slot2 else 'UU(sym)')
            changed = True
            break
        if changed:
            continue
        # (6) f U U at the same point -> f U
        for nf, fac in enumerate(facs):
            if fac[0] != 'f':
                continue
            done = False
            for i, j in itertools.permutations(fac[1:4], 2):
                ui = [(n, x) for n, x in enumerate(facs) if x[0] == 'U' and i in (x[1], x[2])]
                uj = [(n, x) for n, x in enumerate(facs) if x[0] == 'U' and j in (x[1], x[2])]
                if not ui or not uj:
                    continue
                (n1, u1), (n2, u2) = ui[0], uj[0]
                if n1 == n2 or u1[3] != u2[3]:
                    continue
                sign, k = f_bring_front(fac, i, j)
                s1 = 1 if u1[1] == i else 2
                s2 = 1 if u2[1] == j else 2
                d = u1[2] if s1 == 1 else u1[1]
                e = u2[2] if s2 == 1 else u2[1]
                used = set(all_indices(facs))
                hh = fresh(used)[0]
                rest = [x for m, x in enumerate(facs) if m not in (nf, n1, n2)]
                newU = ('U', k, hh, u1[3]) if s1 == 1 else ('U', hh, k, u1[3])
                facs = rest + [('f', d, e, hh), newU]
                coeff = coeff * sign
                rules.append('fUU' if s1 == s2 else 'fUU(sym)')
                done = True
                break
            if done:
                changed = True
                break
        if changed:
            continue
        # (7) f f with two (or three) shared indices
        Fs = [(n, x) for n, x in enumerate(facs) if x[0] == 'f']
        for (n1, f1), (n2, f2) in itertools.combinations(Fs, 2):
            shared = [s for s in f1[1:4] if s in f2[1:4]]
            if len(shared) == 3:
                sign = perm_sign(list(f1[1:4]), list(f2[1:4]))
                rest = [x for k, x in enumerate(facs) if k not in (n1, n2)]
                facs = rest
                coeff = coeff * sign * Nc * (Nc**2 - 1)
                rules.append('ff3')
                changed = True
                break
            if len(shared) == 2:
                i, j = shared
                s1, k1 = f_bring_front(f1, i, j)
                s2, k2 = f_bring_front(f2, i, j)
                rest = [x for k, x in enumerate(facs) if k not in (n1, n2)]
                facs = rest + [('d', k1, k2)]
                coeff = coeff * s1 * s2 * Nc
                rules.append('ff')
                changed = True
                break
    return sp.expand(coeff), facs, rules


def canonical_rename(facs, order_key=None):
    """rename all (dummy) indices to a, b, c, ... in order of appearance."""
    seq = []
    for f in sort_factors(facs):
        for x in indices_of(f):
            if x not in seq:
                seq.append(x)
    mapping = {old: CANON[n] for n, old in enumerate(seq)}
    tmp = [rename_all(f, {k: '#' + v for k, v in mapping.items()}) for f in facs]
    return [rename_all(f, {'#' + v: v for v in mapping.values()}) for f in tmp]


def rename_all(fac, mp):
    if fac[0] == 'U':
        return ('U', mp.get(fac[1], fac[1]), mp.get(fac[2], fac[2]), fac[3])
    if fac[0] == 'f':
        return ('f',) + tuple(mp.get(x, x) for x in fac[1:4])
    if fac[0] == 'd':
        return ('d', mp.get(fac[1], fac[1]), mp.get(fac[2], fac[2]))
    if fac[0] == 'rho':
        return ('rho', mp.get(fac[1], fac[1]), fac[2])


def sort_factors(facs):
    order = {'f': 0, 'd': 1, 'U': 2, 'rho': 3}
    return sorted(facs, key=lambda f: order[f[0]])


# ------------------------------------------------------------------ LaTeX
def tex_pos(p):
    return p


def tex_idx(i):
    return i


def tex_colour(facs, rho_mode='word', keep_order=False):
    """rho_mode: 'word' prints rho's as an ordered product, 'anti' prints the two
    rho's as an anticommutator {rho, rho}."""
    fs = facs if keep_order else sort_factors(facs)
    out = []
    rhos = []
    for f in fs:
        if f[0] == 'f':
            out.append('f^{%s%s%s}' % (f[1], f[2], f[3]))
        elif f[0] == 'd':
            out.append(r'\delta^{%s%s}' % (f[1], f[2]))
        elif f[0] == 'U':
            out.append('U^{%s%s}(%s)' % (f[1], f[2], f[3]))
        elif f[0] == 'rho':
            rhos.append(r'\rho^{%s}(%s)' % (f[1], f[2]))
    if rho_mode == 'anti' and len(rhos) == 2:
        out.append(r'\big\{%s,%s\big\}' % (rhos[0], rhos[1]))
    else:
        out += rhos
    return r'\,'.join(out) if out else '1'


def tex_coeff(c, leading=False):
    """LaTeX for a sympy coefficient; returns (sign_str, body) with body possibly ''."""
    c = sp.nsimplify(sp.expand(c))
    if c == 0:
        return '+', '0'
    # split numeric factor
    num, rest = c.as_coeff_Mul()
    if num < 0:
        sign = '-'
        num = -num
    else:
        sign = '+'
    body = ''
    if rest != 1:
        rest_tex = sp.latex(rest)
        if rest.is_Add:
            rest_tex = r'\big(%s\big)' % rest_tex
        if num == 1:
            body = rest_tex + r'\,'
        else:
            body = sp.latex(num) + r'\,' + rest_tex + r'\,'
    else:
        body = '' if num == 1 else sp.latex(num) + r'\,'
    body = body.replace(r'\frac{1}{4}', r'\tfrac14').replace(r'\frac{1}{2}', r'\tfrac12')
    return sign, body


# --------------------------------------------------------- decomposition
class Rho:
    def __init__(self, idx, pos):
        self.idx, self.pos = idx, pos


def pieces_for(form, word):
    """Return list of (coefficient, X, Y, Z, tag): piece = coeff * {[X,Y], Z}.
    form 'plain': word = [A,B,C]  ->  ABC
    form 'A{BC}': word = [A,B,C]  ->  A{B,C}
    form '{BC}A': word = [B,C,A]  ->  {B,C}A  (A is the single operator on the right)
    """
    if form == 'plain':
        A, B, C = word
        q = sp.Rational(1, 4)
        return [(q, A, B, C, 'AB'), (q, A, C, B, 'AC'), (q, B, C, A, 'BC')]
    if form == 'A{BC}':
        A, B, C = word
        h = sp.Rational(1, 2)
        return [(h, A, B, C, 'AB'), (h, A, C, B, 'AC')]
    if form == '{BC}A':
        B, C, A = word
        h = sp.Rational(1, 2)
        return [(-h, A, B, C, 'AB'), (-h, A, C, B, 'AC')]
    raise ValueError(form)


def make_piece(colour_nonrho, coeff, X, Y, Z):
    """{[X,Y],Z} with [rho^x(p), rho^y(q)] = i f^{xyh} rho^h(p) delta(p-q);
    the delta is integrated over q (the position of the right operator)."""
    used = set(all_indices(colour_nonrho)) | {X.idx, Y.idx, Z.idx}
    hh = fresh(used)[0]
    facs = list(colour_nonrho) + [('f', X.idx, Y.idx, hh), ('rho', hh, X.pos), ('rho', Z.idx, Z.pos)]
    old, new = Y.pos, X.pos
    facs = [subst_pos(f, old, new) for f in facs]
    return coeff * I, facs, (old, new)
