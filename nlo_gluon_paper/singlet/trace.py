"""Colour-singlet projection of a two-rho term and its trace form.

Projection (projectile in a colour-singlet state):
    rho^a(u) rho^b(v)  ->  delta^{ab} <rho^f(u) rho^f(v)> / (N_c^2 - 1).
The remaining colour factor is a closed network of adjoint Wilson lines U^{ij}(p) and structure constants.
It is written with adjoint generators (T^e)_{bc} = -i f^{ebc}, i.e. f^{abc} = i (T^a)_{bc}:
  * no f:      products of traces Tr[U U^dag ...];
  * two f's:   Tr[T^e A T^e B] (theta network) or Tr[A T^e] Tr[B T^e] (dumbbell network),
using  sum_d C^{ad} T^d = C^dag T^a C  for any product C of adjoint Wilson lines.
Traversing U^{ij} from j to i gives U^T = U^dag (U is real orthogonal), and (T^e)^T = -T^e.
Simplifications: U(p)U^dag(p) = 1, Tr 1 = N_c^2-1, T^eT^e = N_c, Tr T^e = 0, Tr[T^e A T^e A^dag] = N_c Tr A.
Words are tuples of letters ('U', p), ('Ud', p), ('T', label).  Coefficients are sympy numbers in N_c."""
from fractions import Fraction
import numpy as np
import sympy as sp

Nc = sp.Symbol('N_c', positive=True)


def rat(c, den=100000):
    c = complex(c)
    re = Fraction(c.real).limit_denominator(den)
    im = Fraction(c.imag).limit_denominator(den)
    assert abs(float(re) - c.real) < 1e-9 and abs(float(im) - c.imag) < 1e-9, c
    return sp.Rational(re.numerator, re.denominator) + sp.I * sp.Rational(im.numerator, im.denominator)


# ---------------------------------------------------------------- words
def tr_letter(l):
    """transpose of one letter: (sign, letter)"""
    if l[0] == 'U':
        return 1, ('Ud', l[1])
    if l[0] == 'Ud':
        return 1, ('U', l[1])
    return -1, l


def transpose_word(w):
    s, out = 1, []
    for l in reversed(w):
        si, li = tr_letter(l)
        s *= si
        out.append(li)
    return s, out


def dagger_word(w):
    """(M1..Mn)^dag for products of Wilson lines only"""
    out = []
    for l in reversed(w):
        assert l[0] in ('U', 'Ud')
        out.append(('Ud', l[1]) if l[0] == 'U' else ('U', l[1]))
    return out


def network(t):
    """singlet projection of a term: returns (extra coefficient, Us, fs) with the two rho indices identified"""
    (i1, _), (i2, _) = t['rs']
    Us, fs = [tuple(u) for u in t['Us']], [tuple(f) for f in t['fs']]
    extra = sp.Integer(1)
    if i1 == i2:
        extra = Nc**2 - 1
    else:
        ren = lambda i: i1 if i == i2 else i
        Us = [(ren(i), ren(j), p) for i, j, p in Us]
        fs = [tuple(ren(x) for x in f) for f in fs]
    return extra, Us, fs


def _cycles(edges):
    """edges: list of (i, j, word); every index appears exactly twice.  Returns (sign, list of words)."""
    occ = {}
    for k, (i, j, w) in enumerate(edges):
        occ.setdefault(i, []).append((k, 0))
        occ.setdefault(j, []).append((k, 1))
    for idx, o in occ.items():
        assert len(o) == 2, (idx, o, edges)
    used = set()
    sign, words = 1, []
    for k0 in range(len(edges)):
        if k0 in used:
            continue
        word, k, enter = [], k0, 0
        while True:
            used.add(k)
            i, j, w = edges[k]
            if enter == 0:
                word += list(w)
                nxt, exit_slot = j, 1
            else:
                s, wt = transpose_word(w)
                sign *= s
                word += wt
                nxt, exit_slot = i, 0
            cand = [(kk, sl) for kk, sl in occ[nxt] if (kk, sl) != (k, exit_slot)]
            assert len(cand) == 1, (nxt, occ[nxt])
            k2, sl2 = cand[0]
            if k2 == k0 and sl2 == 0:
                break
            k, enter = k2, sl2
        words.append(word)
    return sign, words


def _paths_between_f(Us, fs):
    """walk from every f leg along U's until another f leg.  Returns list of paths
    (f0, s0, f1, s1, word) with word the matrix product from leg (f0,s0) to leg (f1,s1), and the U-only rest."""
    occ = {}
    for k, (i, j, p) in enumerate(Us):
        occ.setdefault(i, []).append(('U', k, 0))
        occ.setdefault(j, []).append(('U', k, 1))
    for m, f in enumerate(fs):
        for s, x in enumerate(f):
            occ.setdefault(x, []).append(('f', m, s))
    usedU, done, paths = set(), set(), []
    for m, f in enumerate(fs):
        for s, x in enumerate(f):
            if (m, s) in done:
                continue
            word, cur, frm = [], x, ('f', m, s)
            while True:
                nxt = [o for o in occ[cur] if o != frm]
                assert len(nxt) == 1, (cur, occ[cur])
                o = nxt[0]
                if o[0] == 'f':
                    paths.append((m, s, o[1], o[2], word))
                    done.add((m, s))
                    done.add((o[1], o[2]))
                    break
                k, sl = o[1], o[2]
                usedU.add(k)
                i, j, p = Us[k]
                if sl == 0:
                    word.append(('U', p))
                    cur, frm = j, ('U', k, 1)
                else:
                    word.append(('Ud', p))
                    cur, frm = i, ('U', k, 0)
    rest = [u for k, u in enumerate(Us) if k not in usedU]
    return paths, rest


def gen_edge(f, s, label):
    """f with generator leg s: i (T^{f[s]})_{bc} as an edge (b, c, [T])"""
    a, b, c = f
    if s == 0:
        return (b, c, [('T', label)])
    if s == 1:
        return (c, a, [('T', label)])
    return (a, b, [('T', label)])


def to_traces(t):
    """-> (coefficient (sympy), list of words) for the singlet colour factor of term t (without the 1/(N_c^2-1))"""
    extra, Us, fs = network(t)
    coef = rat(t['c']) * Nc**t.get('ncp', 0) * extra
    if not fs:
        sgn, words = _cycles([(i, j, [('U', p)]) for i, j, p in Us])
        return coef * sgn, words
    assert len(fs) == 2, t
    paths, rest = _paths_between_f(Us, fs)
    bridges = [p for p in paths if p[0] != p[2]]
    selfs = [p for p in paths if p[0] == p[2]]
    if len(bridges) == 3:
        br = min(bridges, key=lambda p: len(p[4]))
    else:
        assert len(bridges) == 1 and len(selfs) == 2, paths
        br = bridges[0]
    m0, s0, m1, s1, C = br
    X = 'e'
    e0 = gen_edge(fs[m0], s0, X)
    e1 = gen_edge(fs[m1], s1, X)
    if C:
        e1 = (e1[0], e1[1], dagger_word(C) + [('T', X)] + list(C))
    edges = [e0, e1]
    # the other paths become matrix edges between the remaining legs (index = leg label)
    for p in paths:
        if p is br:
            continue
        m, s, mm, ss, w = p
        if not w:
            if m == mm:
                return sp.Integer(0), []          # f^{abb} = 0
            continue                              # the two legs share the index directly
        edges.append((fs[m][s], fs[mm][ss], w))
    edges += [(i, j, [('U', p)]) for i, j, p in rest]
    # the generator legs are no longer indices of the network
    sgn, words = _cycles(edges)
    return coef * sgn * sp.I**2, words


# ---------------------------------------------------------------- simplification and canonical form
def _cancel(word):
    w = list(word)
    changed = True
    while changed and w:
        changed = False
        n = len(w)
        for k in range(n):
            a, b = w[k], w[(k + 1) % n]
            if n >= 2 and a[0] in ('U', 'Ud') and b[0] in ('U', 'Ud') and a[1] == b[1] and a[0] != b[0]:
                if n == 2:
                    w = []
                else:
                    idx = sorted({k, (k + 1) % n})
                    w = [l for q, l in enumerate(w) if q not in idx]
                changed = True
                break
    return w


def simplify(coef, words):
    """returns (coef, words) or (0, None)"""
    if coef == 0:
        return sp.Integer(0), None
    out = []
    for w in words:
        w = _cancel(w)
        out.append(w)
    # T^e T^e adjacent in one word -> N_c
    res = []
    for w in out:
        ts = [k for k, l in enumerate(w) if l[0] == 'T']
        n = len(w)
        if len(ts) == 2 and ((ts[1] - ts[0]) == 1 or (ts[0] == 0 and ts[1] == n - 1)):
            coef *= Nc
            w = [l for l in w if l[0] != 'T']
            w = _cancel(w)
        elif len(ts) == 2:
            # Tr[T A T A^dag] = N_c Tr A   (A a product of Wilson lines)
            k0, k1 = ts
            A, B = w[k0 + 1:k1], w[k1 + 1:] + w[:k0]
            if A and dagger_word(A) == B:
                coef *= Nc
                w = _cancel(list(A))
        res.append(w)
    # Tr[T^e] = 0 (a word that is just T)
    for w in res:
        if len(w) == 1 and w[0][0] == 'T':
            return sp.Integer(0), None
    final = []
    for w in res:
        if not w:
            coef *= (Nc**2 - 1)
        else:
            final.append(w)
    return sp.expand(coef), final


def _rotations(w):
    return [w[k:] + w[:k] for k in range(len(w))]


def canon_word(w):
    """minimal representative under rotation and reflection (transpose); returns (sign, word)"""
    best = None
    for s, ww in ((1, list(w)), transpose_word(w)):
        for r in _rotations(ww):
            key = tuple(r)
            if best is None or key < best[1]:
                best = (s, key)
    return best


def canon(coef, words):
    """canonical (coef, tuple of words): words sorted; a generator label shared by two words is fine (label 'e')"""
    if coef == 0:
        return sp.Integer(0), ()
    ws = []
    for w in words:
        s, k = canon_word(w)
        coef *= s
        ws.append(k)
    return sp.expand(coef), tuple(sorted(ws))


def singlet(t):
    c, ws = to_traces(t)
    c, ws = simplify(c, ws)
    if ws is None:
        return sp.Integer(0), ()
    return canon(c, ws)


# ---------------------------------------------------------------- numerics
def word_matrix(w, U, Tg, e):
    M = np.eye(U[next(iter(U))].shape[0], dtype=complex)
    for l in w:
        if l[0] == 'U':
            M = M @ U[l[1]]
        elif l[0] == 'Ud':
            M = M @ U[l[1]].T
        else:
            M = M @ Tg[e]
    return M


def eval_words(coef, words, U, f, nc=3):
    """numerical value of coef * prod Tr[word], summed over the generator label (if any)"""
    Tg = [-1j * f[e] for e in range(f.shape[0])]
    c = complex(coef.subs(Nc, nc))
    hasT = any(l[0] == 'T' for w in words for l in w)
    if not hasT:
        v = c
        for w in words:
            v *= np.trace(word_matrix(w, U, Tg, 0))
        return v
    tot = 0
    for e in range(len(Tg)):
        v = c
        for w in words:
            v *= np.trace(word_matrix(w, U, Tg, e))
        tot += v
    return tot
