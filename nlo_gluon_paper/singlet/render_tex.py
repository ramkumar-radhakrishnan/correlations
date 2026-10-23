"""LaTeX for collected singlet forms: {kernel key: {trace key: coefficient}}."""
import sympy as sp
from trace import Nc

PL = {'x': 'x', 'y': 'y', 'w': 'w', "w'": "w'", 'z': 'z'}
RANK = {'x': 0, 'w': 1, 'z': 2, 'y': 3, "w'": 4}


def kernel_tex(kk):
    lo, zf = [], []
    for f in kk:
        if f[0] == 'inv2':
            (a, b), = f[1:]
            if b == 'z':
                zf.append('K(%s,%s,z)' % (PL[a], PL[a]))
            else:
                lo.append(r'\frac{1}{(%s-%s)^2}' % (PL[a], PL[b]))
        else:
            (a, b), (c, d) = f[1:]
            if b == 'z' and d == 'z':
                zf.append('K(%s,%s,z)' % (PL[a], PL[c]))
            else:
                lo.append(r'\mathcal K(%s-%s)\cdot\mathcal K(%s-%s)' % (PL[a], PL[b], PL[c], PL[d]))
    return r'\,'.join(lo + zf)


def letter_tex(l):
    if l[0] == 'U':
        return 'U(%s)' % PL[l[1]]
    if l[0] == 'Ud':
        return r'U^\dagger(%s)' % PL[l[1]]
    return 'T^e'


def word_tex(w):
    if len(w) == 2 and {w[0][0], w[1][0]} == {'U', 'Ud'}:
        p, q = sorted([w[0][1], w[1][1]], key=lambda s: RANK[s])
        return r'{\rm Tr}\{U^\dagger(%s)U(%s)\}' % (PL[p], PL[q])
    return r'{\rm Tr}\{' + ''.join(letter_tex(l) for l in w) + r'\}'


def coef_tex(c, first):
    c = sp.factor(c)
    neg = (c.could_extract_minus_sign())
    a = -c if neg else c
    s = sp.latex(a).replace('N_{c}', 'N_c').replace(r'\left(N_c - 1\right) \left(N_c + 1\right)', r'\left(N_c^{2} - 1\right)')
    if a == 1:
        s = ''
    sign = '-' if neg else ('' if first else '+')
    return sign + s


def mono_width(ws):
    return sum(len(w) + 2 for w in ws) if ws else 1


def block_tex(coll, maxw=17, indent=r'\quad'):
    """align body (no environment) for one collected dict"""
    lines = []
    for n, (kk, d) in enumerate(coll.items()):
        head = ('' if n == 0 else '+') + kernel_tex(kk) + r'\Big['
        parts, cur, width, first = [], head, 0, True
        for ws, c in d.items():
            if not ws:
                t = coef_tex(c, first) or '1'
                if t in ('+', '-'):
                    t += '1'
            else:
                t = coef_tex(c, first) + r'\,'.join(word_tex(w) for w in ws)
            wdt = mono_width(ws)
            if width + wdt > maxw and not first:
                parts.append(cur)
                cur, width = indent + r'\quad', 0
            cur += t
            width += wdt
            first = False
        cur += r'\Big]'
        parts.append(cur)
        lines.append(('\\\\\n&' + indent).join(parts) if len(parts) > 1 else parts[0])
    return '&' + '\\\\\n&'.join(lines)


# ---------------------------------------------------------------- dipole factoring
CANDS = ['x', 'y', 'z']


def _dip_key(ws):
    """('D', frozenset) for Tr{U^dag(p) U(q)}, ('C',) for the empty monomial, None otherwise"""
    if not ws:
        return ('C',)
    if len(ws) == 1 and len(ws[0]) == 2 and {ws[0][0][0], ws[0][1][0]} == {'U', 'Ud'}:
        return ('D', frozenset([ws[0][0][1], ws[0][1][1]]))
    return None


def factor_dipoles(d):
    """d: {trace key: coef}.  Returns (list of (alpha, a, b) for alpha Tr{(U^dag(a)-U^dag(w))(U(b)-U(w'))},
    remaining dict).  D(a,a) is taken from the constant (Tr 1 = N_c^2-1)."""
    dip = {}
    rest = {}
    for ws, c in d.items():
        k = _dip_key(ws)
        if k is None:
            rest[ws] = c
        else:
            dip[k] = dip.get(k, 0) + c
    def get(p, q):
        if p == q:
            return sp.expand(sp.cancel(dip.get(('C',), 0) / (Nc**2 - 1)))
        return dip.get(('D', frozenset([p, q])), 0)
    def add(p, q, v):
        if p == q:
            dip[('C',)] = sp.expand(dip.get(('C',), 0) + v * (Nc**2 - 1))
        else:
            k = ('D', frozenset([p, q]))
            dip[k] = sp.expand(dip.get(k, 0) + v)
    facs = []
    changed = True
    while changed:
        changed = False
        for a in CANDS:
            for b in CANDS:
                al = sp.expand(get(a, b))
                if al == 0:
                    continue
                if sp.expand(get(a, "w'") + al) == 0 and sp.expand(get('w', b) + al) == 0:
                    facs.append((al, a, b))
                    add(a, b, -al)
                    add(a, "w'", al)
                    add('w', b, al)
                    add('w', "w'", -al)
                    changed = True
    remaining = dict(rest)
    for k, v in dip.items():
        v = sp.expand(v)
        if v == 0:
            continue
        if k == ('C',):
            remaining[()] = v
        else:
            p, q = sorted(k[1], key=lambda s: RANK[s])
            remaining[((('U', p), ('Ud', q)),)] = v
    return facs, remaining


def block_tex_f(coll, maxw=17, indent=r'\quad'):
    """like block_tex, with the dipoles factored into Tr{(U^dag(a)-U^dag(w))(U(b)-U(w'))}"""
    lines = []
    for n, (kk, d) in enumerate(coll.items()):
        facs, rem = factor_dipoles(d)
        items = []
        for al, a, b in facs:
            items.append((al, r"{\rm Tr}\{(U^\dagger(%s)-U^\dagger(w))(U(%s)-U(w'))\}" % (PL[a], PL[b]), 7))
        for ws, c in sorted(rem.items(), key=lambda it: (len(it[0]), sum(len(w) for w in it[0]), str(it[0]))):
            if not ws:
                items.append((c, '', 1))
            else:
                items.append((c, r'\,'.join(word_tex(w) for w in ws), mono_width(ws)))
        head = ('' if n == 0 else '+') + kernel_tex(kk) + r'\Big['
        parts, cur, width, first = [], head, 6, True
        for c, tx, wdt in items:
            t = coef_tex(c, first)
            if not tx:
                t = (t if t not in ('', '+', '-') else t + '1')
            t += tx
            if width + wdt > maxw and not first:
                parts.append(cur)
                cur, width = indent + r'\quad', 0
            cur += t
            width += wdt
            first = False
        cur += r'\Big]'
        parts.append(cur)
        lines.append(('\\\\\n&' + indent).join(parts) if len(parts) > 1 else parts[0])
    return '&' + '\\\\\n&'.join(lines)


def unfactor(facs, rem):
    """expand the factored form back to {('D', pair) | ('C',) | trace key: coef} for checking"""
    out = {}
    def add(k, v):
        out[k] = sp.expand(out.get(k, 0) + v)
    def D(p, q, v):
        if p == q:
            add(('C',), v * (Nc**2 - 1))
        else:
            add(('D', frozenset([p, q])), v)
    for al, a, b in facs:
        D(a, b, al); D(a, "w'", -al); D('w', b, -al); D('w', "w'", al)
    for ws, c in rem.items():
        k = _dip_key(ws)
        add(k if k is not None else ws, c)
    return {k: v for k, v in out.items() if sp.expand(v) != 0}


def check_factoring(coll):
    """True if every kernel block expands back to the original monomials"""
    for kk, d in coll.items():
        facs, rem = factor_dipoles(d)
        orig = {}
        for ws, c in d.items():
            k = _dip_key(ws)
            k = k if k is not None else ws
            orig[k] = sp.expand(orig.get(k, 0) + c)
        orig = {k: v for k, v in orig.items() if v != 0}
        back = unfactor(facs, rem)
        if set(orig) != set(back) or any(sp.simplify(orig[k] - back[k]) != 0 for k in orig):
            return False
    return True
