"""Wilson-line chains of a term before the projection: walk from every rho index and every f leg along U's.
pattern: endpoint types and, for every U on the way, '+' if traversed from its first to its second index."""
from collections import Counter
import env  # noqa: F401


def chains(t):
    occ = {}
    for k, (i, j, p) in enumerate(t['Us']):
        occ.setdefault(i, []).append(('U', k, 0))
        occ.setdefault(j, []).append(('U', k, 1))
    for m, f in enumerate(t['fs']):
        for s, x in enumerate(f):
            occ.setdefault(x, []).append(('f', m, s))
    for r, (i, p) in enumerate(t['rs']):
        occ.setdefault(i, []).append(('r', r, 0))
    seen, out = set(), []
    starts = [('r', r, 0) for r in range(len(t['rs']))] + [('f', m, s) for m in range(len(t['fs'])) for s in range(3)]
    for st in starts:
        if st in seen:
            continue
        idx = t['rs'][st[1]][0] if st[0] == 'r' else t['fs'][st[1]][st[2]]
        frm, cur, pat, pos = st, idx, '', []
        seen.add(st)
        while True:
            nxt = [o for o in occ[cur] if o != frm]
            assert len(nxt) == 1, (cur, occ[cur], t)
            o = nxt[0]
            if o[0] in ('r', 'f'):
                seen.add(o)
                out.append((st[0], o[0], pat, tuple(pos)))
                break
            k, sl = o[1], o[2]
            i, j, p = t['Us'][k]
            pos.append(p)
            if sl == 0:
                pat += '+'
                cur, frm = j, ('U', k, 1)
            else:
                pat += '-'
                cur, frm = i, ('U', k, 0)
    return out


def norm(ch):
    """orient rho-f chains as f -> rho; rho-rho chains: both readings"""
    a, b, pat, pos = ch
    if a == 'r' and b == 'f':
        rev = ''.join('+' if c == '-' else '-' for c in reversed(pat))
        return ('f', 'r', rev, tuple(reversed(pos)))
    return ch


if __name__ == '__main__':
    import sets as S2
    import target as TG
    from parse_user import sec2_terms
    import os
    EVO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'final_check', 'input', 'Evolution_terms.tex')
    A = S2.A_sets()
    L = [t for _, t in sec2_terms(EVO)] + [t for k in ['g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III', '4rho'] for t in A[k]]
    for name, ts in (('JIMWLK target', TG.jimwlk()), ('lists', L)):
        C = Counter()
        for t in ts:
            for ch in chains(t):
                a, b, pat, pos = norm(ch)
                C[(a + '-' + b, pat)] += 1
        print(name)
        for k, v in sorted(C.items()):
            print('   %-6s %-6s %d' % (k[0], k[1], v))
