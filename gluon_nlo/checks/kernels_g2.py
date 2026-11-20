"""Group II, Row I second part, eqs (7.3)-(7.8): all six share one colour factor
C = [f f U U rho rho - ... + N_c U U rho rho].  Extract the scalar kernels mechanically
from the TeX and check that their SUM is finite on the surfaces where individual
denominators vanish."""
import re, numpy as np, texeval as te

B = te.blocks('../Final_terms_original.tex')
COL = r'\Bigg[f^{c'

def bal(s, i):
    """index after the balanced {} group starting at s[i]=='{'"""
    d = 0
    for j in range(i, len(s)):
        d += s[j] == '{'; d -= s[j] == '}'
        if d == 0:
            return j + 1
    raise ValueError

def frac(s):
    while r'\frac' in s:
        i = s.index(r'\frac'); a0 = i + 5; a1 = bal(s, a0); b1 = bal(s, a1)
        s = s[:i] + '((' + s[a0 + 1:a1 - 1] + ')/(' + s[a1 + 1:b1 - 1] + '))' + s[b1:]
    return s

def to_py(k):
    for t in [r'\notag', '\\\\', '&', r'\,', r'\left', r'\right', r'\Bigg', r'\Big', r'\bigg', r'\big']:
        k = k.replace(t, '')
    k = k.replace(r'\times', '*').replace(r'\log', 'log').replace('d_{\\perp}', 'dperp')
    k = re.sub(r'\s+', '', k)
    k = k.replace(r'\{', '').replace(r'\}', '')
    k = frac(k)
    k = k.replace('{', '(').replace('}', ')').replace('[', '(').replace(']', ')')
    P = r"([a-z]'?)-([a-z]'?)"
    k = re.sub(r'\(' + P + r'\)\\cdot\(' + P + r'\)', r'dot(\1,\2,\3,\4)', k)
    k = re.sub(r'\(' + P + r'\)\^\((\d)\)', r'pw(\1,\2,\3)', k)
    k = k.replace("'", 'p')
    k = re.sub(r'\)(?=[\(a-z])', ')*', k)          # implicit products
    k = k.replace('log*', 'log').replace('dot*', 'dot').replace('pw*', 'pw')
    k = re.sub(r'(\d)(?=[a-z\(])', r'\1*', k)
    k = k.rstrip('*')
    assert k.count('(') == k.count(')'), k
    return k

def terms(block):
    """split block into (sign*relcoef, kernel) pieces, one per colour bracket"""
    out = []
    parts = block.split(COL)[:-1]
    for p in parts:
        # the piece of text since the previous colour bracket ended
        p = p.split(r'\Bigg]', 1)[-1] if r'\rho^{e}(x)\Bigg]' in p else p
        head = p
        pre = head.split(r'\frac{1}{(2\pi)^{3}}')[0]
        for t in [r'\notag', '\\\\', '&', r'\,']:
            pre = pre.replace(t, '')
        pre = re.sub(r'\s+', '', pre)
        sign = -1 if pre.endswith('-') else 1
        coef = 1.0
        if 'd_{\\perp}}{32' in head: coef = 0.25 * 2   # d_perp -> 2          # g^4 d/(32 pi^5) relative to g^4/(8 pi^5), d kept below
        m = re.search(r"e\^\{-ik\\cdot\s*\(w-w'\)\}", head)
        ker = head[m.end():]
        out.append((sign * coef, ker, head[:m.start()][-160:]))
    return out

def dot(a, b, c, d): return np.sum((a - b) * (c - d), -1)
def pw(a, b, n): return np.sum((a - b) ** 2, -1) ** (n // 2)
log = np.log

pieces = []
for bi in range(23, 29):
    for s, k, head in terms(B[bi]):
        pieces.append((bi, s, to_py(k)))
for bi, s, k in (pieces if __name__ == '__main__' else []):
    print(f'eq block {bi}: sign*coef={s:+}  kernel = {k[:150]}...')
if __name__ == '__main__':
    np.save('g2_pieces.npy', np.array(pieces, dtype=object), allow_pickle=True)
