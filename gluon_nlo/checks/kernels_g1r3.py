"""Group I row III, eqs (6.4)-(6.7): common colour factor
   i f^{abd} [U^{ab'}(x') - U^{ab'}(w')] rho^{b'}(x') x [U^{de}(x)U^{bc}(z) - U^{de}(x)U^{bc}(y)] rho^e(x) rho^c(y).
Extract the scalar kernels mechanically, sum them (with their g^4/(n pi^5) coefficients), and probe
(i) the surface (x-z)^2=(x-w)^2, (ii) collinear x,w,z (A0), (iii) the imaginary parts of logs of
negative arguments (taken on the principal branch)."""
import re, numpy as np, texeval as te
from kernels_g2 import frac

B = te.blocks('../Final_terms_original.tex')
COEF = re.compile(r"([+\-]?)\s*(?:\\,)?\s*\\frac\{ig\^\{4\}\s*(?:\\,)?\s*f\^\{abd\}\}\{(\d+)\\pi\^\{5\}\}")

def to_py(k):
    for t in [r'\notag', '\\\\', '&', r'\,', r'\left', r'\right', r'\Bigg', r'\Big', r'\bigg', r'\big']:
        k = k.replace(t, '')
    k = k.replace(r'\times', '*').replace(r'\log', 'log').replace(r'\mathcal{A}_{0}', 'A0(x,z,w)')
    k = re.sub(r'\s+', '', k)
    k = frac(k)
    k = k.replace('{', '(').replace('}', ')').replace('[', '(').replace(']', ')')
    P = r"([a-z]'?)-([a-z]'?)"
    k = re.sub(r'\(' + P + r'\)\\cdot\(' + P + r'\)', r'dot(\1,\2,\3,\4)', k)
    k = re.sub(r'\(' + P + r'\)\^\((\d)\)', r'pw(\1,\2,\3)', k)
    k = re.sub(r'\)\^\((\d)\)', r')**\1', k)
    k = k.replace("'", 'p')
    k = re.sub(r'\)(?=[\(a-zA-Z])', ')*', k)
    k = k.replace('A0*(', 'A0(')
    for f in ('log', 'dot', 'pw', 'A0'):
        k = k.replace(f + '*', f)
    k = re.sub(r'(\d)(?=[a-z\(])', r'\1*', k).rstrip('*')
    k = k.replace('A0*(', 'A0(')
    assert k.count('(') == k.count(')'), k
    return k

pieces = []
for bi in (13, 14, 15, 16):
    blk = B[bi]
    head = blk[:blk.index(r'\int_{w')]
    bsign = -1 if re.search(r'=\s*-', head.replace(r'\,', '')) else 1
    ms = list(COEF.finditer(blk))
    for j, m in enumerate(ms):
        end = ms[j + 1].start() if j + 1 < len(ms) else len(blk)
        seg = blk[m.end():end]
        seg = seg.split(r"\int_{x,y,z,x'}", 1)[1]
        seg = re.split(r"\\(?:left|Big)\[U\^\{ab'\}", seg)[0]
        s = (-1 if m.group(1) == '-' else 1) * bsign * 16.0 / float(m.group(2))   # relative to g^4/(16 pi^5)
        pieces.append((bi, s, to_py(seg)))

def dot(a, b, c, d): return np.sum((a - b) * (c - d), -1)
def pw(a, b, n): return np.sum((a - b) ** 2, -1) ** (n // 2)
log = lambda z: np.log(np.asarray(z, dtype=complex))        # principal branch: exposes log(negative)
def A0(x, z, w):
    a, b = x - w, x - z
    S = np.sqrt(dot(x, w, x, w) * dot(x, z, x, z) - dot(x, w, x, z) ** 2)
    return np.arccos(dot(x, w, x, w - (x - z) + x - x) / np.sqrt(dot(x, w, x, w) * dot(w, z, w, z))) / S if False else \
        np.arccos((dot(x, w, x, w) - dot(x, w, x, z)) / np.sqrt(dot(x, w, x, w) * dot(w, z, w, z))) / S
def A0_arctan(x, z, w):
    S = np.sqrt(dot(x, w, x, w) * dot(x, z, x, z) - dot(x, w, x, z) ** 2)
    return (np.arctan((dot(x, z, x, z) - dot(x, w, x, z)) / S) + np.arctan(dot(x, w, x, z) / S)) / S

fns = [(bi, s, k, eval('lambda x,y,z,xp,w,wp: ' + k)) for bi, s, k in pieces]
def total(P, parts=False):
    v = [s * f(*P) for bi, s, k, f in fns]
    return v if parts else sum(v)

if __name__ == '__main__':
    for bi, s, k, f in fns:
        print(f'block {bi}: coef {s:+.0f}  kernel {k[:110]}...')
