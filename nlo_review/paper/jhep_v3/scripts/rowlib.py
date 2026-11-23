# rowlib.py -- helpers to transcribe the equations of Sec. 8 of
# nlo_review/final_explicit/Final_terms_explicit_2_fixed.tex into the paper mechanically.
import re, os
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'final_explicit', 'Final_terms_explicit_2_fixed.tex')

def load_sec8():
    doc = open(SRC).read()
    F1 = doc[doc.index(r'\label{sec:F1}'):doc.index(r'\label{sec:F2}')]
    F2 = doc[doc.index(r'\label{sec:F2}'):doc.index(r'\label{sec:F3}')]
    return F1, F2

def aligns(T):
    return re.findall(r'\\begin\{align\}(.*?)\\end\{align\}', T, flags=re.S)

def split_equations(body):
    """split an align body into equations: rows ending with \\ not preceded by \notag end an equation"""
    rows = re.split(r'\\\\', body)
    eqs, cur = [], []
    for i, r in enumerate(rows):
        cur.append(r)
        last = (i == len(rows) - 1)
        if last or not re.search(r'\\notag\s*$', r):
            eqs.append(' '.join(cur)); cur = []
    return [e for e in eqs if e.strip()]

def flat(s):
    """one-line form: remove row/alignment/label markup"""
    lab = re.findall(r'\\label\{([^}]*)\}', s)
    s = re.sub(r'\\label\{[^}]*\}', '', s)
    s = s.replace(r'\notag', ' ').replace('&', ' ')
    s = re.sub(r'\\qquad|\\quad', ' ', s)
    s = s.replace(r'^{\prime}', "'")
    s = re.sub(r'\s+', ' ', s).strip()
    return s, (lab[0] if lab else None)

def match_brace(s, i):
    """s[i]=='{' (not escaped) -> index of matching '}'"""
    d = 0; j = i
    while j < len(s):
        c = s[j]
        if c == '\\': j += 2; continue
        if c == '{': d += 1
        elif c == '}':
            d -= 1
            if d == 0: return j
        j += 1
    raise ValueError('unbalanced')

def brace_pairs(s):
    st = []; op = {}; i = 0
    while i < len(s):
        c = s[i]
        if c == '\\': i += 2; continue
        if c == '{': st.append(i)
        elif c == '}': q = st.pop(); op[i] = q
        i += 1
    assert not st
    return op

def strip_wrappers(s, dbg=None):
    r"""remove redundant TeX groups {...} that are not arguments (left over from colour markup),
    including the braces of \cdot{(..)}; returns (string, number removed)"""
    op = brace_pairs(s)
    def is_arg(i):
        prev = s[:i].rstrip()
        if prev.endswith(('^', '_')): return True
        m = re.search(r'\\([A-Za-z]+)$', prev)
        if m is not None: return m.group(1) not in ('cdot', 'times')
        if prev.endswith('}') and not prev.endswith(r'\}'):
            q = op[len(prev) - 1]
            return s[:q].rstrip().endswith(r'\frac')   # second argument of \frac
        return False
    n = 0; out = []; i = 0; closers = set()
    while i < len(s):
        c = s[i]
        if c == '\\':
            out.append(s[i:i+2]); i += 2; continue
        if c == '{' and (not is_arg(i)) and s[i+1:i+4] != r'\rm':
            j = match_brace(s, i); closers.add(j); n += 1
            if dbg is not None: dbg.append((s[:i].rstrip()[-10:], s[i:i+30]))
            i += 1; continue
        if c == '}' and i in closers:
            i += 1; continue
        out.append(c); i += 1
    return ''.join(out), n

def norm(s):
    """normal form for comparison: no spacing/row markup, no braces, unified macros"""
    s, _ = flat(s)
    s = s.replace(r'\Kv', r'\mathcal{K}').replace(r'\Tr', r'{\rm Tr}')
    s = re.sub(r'\s+', '', s)
    s = s.replace('{', '').replace('}', '')
    return s

def top_eq_split(s):
    """split 'LHS = RHS' at the first '=' at brace depth 0"""
    d = 0; i = 0
    while i < len(s):
        c = s[i]
        if c == '\\': i += 2; continue
        if c == '{': d += 1
        elif c == '}': d -= 1
        elif c == '=' and d == 0: return s[:i].strip(), s[i+1:].strip()
        i += 1
    raise ValueError('no =')

def block(s, opener, closer, nth=0, start=0):
    """nth balanced opener...closer block in s (tokens opener/closer counted), returns (i,j) with s[i:j] the block incl. delimiters"""
    pos = start
    for _ in range(nth + 1):
        i = s.index(opener, pos); d = 0; j = i
        while True:
            a = s.find(opener, j); b = s.find(closer, j)
            if b < 0: raise ValueError('unbalanced block')
            if 0 <= a < b: d += 1; j = a + len(opener)
            else:
                d -= 1; j = b + len(closer)
                if d == 0: break
        pos = j
    return i, j
