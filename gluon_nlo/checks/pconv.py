"""TeX -> python for kernels that carry a longitudinal integral (Group I row V, Group III).
Handles free vector indices (summed explicitly over 2 transverse components), Kronecker deltas,
plus-distributions (replaced by symbols PL0, PL1, PLQ, PLQK, PLG whose subtraction is done at
integration time), and k-dependent phases."""
import re
from kernels_g2 import frac

PLUS = [  # (regex on cleaned TeX, symbol)
    (r'\[\\frac\{k\^\{\+\}\}\{p\^\{\+\}\(p\^\{\+\}\+k\^\{\+\}\)\}\]_\{\+\}', 'PLG'),
    (r'\[\\frac\{k\^\{\+\}\}\{k\^\{\+\}-p\^\{\+\}\}\]_\{\+\}', 'PL1'),
    (r'\[\\frac\{k\^\{\+\}\}\{p\^\{\+\}\}\]_\{\+\}', 'PL0'),
    (r'\[\\frac\{k\^\{\+\}\}\{q\^\{\+\}\}\]_\{\+\}', 'PLQK'),
    (r'\[\\frac\{1\}\{q\^\{\+\}\}\]_\{\+\}', 'PLQ'),
]
IDX = {"i": "ii", "j": "jj", "k": "kk", "m": "mm", "k'": "kq", "m'": "mq"}

def clean(k):
    for t in [r'\notag', '\\\\', '&', r'\,', r'\left', r'\right', r'\Bigg', r'\Big', r'\bigg', r'\big', r'\color{red}']:
        k = k.replace(t, '')
    k = k.replace(r'w^\prime', "w'").replace(r'w^{\prime}', "w'")
    return re.sub(r'\s+', '', k)

def to_py(k):
    k = clean(k)
    for pat, sym in PLUS:
        k = re.sub(pat, sym, k)
    P = r"([a-z]'?)-([a-z]'?)"
    # phases
    k = re.sub(r"e\^\{-ik\\cdot\(" + P + r"\)\\frac\{p\^\{\+\}\}\{k\^\{\+\}\}\}", r"EXP(\1,\2,p/K)", k)
    k = re.sub(r"e\^\{-ik\\cdot\(" + P + r"\)\\frac\{q\^\{\+\}\}\{k\^\{\+\}\}\}", r"EXP(\1,\2,q/K)", k)
    k = re.sub(r"e\^\{-ik\\cdot\\frac\{p\^\{\+\}\}\{k\^\{\+\}\}\(" + P + r"\)\}", r"EXP(\1,\2,p/K)", k)
    k = re.sub(r"e\^\{-ik\\cdot\(" + P + r"\)\}", r"EXP(\1,\2,1)", k)
    k = re.sub(r"e\^\{-ik\(" + P + r"\)\}", r"EXP(\1,\2,1)", k)
    assert 'e^' not in k, k[k.index('e^'):k.index('e^') + 60]
    k = k.replace('k^{+2}', 'K**2').replace('k^{+2}', 'K**2').replace(r'k^{+\,2}', 'K**2').replace('k^{+}', 'K')
    k = k.replace('p^{+}', 'p').replace('q^{+}', 'q').replace(r'd_{\perp}', 'DPERP').replace(r'\times', '*')
    # deltas
    k = re.sub(r"\\delta[_^]\{([a-z]'?)([a-z]'?)\}", lambda m: f"D({IDX[m.group(1)]},{IDX[m.group(2)]})", k)
    k = frac(k)
    k = re.sub(r'\(' + P + r'\)\\cdot\(' + P + r'\)', r'DOT(\1,\2,\3,\4)', k)
    k = re.sub(r'\(' + P + r'\)\^\{2\}', r'PW(\1,\2,2)', k)
    k = re.sub(r'\(' + P + r"\)\^\{([a-z]'?)\}", lambda m: f"V({m.group(1)},{m.group(2)},{IDX[m.group(3)]})", k)
    k = k.replace(r'\{', '(').replace(r'\}', ')')
    k = k.replace('{', '(').replace('}', ')').replace('[', '(').replace(']', ')')
    k = k.replace("'", 'p')
    k = re.sub(r'\)(?=[\(a-zA-Z])', ')*', k)
    for f in ('EXP', 'dot', 'pw', 'V', 'D', 'PL0', 'PL1', 'PLQ', 'PLQK', 'PLG', 'DPERP', 'K', 'p', 'q'):
        pass
    k = re.sub(r'(?<![A-Za-z_])(PLQK|PLQ|PLG|PL0|PL1|DPERP|p|q|K)(?=[\(A-Za-z])', r'\1*', k)
    for f in ('EXP', 'DOT', 'PW', 'V', 'D'):
        k = k.replace(f + '*(', f + '(')
    k = re.sub(r'\)\^\((\d)\)', r')**\1', k)
    k = re.sub(r'(\d)(?=[A-Za-z\(])', r'\1*', k)
    k = k.strip('*')
    assert k.count('(') == k.count(')'), k
    return k

def indices(code):
    return sorted(set(re.findall(r'\b(ii|jj|kk|mm|kq|mq)\b', code)))
