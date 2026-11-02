r"""Your Other_terms.tex (input/), split into equations and additive terms, translated for pointwise evaluation.

Every equation (align block that is not inside \begin{comment}) is split at its top level (inside its outer
\Bigg\{ ... \Bigg\} when it has one) into the terms that start a line with + or -.  The outer factors are kept on every
term, so  sum(terms) = equation.  p+ / q+ integrals are kept as a parameter: \int dp^+ is dropped, p^+ -> PP, q^+ -> QQ,
[k^+/p^+]_+ -> PLUSK, [k^+/q^+]_+ -> PLUSKQ, [1/q^+]_+ -> PLUS1Q (evaluated with the subtraction at p^+ = 0, uvscan.py)."""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'llsplit'))
sys.path.insert(0, os.path.join(HERE, '..', 'bfkl_tex_check'))
import texeval as TE  # noqa: E402

SRC = open(os.path.join(HERE, 'input', 'Other_terms.tex')).read()


def live_text():
    t = re.sub(r'\\begin\{comment\}.*?\\end\{comment\}', '', SRC, flags=re.S)
    return t[t.index(r'\begin{document}'):]


def headings_and_blocks():
    """[(heading path, block text)] in document order"""
    t = live_text()
    out, path = [], {}
    pos = 0
    for m in re.finditer(r'\\(section|subsection|subsubsection)\*?\{([^}]*)\}|\\begin\{align\}(.*?)\\end\{align\}', t, re.S):
        if m.group(1):
            lvl = m.group(1)
            path[lvl] = m.group(2).replace(r'\,', ' ').strip()
            if lvl == 'section':
                path.pop('subsection', None)
                path.pop('subsubsection', None)
            if lvl == 'subsection':
                path.pop('subsubsection', None)
        else:
            name = ' / '.join(path[k] for k in ('section', 'subsection', 'subsubsection') if k in path)
            out.append((name, m.group(3)))
    return out


OPEN = {'{': '}', '(': ')', '[': ']'}


def depth_scan(t):
    """depth before each character; \\{ \\} count as brackets, \\left/\\right ignored (their delimiters count)"""
    d, ds, i = 0, [], 0
    while i < len(t):
        ds.append(d)
        c = t[i]
        if c == '\\' and i + 1 < len(t) and t[i + 1] in '{}':
            d += 1 if t[i + 1] == '{' else -1
            ds.append(d)
            i += 2
            continue
        if c in '{([':
            d += 1
        elif c in '})]':
            d -= 1
        i += 1
    return ds


def top_split(body):
    """split at line starts (after \\\\ and &) whose first sign is + or -, at depth 0 of body"""
    ds = depth_scan(body)
    cuts = [0]
    for m in re.finditer(r'\\\\\s*&?\s*(?=[+-])', body):
        k = m.end()
        if ds[k] == 0:
            cuts.append(k)
    # also a sign right after the opening (first term) is fine
    parts = [body[a:b] for a, b in zip(cuts, cuts[1:] + [len(body)])]
    return [p for p in parts if re.sub(r'\\notag|\\\\|&|\s', '', p)]


def strip_lhs(b):
    i = b.index('=')
    return b[i + 1:]


def outer_group(e):
    """(prefix, inner, suffix) for the first top-level \\Bigg\\{ ... \\Bigg\\} that spans several lines, else None"""
    m = re.search(r'\\(Bigg|bigg)\\\{', e)
    if not m:
        return None
    ds = depth_scan(e)
    if ds[m.start()] != 0:
        return None
    k = m.end()
    start_depth = ds[k]
    j = k
    while j < len(ds) and ds[j] >= start_depth:
        j += 1
    close = j - 1
    assert e.startswith(r'\}', close), e[close - 20:close + 5]
    inner = e[k:close]
    pre = e[:m.start()]
    suf = e[close + 2:]
    # remove the size command before \}
    pre_inner = re.sub(r'\\(Bigg|bigg)\s*$', '', inner)
    return pre, pre_inner, suf


def equation_terms(block):
    e = strip_lhs(block)
    parts = top_split(e)
    if len(parts) > 1:
        return parts
    g = outer_group(e)
    if g is None:
        return parts
    pre, inner, suf = g
    sub = top_split(inner)
    return [pre + r' \Bigg\{ ' + s + r' \Bigg\} ' + suf for s in sub]


PLUS = [(r'\left[\frac{k^{+}}{p^{+}}\right]_{+}', ' PLUSK '), (r'\left[\frac{k^{+}}{q^{+}}\right]_{+}', ' PLUSKQ '),
        (r'\left[\frac{1}{q^{+}}\right]_{+}', ' PLUS1Q '),
        (r'\left[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\right]_{+}', ' PLUSKP ')]


def prep(t):
    """longitudinal integrals -> parameters"""
    for a, b in PLUS:
        t = t.replace(a, b)
    t = re.sub(r'\\int_\{[^{}]*\}\^\{[^{}]*\}\s*\\?,?\s*\\frac\{d[pq]\^\{\+\}\}\{2\\pi\}', ' ', t)
    t = re.sub(r'\\int\^\{[^{}]*\}_\{[^{}]*\}\s*\\?,?\s*\\frac\{d[pq]\^\{\+\}\}\{2\\pi\}', ' ', t)
    t = re.sub(r'\\int_\{[^{}]*(\{[^{}]*\}[^{}]*)*\}\^\{[^{}]*(\{[^{}]*\}[^{}]*)*\}\s*\\?,?\s*\\frac\{d[pq]\^\{\+\}\}\{2\\pi\}', ' ', t)
    t = t.replace('p^{+}', ' PP ').replace('q^{+}', ' QQ ')
    # e^{-ik.(PP/k+)(a-b)} -> e^{-ik.(a-b) PP/k+}
    t = re.sub(r'e\^\{-ik\s*\\cdot\s*(\\frac\{\s*PP\s*\}\{k\^\{\+\}\})\s*(\\,)?\s*(\\left)?\(([^()]*?)(\\right)?\)\s*\}',
               lambda m: 'e^{-ik\\cdot (%s)%s}' % (m.group(4), m.group(1)), t)
    return t


def code_of(t):
    return TE.translate(prep(t))


def all_equations():
    """[(name, n, [term texts])]"""
    out = []
    cnt = {}
    for name, b in headings_and_blocks():
        if r'\mathcal{A}_{0} & =' in b or b.strip().startswith('R^{i}'):
            continue
        cnt[name] = cnt.get(name, 0) + 1
        out.append((name, cnt[name], equation_terms(b)))
    return out


if __name__ == '__main__':
    for name, n, ts in all_equations():
        ok = []
        for t in ts:
            try:
                code_of(t)
                ok.append('ok')
            except Exception as ex:
                ok.append('ERR:' + str(ex)[:60])
        print(name, n, len(ts), ok)
