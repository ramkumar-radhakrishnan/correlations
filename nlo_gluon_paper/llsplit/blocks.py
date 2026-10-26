r"""Your final row expressions (Group_I/II/III.tex), split into terms.

A block is the right-hand side of one align environment (after '& =').  If the terms are enclosed in \Bigg\{ ... \Bigg\}
the text before is the prefix (common factor) and the text after the closing \Bigg\} the suffix (common colour factor).
Inside, the text is cut at the LaTeX line breaks; a term starts at a line that begins with an optional sign and a
prefactor (\frac{1}{(2\pi)^{3}}, \frac{ig^{4}..}, \frac{g^{4}..}, \frac{g^{2}}{4\pi^{2}}).
Value of the block = sum over terms of  prefix x term x suffix."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n')
         for g in ('Group_I', 'Group_II', 'Group_III', '4rho')}

# final blocks: (group, first line, last line) of the align environment (1-based, inclusive)
FINAL = {
    'G1 II(1)': ('Group_I', 49, 51), 'G1 II(2a)': ('Group_I', 68, 73), 'G1 II(2b)': ('Group_I', 74, 79),
    'G1 III N1': ('Group_I', 95, 97), 'G1 III N2': ('Group_I', 98, 100), 'G1 III N3a': ('Group_I', 101, 120),
    'G1 III N3b': ('Group_I', 121, 129),
    'G1 IV': ('Group_I', 142, 145),
    'G1 V(1)': ('Group_I', 244, 246), 'G1 V(2a)': ('Group_I', 248, 250), 'G1 V(2b UV)': ('Group_I', 251, 255),
    'G1 V(2b fin)': ('Group_I', 256, 263),
    'G2 I(1)': ('Group_II', 49, 51),
    'G2 I(2a)': ('Group_II', 66, 68), 'G2 I(2b)': ('Group_II', 69, 71), 'G2 I(2c)': ('Group_II', 72, 74),
    'G2 I(2d)': ('Group_II', 75, 77), 'G2 I(2e)': ('Group_II', 78, 80), 'G2 I(2f)': ('Group_II', 81, 83),
    'G2 IV': ('Group_II', 132, 134),
    'G3 I': ('Group_III', 47, 49), 'G3 II': ('Group_III', 58, 60), 'G3 III(1)': ('Group_III', 75, 77),
    'G3 III(2)': ('Group_III', 89, 93), 'G3 V': ('Group_III', 112, 114), 'G3 VI': ('Group_III', 127, 129),
}

TERM_START = re.compile(r"^\s*(?:\{\\color\{red\}\s*)?(?:[+-])?\s*(?:\\,)?\s*"
                        r"(?=\\frac\{1\}\{\(2\\pi\)\^\{3\}\}|\\frac\{i\s*g\^\{4\}|\\frac\{g\^\{4\}|\\frac\{g\^\{2\}\}\{4\\pi\^\{2\}\}"
                        r"|\\frac\{i\\,g\^\{4\}|\\frac\{ig\^\{4\})")


def align_text(name):
    g, a, b = FINAL[name]
    t = '\n'.join(FILES[g][a - 1:b])
    i = t.index(r'\begin{align}') + len(r'\begin{align}')
    j = t.rindex(r'\end{align}')
    return t[i:j]


def split_lines(t):
    parts = re.split(r'(\\\\)', t)
    out, cur = [], ''
    for p in parts:
        cur += p
        if p == '\\\\':
            out.append(cur)
            cur = ''
    out.append(cur)
    return out


def cell_text(line):
    s = re.sub(r'^\s*&', '', line)
    return re.sub(r'(\\notag|\\nonumber)?\s*\\\\\s*$', '', s)


def bigg_pair(t):
    """position of the first \\Bigg\\{ and of its matching \\Bigg\\}"""
    i = t.find(r'\Bigg\{')
    if i < 0:
        return None
    d, k = 0, i
    while k < len(t):
        if t.startswith(r'\Bigg\{', k):
            d += 1
            k += 7
            continue
        if t.startswith(r'\Bigg\}', k):
            d -= 1
            if d == 0:
                return i, k
            k += 7
            continue
        k += 1
    return i, None


def pure_prefix(t):
    """only constants, measures and the phase (no kernel, colour or log)"""
    while 'e^{' in t:
        i = t.index('e^{')
        d, k = 0, i + 2
        while True:
            d += {'{': 1, '}': -1}.get(t[k], 0)
            k += 1
            if d == 0:
                break
        t = t[:i] + t[k:]
    return not re.search(r'\\log|U\^|\\rho|\((x|y|z|w)|\\left\[|\\Bigg\[|\\bigg\[', t)


class Block:
    def __init__(self, name, text=None):
        self.name = name
        t = text if text is not None else align_text(name)
        m = re.search(r'&\s*=', t)
        self.lhs = t[:m.end()]
        rhs = t[m.end():]
        bp = bigg_pair(rhs)
        first_term = self._first_term_pos(rhs)
        if bp and pure_prefix(rhs[:bp[0]]):
            self.prefix = rhs[:bp[0]]
            if bp[1] is None:
                inner, self.suffix = rhs[bp[0] + 7:], ''
            else:
                inner = rhs[bp[0] + 7:bp[1]]
                self.suffix = rhs[bp[1] + 7:]
            self.braced = True
        else:
            self.prefix, inner, self.suffix = rhs[:first_term], rhs[first_term:], ''
            self.braced = False
        self.inner = inner
        self.terms = self._split_terms(inner)

    @staticmethod
    def _first_term_pos(rhs):
        pos = 0
        for line in split_lines(rhs):
            s = cell_text(line)
            off = len(line) - len(line.lstrip())
            m = TERM_START.match(s)
            if m:
                return pos + line.index(s.lstrip()[:20]) if s.strip() else pos
            pos += len(line)
        return len(rhs)

    @staticmethod
    def _split_terms(inner):
        lines = split_lines(inner)
        starts = [k for k, ln in enumerate(lines) if TERM_START.match(cell_text(ln))]
        first = next((k for k, ln in enumerate(lines) if cell_text(ln).strip()), 0)
        if not starts or starts[0] != first:
            starts = [first] + starts
        out = []
        for n, k in enumerate(starts):
            e = starts[n + 1] if n + 1 < len(starts) else len(lines)
            out.append(''.join(lines[k:e]))
        out = [t for t in out if cell_text(t).strip()]
        # a leading line with only \notag \\ belongs to nothing
        return out

    def term_full(self, n, text=None):
        """prefix x term x suffix, as one LaTeX expression"""
        t = self.terms[n] if text is None else text
        return self.prefix + ' (' + t + ') ' + self.suffix


def logs_in(t):
    """the \\log factors in a string (brace-matched), as written"""
    out = []
    for m in re.finditer(r'\\log', t):
        i = m.start()
        j = i + 4
        while j < len(t) and t[j] == ' ':
            j += 1
        if t.startswith('^{2}', j):
            j += 4
        if t[j] == '{':
            d, k = 1, j + 1
            while d:
                d += {'{': 1, '}': -1}.get(t[k], 0)
                k += 1
            out.append(t[i:k])
        else:
            for o, c in ((r'\left(', r'\right)'), (r'\left[', r'\right]'), (r'\Big(', r'\Big)'), ('(', ')')):
                if t.startswith(o, j):
                    d, k = 1, j + len(o)
                    while d:
                        if t.startswith(o, k):
                            d += 1
                            k += len(o)
                            continue
                        if t.startswith(c, k):
                            d -= 1
                            k += len(c)
                            continue
                        k += 1
                    out.append(t[i:k])
                    break
            else:
                out.append(t[i:j + 40])
    return out


def g2_i1_parts():
    """Group II Row I part 1 has two blocks, each with its own \\Bigg\\{..\\Bigg\\} and colour"""
    t = align_text('G2 I(1)')
    m = re.search(r"&\s*\+\s*(?=\\frac\{1\}\{\(2\\pi\)\^\{3\}\}\\,\\int_\{x,y,z,w,x',w\^\{\\prime\}\})", t)
    a = t[:m.start()]
    b = r'\frac{d^{3}N}{d^{3}k} & = ' + t[m.end():]
    return Block('G2 I(1) A', a), Block('G2 I(1) B', b)


def blocks():
    out = {}
    for name in FINAL:
        if name == 'G2 I(1)':
            A, B = g2_i1_parts()
            out[A.name], out[B.name] = A, B
        else:
            out[name] = Block(name)
    return out


if __name__ == '__main__':
    for name, B in blocks().items():
        print('=' * 10, name, '| terms', len(B.terms), '| braced', B.braced, '| prefix:', re.sub(r'\s+', ' ', B.prefix)[:60],
              '| suffix:', re.sub(r'\s+', ' ', B.suffix)[:50])
        for n, tt in enumerate(B.terms):
            lg = [re.sub(r'\s+', '', x)[:60] for x in logs_in(tt)]
            print('  %2d %-60s %s' % (n, re.sub(r'\s+', ' ', cell_text(tt))[:60], lg))
