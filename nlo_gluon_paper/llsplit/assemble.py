r"""Rebuild align environments from selected (and modified) terms of a block.

assemble(B, items): items = list of (text, style), style None | 'red' | 'blue'.
Returns the LaTeX of one align environment: your left-hand side, the block prefix, the terms (in their order),
the closing \Bigg\} and the common colour factor if the block has one."""
import re

import gen as GN


def mark(text, cmd):
    """put \\rc or \\nc at the start of every align cell of a term"""
    parts = re.split(r'(?<!\\)&', text)
    out = [(cmd + ' ' + parts[0]) if parts[0].strip() else parts[0]]
    for p in parts[1:]:
        out.append(' ' + cmd + ' ' + p)
    return '&'.join(out)


def clean(t):
    return GN.strip_user_colours(t)


def assemble(B, items, lhs=None):
    if not items:
        return ''
    lhs = B.lhs if lhs is None else lhs
    texts = []
    for k, (t, style) in enumerate(items):
        t = clean(t).rstrip()
        t = re.sub(r'(\\notag|\\nonumber)?\s*\\\\\s*$', '', t)
        s = t.lstrip()
        first_inline = (k == 0 and not B.braced)
        if first_inline:
            if s.startswith('&'):
                s = s[1:]
        else:
            if not s.startswith('&'):
                s = '& ' + s
        if style == 'red':
            s = mark(s, r'\rc')
        elif style == 'blue':
            s = mark(s, r'\nc')
        texts.append(s)
    body = (r' \notag \\' + '\n').join(texts)
    prefix = clean(B.prefix)
    suffix = clean(B.suffix)
    if B.braced:
        out = lhs + ' ' + prefix + r'\Bigg\{\notag \\' + '\n' + body + r'\Bigg\}' + suffix
    else:
        out = lhs + ' ' + prefix + body
    return '\\begin{align}\n' + out.rstrip() + '\n\\end{align}\n'
