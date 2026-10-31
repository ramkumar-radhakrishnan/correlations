r"""The terms (B) of Group III Rows III(1), III(2), V, VI as LaTeX (../bfkl_final/static.tex, 'The new terms').
(B) has y' -> z in the colour and K(y'-z) -> K(w'-z); the y' left in the text there is w'."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
LVK = r'\log{\left(\frac{\vee}{k^{+}}\right)}'


def hidden_B_text():
    S = open(os.path.join(HERE, '..', 'bfkl_final', 'static.tex')).read()
    i = S.index(r'\paragraph{The new terms}')
    blk = S[S.index(r'\begin{align}', i) + len(r'\begin{align}'):S.index(r'\end{align}', i)]
    parts = re.split(r'\\text\{Row ([IV()12]+):\}\\quad&', blk)[1:]
    out = {}
    for k in range(0, len(parts), 2):
        t = parts[k + 1].replace(r'\Lv', LVK).replace(r'\nc', '')
        out['G3 ' + parts[k]] = re.sub(r'[,.]\s*(\\\\)?\s*$', '', t.strip()).replace("y'", "w'")
    return out
