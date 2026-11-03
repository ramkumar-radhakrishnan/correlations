"""LaTeX of singlet3.pkl -> singlet3_frag.tex (one macro-free block per piece)."""
import os
import pickle
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'finite_singlet'))
sys.path.insert(0, os.path.join(HERE, '..', 'singlet'))
sys.modules.pop('trace', None)
import srender as SR  # noqa: E402

st = pickle.load(open(os.path.join(HERE, 'singlet3.pkl'), 'rb'))


def out():
    lines = []
    for name, d in st['two'].items():
        lines.append('%% TWO %s\n' % name + SR.bracket(d['comb']))
    for name, g in st['three'].items():
        for pair, d in g.items():
            lines.append('%% THREE %s %s<-%s\n' % (name, pair[0], pair[1]) + SR.bracket(d['comb']))
    return '\n\n'.join(lines)


if __name__ == '__main__':
    s = out()
    open(os.path.join(HERE, 'singlet3_frag.tex'), 'w').write(s)
    print(s)
