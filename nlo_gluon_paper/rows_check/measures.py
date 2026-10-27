r"""Integration measures: every transverse point that appears in an equation should be integrated, and every integrated
point should appear.  Amplitude-level equations (with \mf{...}) are skipped.  Output: measures.json"""
import json
import re
import pe

PT = r"(x'|y'|w'|z'|x|y|z|w)"


def norm(t):
    return t.replace(r'^{\prime}', "'").replace(r'^\prime', "'")


def used(t):
    t = norm(t)
    pts = set()
    for a, b in re.findall(r"\(\s*" + PT + r"\s*-\s*" + PT + r"\s*\)", t):
        pts |= {a, b}
    for a, b in re.findall(r"\\left\(\s*" + PT + r"\s*-\s*" + PT + r"\s*\\right\)", t):
        pts |= {a, b}
    for p in re.findall(r"(?:U|\\rho)\^\{[^}]*\}\s*\(\s*" + PT + r"\s*\)", t):
        pts.add(p)
    return pts


def measured(t):
    t = norm(t)
    pts = set()
    for sub in re.findall(r"\\int\s*_\{([^{}]*)\}", t):
        if re.fullmatch(r"[\sxyzw',\\,]*", sub):
            pts |= set(re.findall(r"[xyzw]'?", sub.replace(' ', '')))
    for sub in re.findall(r"\\int\s*_([xyzw])(?![a-z])", t):
        pts.add(sub)
    return pts


def equations(g):
    L = pe.LINES[g]
    i = 0
    while i < len(L):
        if r'\begin{align}' in L[i] and not L[i].lstrip().startswith('%'):
            j = i
            while r'\end{align}' not in L[j]:
                j += 1
            yield i + 1, j + 1, '\n'.join(L[i:j + 1])
            i = j
        i += 1


if __name__ == '__main__':
    out = {}
    for g in ('Group_I', 'Group_II', 'Group_III'):
        for a, b, t in equations(g):
            if r'\mf{' in t or r'R^{i} =' in t or '=' not in t:
                continue
            t2 = re.sub(r'(?<!\\)%[^\n]*', '', t)
            u, m = used(t2), measured(t2)
            if u != m:
                out['%s %d-%d' % (g, a, b)] = dict(not_integrated=sorted(u - m), integrated_not_used=sorted(m - u))
                print(g, a, b, 'missing:', sorted(u - m), 'extra:', sorted(m - u))
    json.dump(out, open('measures.json', 'w'), indent=1)
