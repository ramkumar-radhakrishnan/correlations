r"""Integration measures of every displayed equation: points in the measures vs points that appear in the integrand."""
import re
from rows_user import rows

P = r"(x'|y'|w'|z'|x\^\{\\prime\}|w\^\{\\prime\}|w\^\\prime|y\^\{\\prime\}|x|y|z|w)"


def norm(p):
    return p.replace(r'^{\prime}', "'").replace(r'^\prime', "'")


def scan(f, stop=None):
    out = []
    for g, r, p, ln, t in rows(f, stop):
        t = t.replace(r'^{\prime}', "'").replace(r'^\prime', "'")
        meas = re.findall(r'\\int_\{([^}]*(?:\{\\prime\}[^}]*)*)\}', t)
        ms = []
        for m in meas:
            m = m.replace(r'^{\prime}', "'").replace(r'^\prime', "'")
            ms.append(m)
        body = re.sub(r'\\int_\{[^}]*(?:\{\\prime\}[^}]*)*\}', ' ', t)
        body = re.sub(r'\\frac\{dp\^\{\+\}\}', ' ', body)
        used = set(norm(q) for q in re.findall(r'[\(,]\s*' + P + r"\s*[-)]", body)) | \
            set(norm(q) for q in re.findall(r'-\s*' + P + r'\s*\)', body))
        listed = set()
        bad = []
        for m in ms:
            items = [s.strip() for s in m.split(',')]
            if '' in items:
                bad.append('empty item in \\int_{%s}' % m)
            for s in items:
                if s and s not in ("x", "y", "z", "w", "x'", "y'", "w'", "z'"):
                    bad.append('item "%s" in \\int_{%s}' % (s, m))
            listed |= set(s for s in items if s)
        miss, extra = used - listed, listed - used
        if bad or miss or extra:
            out.append((ln, g, r, p, sorted(miss), sorted(extra), bad))
    return out


if __name__ == '__main__':
    for f, s in (('Evolution_two_rho.tex', None), ('Evolution_three_rho.tex', None)):
        for e in scan(f, s):
            print(f[10:13], e)
