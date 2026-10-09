"""Explicit kernels: substitution of positions and LaTeX rendering."""
import sympy as sp


def subst_label(lab, old, new):
    return new if lab == old else lab


def subst_kernel(kexpl, old, new):
    out = []
    for c, facs in kexpl:
        nf = []
        for f in facs:
            if f[0] == 'KK':
                nf.append(('KK',) + tuple(subst_label(x, old, new) for x in f[1:]))
            else:
                nf.append(('raw', f[1], [subst_label(x, old, new) for x in f[2]]))
        out.append((c, nf))
    return out


def tex_factor(f):
    if f[0] == 'KK':
        p, q, r, s = f[1:]
        if (p, q) == (r, s):
            return r'\frac{1}{(%s-%s)^2}' % (p, q)
        if (p, q) == (s, r):
            return r'\Big(-\frac{1}{(%s-%s)^2}\Big)' % (p, q)
        return r'\KK(%s-%s)\cdot\KK(%s-%s)' % (p, q, r, s)
    out = f[1]
    for n, a in enumerate(f[2]):
        out = out.replace('@%d' % n, a)
    return out


def tex_coef(c):
    c = sp.nsimplify(c)
    if c == 1:
        return ''
    if c == sp.Rational(1, 2):
        return r'\tfrac12\,'
    return sp.latex(c) + r'\,'


def tex_kernel(kexpl):
    """render sum_k c_k prod(factors), pulling out factors common to all terms."""
    keys = [[tex_factor(f) for f in facs] for _, facs in kexpl]
    common = []
    if len(keys) > 1:
        for k in keys[0]:
            if all(k in kk for kk in keys[1:]):
                common.append(k)
    parts = []
    for (c, _), kk in zip(kexpl, keys):
        rest = list(kk)
        for k in common:
            rest.remove(k)
        body = r'\,'.join(rest) if rest else '1'
        parts.append(tex_coef(c) + body)
    if len(parts) == 1:
        return parts[0]
    inner = '+'.join(parts)
    return r'\,'.join(common) + r'\,\Big[' + inner + r'\Big]'
