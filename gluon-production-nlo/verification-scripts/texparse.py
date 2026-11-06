"""Parser for the canonical-form LaTeX terms used in the Evolution_* files.

A 'block' is the body of an align environment after its header (prefactor etc.).
Terms are split at top-level +/- signs. Each term may contain:
  numeric coefficient (integer or \frac{1}{2}, \frac{3}{2}...)
  N_c factors: N_{c}
  \mathcal{K}(A-B)\cdot\mathcal{K}(C-D)  (possibly chained with \cdot typos)
  \frac{1}{(A-B)^{2}}
  f^{abc}
  U^{ij}(pos)
  \rho^{i}(pos),  \big\{\rho^{i}(p),\rho^{j}(q)\big\}  (also \Big\{..\Big\})
Positions: x,y,z,u,w,x',y',w' -> x,y,z,u,w,xp,yp,wp
Indices: letters with optional prime -> letter or letter+'p'
"""
import re
from engine import Term

def _normpos(p):
    p = p.strip().replace("\\,", "").replace(" ", "")
    return {"x'": "xp", "y'": "yp", "w'": "wp", "u'": "up"}.get(p, p)

def _normidx(s):
    s = s.strip()
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if not c.isalpha():
            i += 1
            continue
        if i + 1 < len(s) and s[i + 1] == "'":
            out.append(c + "p")
            i += 2
        else:
            out.append(c)
            i += 1
    return out

def clean(tex):
    tex = tex.replace("\n", " ")
    tex = re.sub(r"\\notag|\\nonumber|\\\\|&", " ", tex)
    tex = re.sub(r"\\color\{\w+\}", "", tex)
    tex = re.sub(r"\\(alpha|beta)'?", "", tex)
    tex = tex.replace("\\quad", " ").replace("\\;", " ")
    tex = tex.replace("\\Bigg[", " ").replace("\\Bigg]", " ")
    tex = tex.replace("\\Big[", " ").replace("\\Big]", " ")
    # keep \big\{ markers for anticommutators
    return tex

def split_terms(body):
    """split at top-level + / - (depth 0 wrt {} and ())"""
    body = clean(body)
    terms, cur, depth = [], "", 0
    i = 0
    while i < len(body):
        c = body[i]
        if c in "{(":
            depth += 1
        elif c in ")}":
            depth -= 1
        if c in "+-" and depth == 0:
            # not inside ^{...} etc. since depth==0; also skip e^{-ik...} (none in body)
            if cur.strip():
                terms.append(cur)
            cur = c
        else:
            cur += c
        i += 1
    if cur.strip():
        terms.append(cur)
    # drop empties / stray
    return [t for t in terms if re.search(r"\\rho", t)]

FRAC = re.compile(r"\\frac\{(\d+)\}\{(\d+)\}")
KPAT = re.compile(r"\\mathcal\{K\}\(([^)]*)\)")
INV2 = re.compile(r"\\frac\{1\}\{\(([^)]*)\)\^\{2\}\}")
FPAT = re.compile(r"f\^\{([^}]*)\}")
UPAT = re.compile(r"U\^\{([^}]*)\}\(([^)]*)\)")
ACOM = re.compile(r"\\[bB]ig\\\{\s*\\rho\^\{([^}]*)\}\(([^)]*)\)\s*,\s*\\rho\^\{([^}]*)\}\(([^)]*)\)\s*\\[bB]ig\\\}")
RPAT = re.compile(r"\\rho\^\{([^}]*)\}\(([^)]*)\)")
NCPAT = re.compile(r"N_\{c\}")

def _splitdiff(arg):
    # "x'-w'" -> ("xp","wp")
    parts = arg.split("-")
    assert len(parts) == 2, arg
    return _normpos(parts[0]), _normpos(parts[1])

def parse_term(ttex, extra_coef=1.0, tag=""):
    """Returns list of Term (anticommutators expanded). rho order = order of appearance."""
    s = ttex
    sign = 1.0
    ms = re.match(r"\s*([+-])", s)
    if ms:
        if ms.group(1) == "-":
            sign = -1.0
        s = s[ms.end():]
    coef = sign * extra_coef
    mfrac = FRAC.search(s)
    if mfrac:
        coef *= int(mfrac.group(1)) / int(mfrac.group(2))
        s = s[:mfrac.start()] + " " + s[mfrac.end():]
    else:
        mint = re.match(r"\s*(\d+)\s*\\,", s)
        if mint:
            coef *= int(mint.group(1))
            s = s[mint.end():]
    ncp = len(NCPAT.findall(s))
    s = NCPAT.sub(" ", s)

    kern = []
    kargs = KPAT.findall(s)
    assert len(kargs) % 2 == 0, f"odd K count in: {ttex}"
    for ia in range(0, len(kargs), 2):
        a, b = _splitdiff(kargs[ia])
        c, d = _splitdiff(kargs[ia + 1])
        kern.append(('K', a, b, c, d))
    for m in INV2.findall(s):
        a, b = _splitdiff(m)
        kern.append(('inv2', a, b))

    color = []
    for m in FPAT.findall(s):
        idx = _normidx(m)
        assert len(idx) == 3, (m, ttex)
        color.append(('f', *idx))
    for m in UPAT.findall(s):
        idx = _normidx(m.group(1) if hasattr(m, 'group') else m[0])
        pos = _normpos(m[1] if isinstance(m, tuple) else m.group(2))
        assert len(idx) == 2, (m, ttex)
        color.append(('U', idx[0], idx[1], pos))

    # rho structure in order of appearance: build event list over the string
    events = []
    for m in ACOM.finditer(s):
        i1 = _normidx(m.group(1)); p1 = _normpos(m.group(2))
        i2 = _normidx(m.group(3)); p2 = _normpos(m.group(4))
        events.append((m.start(), 'A', (i1[0], p1), (i2[0], p2)))
    s_noac = ACOM.sub(" ", s)
    for m in RPAT.finditer(s_noac):
        idx = _normidx(m.group(1)); pos = _normpos(m.group(2))
        events.append((m.start(), 'R', (idx[0], pos)))
    # positions in s vs s_noac differ slightly; order is still right because ACOM
    # removal replaces with single space -- but offsets shift. Redo properly:
    events = []
    masked = s
    spans = []
    for m in ACOM.finditer(s):
        spans.append((m.start(), m.end(),
                      ('A', (_normidx(m.group(1))[0], _normpos(m.group(2))),
                            (_normidx(m.group(3))[0], _normpos(m.group(4))))))
    def inspan(i):
        return any(a <= i < b for a, b, _ in spans)
    for m in RPAT.finditer(s):
        if not inspan(m.start()):
            spans.append((m.start(), m.end(), ('R', (_normidx(m.group(1))[0], _normpos(m.group(2))))))
    spans.sort()
    seq = [ev for _, _, ev in spans]

    # handle repeated indices between rhos (e.g. rho^a(x) rho^a(x')) -> rename + delta
    out_terms = []
    def build(rholist, cf):
        rl, colr = [], list(color)
        seen = {}
        cnt = [0]
        for (ix, p) in rholist:
            if ix in seen or any(ix in c[1:-1] if c[0] == 'U' else ix in c[1:] for c in []):
                pass
            if ix in seen:
                nix = ix + f"r{cnt[0]}"; cnt[0] += 1
                colr.append(('delta', ix, nix))
                rl.append((nix, p))
            else:
                seen[ix] = True
                rl.append((ix, p))
        return Term(cf, kern, colr, rl, Nc_pow=ncp, tag=tag)

    # expand anticommutators
    combos = [([], coef)]
    for ev in seq:
        if ev[0] == 'R':
            combos = [(lst + [ev[1]], cf) for lst, cf in combos]
        else:
            new = []
            for lst, cf in combos:
                new.append((lst + [ev[1], ev[2]], cf))
                new.append((lst + [ev[2], ev[1]], cf))
            combos = new
    # NOTE: {A,B} = AB + BA exactly (no 1/2) -- any explicit 1/2 is already in the tex
    for lst, cf in combos:
        out_terms.append(build(lst, cf))
    return out_terms

def parse_block(body, extra_coef=1.0, tag=""):
    out = []
    for i, t in enumerate(split_terms(body)):
        out.extend(parse_term(t, extra_coef=extra_coef, tag=f"{tag}#{i}"))
    return out
