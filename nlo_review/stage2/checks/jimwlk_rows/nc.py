# Non-commutative rho algebra helper: parse LaTeX monomials (kernels * f's * U's * ordered rho's),
# evaluate numerically at assignments of named points to physical points, extract contact (two-rho) terms.
import re, itertools, numpy as np
from scipy.linalg import expm
UP = "/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
LET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

def su(N=3):
    gens = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex); m[i, j] = m[j, i] = 0.5; gens.append(m)
            m = np.zeros((N, N), complex); m[i, j] = -0.5j; m[j, i] = 0.5j; gens.append(m)
    for k in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(k): m[i, i] = 1
        m[k, k] = -k; gens.append(m / np.sqrt(2 * k * (k + 1)))
    t = np.array(gens); n = len(t)
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            c = t[a] @ t[b] - t[b] @ t[a]
            for d in range(n): f[a, b, d] = np.real(-2j * np.trace(c @ t[d]))
    return t, f

class World:
    """physical points: 'z','w',"w'" and P0..P3 ; each has a position and an adjoint Wilson line"""
    def __init__(self, seed=0, N=3, keys=("z", "w", "w'", "P0", "P1", "P2", "P3")):
        rng = np.random.default_rng(seed)
        self.N = N; self.t, self.f = su(N); self.n = len(self.t)
        self.X = {}; self.U = {}
        for k in keys:
            self.X[k] = rng.normal(size=2)
            H = sum(rng.normal() * self.t[a] for a in range(self.n)); V = expm(1j * H)
            self.U[k] = np.real(np.array([[2 * np.trace(self.t[a] @ V @ self.t[b] @ V.conj().T)
                                           for b in range(self.n)] for a in range(self.n)]))
    def swapped_ww(self):
        """same world with w <-> w' interchanged (for complex conjugation)"""
        o = World.__new__(World); o.N = self.N; o.t = self.t; o.f = self.f; o.n = self.n
        o.X = dict(self.X); o.U = dict(self.U)
        o.X["w"], o.X["w'"] = self.X["w'"], self.X["w"]; o.U["w"], o.U["w'"] = self.U["w'"], self.U["w"]
        return o

def _vec(expr):
    out = []
    for sgn, name in re.findall(r"([+-]?)\s*([a-z]'?)", expr):
        out.append((-1 if sgn == '-' else 1, name))
    return out

class Mono:
    __slots__ = ("coef", "dots", "invs", "kfun", "fs", "Us", "rhos", "ncp", "tag", "src")
    def __init__(self):
        self.coef = 1.0; self.dots = []; self.invs = []; self.kfun = None; self.fs = []; self.Us = []
        self.rhos = []; self.ncp = 0; self.tag = None; self.src = ""
    def copy(self):
        m = Mono()
        for s in self.__slots__: setattr(m, s, getattr(self, s) if not isinstance(getattr(self, s), list) else list(getattr(self, s)))
        return m
    def pts(self):
        return tuple(p for _, p in self.rhos)
    def kern(self, world, assign):
        pos = lambda nm: world.X[assign.get(nm, nm)]
        def v(e):
            r = np.zeros(2)
            for s, nm in e: r = r + s * pos(nm)
            return r
        val = 1.0
        for a, b in self.dots:
            va, vb = v(a), v(b); val *= (va @ vb) / ((va @ va) * (vb @ vb))
        for a in self.invs:
            va = v(a); val *= 1.0 / (va @ va)
        if self.kfun is not None:
            val *= self.kfun(lambda nm: pos(nm))
        return val
    def colour(self, world, assign):
        ops = []; subs = []
        for tr in self.fs: ops.append(world.f); subs.append(list(tr))
        for i, j, p in self.Us: ops.append(world.U[assign.get(p, p)]); subs.append([i, j])
        out = []
        cnt = {}
        for s in subs:
            for i in s: cnt[i] = cnt.get(i, 0) + 1
        seen = {}
        for k, (i, p) in enumerate(self.rhos):
            if i in seen or (cnt.get(i, 0) == 0 and sum(1 for j, _ in self.rhos if j == i) == 2 and i in [q for q, _ in self.rhos[:k]]):
                new = i + "#" + str(k); ops.append(np.eye(world.n)); subs.append([i, new]); out.append(new)
            else:
                out.append(i)
            seen[i] = 1
        # sanity: every index appears exactly twice overall (counting output once)
        allc = {}
        for s in subs:
            for i in s: allc[i] = allc.get(i, 0) + 1
        for i in out: allc[i] = allc.get(i, 0) + 1
        bad = [i for i, c in allc.items() if c != 2]
        assert not bad, ("index count", bad, self.src)
        L = {}
        def l(i):
            if i not in L: L[i] = LET[len(L)]
            return L[i]
        spec = ",".join("".join(l(i) for i in s) for s in subs) + "->" + "".join(l(i) for i in out)
        return np.einsum(spec, *ops, optimize=True)
    def eval(self, world, assign):
        return self.coef * (world.N ** self.ncp) * self.kern(world, assign) * self.colour(world, assign)

def _strip_wrappers(s):
    # remove {\color{name} ... } wrappers keeping content
    while True:
        m = re.search(r"\{\\color\{\w+\}", s)
        if not m: break
        i = m.end(); depth = 1; j = i
        while depth:
            c = s[j]
            if c == '{' and s[j - 1] != '\\': depth += 1
            elif c == '}' and s[j - 1] != '\\': depth -= 1
            j += 1
        s = s[:m.start()] + " " + s[i:j - 1] + " " + s[j:]
    return s

def clean(s):
    s = s.replace(r"^{\prime}", "'")
    s = _strip_wrappers(s)
    for x in (r"\notag", r"\nonumber", r"\\", "&", r"\quad", r"\,", r"\;", r"\!", r"\left", r"\right", r"\times",
              r"\Bigg", r"\bigg", r"\Big", r"\big"):
        s = s.replace(x, " ")
    return s

IDX = r"[a-z]'*"
def parse_mono(s, extra_coef=1.0):
    """returns list of Mono (anticommutators expanded). numeric factors anywhere are multiplied."""
    orig = s
    s = clean(s)
    m0 = Mono(); m0.src = orig.strip()[:200]
    mt = re.search(r"\\(alpha|beta)\s*('?)", s)
    if mt:
        m0.tag = mt.group(1) + mt.group(2); s = s[:mt.start()] + s[mt.end():]
    # kernels
    Ks = re.findall(r"\\mathcal\s*\{?K\}?\s*\(([^()]*)\)", s)
    assert len(Ks) % 2 == 0, ("odd K", orig)
    for a, b in zip(Ks[0::2], Ks[1::2]): m0.dots.append((_vec(a), _vec(b)))
    s = re.sub(r"\\mathcal\s*\{?K\}?\s*\(([^()]*)\)", " ", s).replace(r"\cdot", " ")
    for a in re.findall(r"\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}", s): m0.invs.append(_vec(a))
    s = re.sub(r"\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}", " ", s)
    for g in re.findall(r"f\^\{([^{}]*)\}", s):
        tr = re.findall(IDX, g); assert len(tr) == 3, (g, orig); m0.fs.append(tuple(tr))
    s = re.sub(r"f\^\{([^{}]*)\}", " ", s)
    for g, p in re.findall(r"U\s*\^\{([^{}]*)\}\s*\(([a-z]'?)\)", s):
        ij = re.findall(IDX, g.replace(r"\dagger", "")); assert len(ij) == 2, (g, orig); m0.Us.append((ij[0], ij[1], p))
    s = re.sub(r"U\s*\^\{([^{}]*)\}\s*\(([a-z]'?)\)", " ", s)
    # rho items in order
    R1 = r"\\rho\^\{(" + IDX + r")\}\s*\(([a-z]'?)\)"
    items = []
    pat = re.compile(r"\\\{\s*" + R1 + r"\s*,\s*" + R1 + r"\s*\\\}|" + R1)
    for m in pat.finditer(s):
        g = m.groups()
        if g[0] is not None: items.append(("anti", (g[0], g[1]), (g[2], g[3])))
        else: items.append(("one", (g[4], g[5])))
    s = pat.sub(" ", s)
    while "N_{c}" in s:
        m0.ncp += 1; s = s.replace("N_{c}", " ", 1)
    c = extra_coef
    for a, b in re.findall(r"\\frac\{(\d+)\}\{(\d+)\}", s): c *= int(a) / int(b)
    s = re.sub(r"\\frac\{(\d+)\}\{(\d+)\}", " ", s)
    for a in re.findall(r"\d+", s): c *= int(a)
    s = re.sub(r"\d+", " ", s)
    c *= (-1) ** s.count("-")
    left = s.replace("-", "").replace("+", "").replace(r"\-", "").replace("\\", "").strip()
    assert left == "", ("leftover", repr(left), orig)
    m0.coef = c
    outs = [m0]
    # expand rho orderings
    seqs = [[]]
    for it in items:
        if it[0] == "one": seqs = [q + [it[1]] for q in seqs]
        else: seqs = [q + [it[1], it[2]] for q in seqs] + [q + [it[2], it[1]] for q in seqs]
    res = []
    for q in seqs:
        m = m0.copy(); m.rhos = q; res.append(m)
    return res

# ---------- bracket expansion: [A - B] x [C + D(...)] ----------
def expand(s):
    r"""expand products of bracketed sums into a list of (sign, monomial text).
       groups: [ ... ]  and  <( ... )>  (use <( )> for inner \bigg( \bigg) groups).  \{ \} are anticommutators (kept)."""
    s = s.replace(r"\{", "\x01").replace(r"\}", "\x02")
    pos = 0
    def parse_expr(stop):
        nonlocal pos
        terms = []; cur = [(1, "")]; depth = 0; sign = 1
        def flush():
            nonlocal cur, sign
            for sg, t in cur: terms.append((sg * sign, t))
            cur = [(1, "")]; sign = 1
        while pos < len(s):
            ch = s[pos]
            if depth == 0 and ((stop == "]" and ch == "]") or (stop == ")>" and s.startswith(")>", pos))):
                break
            if ch in "({": depth += 1
            if ch in ")}": depth -= 1
            if depth == 0 and ch == "[":
                pos += 1; sub = parse_expr("]"); assert s[pos] == "]"; pos += 1
                cur = [(a * b, t + " " + u) for a, t in cur for b, u in sub]; continue
            if depth == 0 and s.startswith("<(", pos):
                pos += 2; sub = parse_expr(")>"); assert s.startswith(")>", pos); pos += 2
                cur = [(a * b, t + " " + u) for a, t in cur for b, u in sub]; continue
            if depth == 0 and ch in "+-":
                if all(t.strip() == "" for _, t in cur) and len(cur) == 1:
                    if ch == "-": sign = -sign
                else:
                    flush()
                    if ch == "-": sign = -1
                pos += 1; continue
            cur = [(a, t + ch) for a, t in cur]; pos += 1
        flush()
        return [(a, t) for a, t in terms if t.strip() != ""]
    out = parse_expr(None)
    return [(a, t.replace("\x01", r"\{").replace("\x02", r"\}")) for a, t in out]

def monos_from(s, coef=1.0, kfun=None, dots=None):
    """bracketed expression -> list of Mono ; optional common kernel"""
    res = []
    for sg, t in expand(s):
        for m in parse_mono(t, extra_coef=sg * coef):
            if kfun is not None: m.kfun = kfun
            if dots: m.dots = m.dots + [(_vec(a), _vec(b)) for a, b in dots]
            res.append(m)
    return res

def cc(monos):
    """Hermitian conjugate of an operator list, relabelled so that the phase is again e^{-ik(w'-w)}:
       coefficient functions evaluated with w<->w' swapped, rho order reversed, overall i -> -i handled by caller
       (returned Monos carry flag via attribute name 'CC' in src so that evaluation uses the swapped world)."""
    out = []
    for m in monos:
        c = m.copy(); c.rhos = list(reversed(m.rhos)); c.src = "CC:" + m.src; out.append(c)
    return out

def ev(m, world, assign):
    if m.src.startswith("CC:"):
        return m.eval(world.swapped_ww(), assign)
    return m.eval(world, assign)

# ---------- canonical forms ----------
def canon_ordered(monos, world, phys):
    """dict: ordered tuple of physical points -> coefficient tensor (sum over bijections names->phys)"""
    D = {}
    for m in monos:
        names = sorted(set(m.pts()))
        assert len(names) == len(m.rhos) == len(phys), (names, m.src)
        for perm in itertools.permutations(phys):
            a = dict(zip(names, perm))
            T = ev(m, world, a); key = tuple(a[p] for p in m.pts())
            D[key] = D.get(key, 0) + T
    return D

def canon_sym(monos, world, phys):
    """commuting rho's: single tensor with axes ordered as phys"""
    tot = 0
    for key, T in canon_ordered(monos, world, phys).items():
        tot = tot + np.transpose(T, [key.index(p) for p in phys])
    return tot

def two_from_three(monos, world, phys=("P0", "P1")):
    """contact two-rho part of ordered three-rho monomials under full symmetrisation (no overall i):
       rho1 rho2 rho3 -> 1/2 sum_{i<j} [rho_i,rho_j] rho_k ,  [rho^a(p),rho^b(q)] = i f^{abc} rho^c(p) delta(p-q).
       returns M[a,b] = coefficient of rho^a(phys0) rho^b(phys1)"""
    M = np.zeros((world.n, world.n), complex)
    for m in monos:
        pts = m.pts(); assert len(pts) == 3 and len(set(pts)) == 3, m.src
        for i, j in ((0, 1), (0, 2), (1, 2)):
            k = 3 - i - j
            for pm, pn in (phys, phys[::-1]):
                a = {pts[i]: pm, pts[j]: pm, pts[k]: pn}
                T = ev(m, world, a)
                Tm = np.moveaxis(T, (i, j, k), (0, 1, 2))
                X = 0.5 * 1j * np.einsum("abc,abk->ck", world.f, Tm)
                M += X if pm == phys[0] else X.T
    return M

def two_list(monos, world, phys=("P0", "P1")):
    M = np.zeros((world.n, world.n), complex)
    for m in monos:
        pts = m.pts(); assert len(pts) == 2 and len(set(pts)) == 2, m.src
        for pm, pn in (phys, phys[::-1]):
            a = {pts[0]: pm, pts[1]: pn}
            T = ev(m, world, a)
            M += T if pm == phys[0] else T.T
    return M

def two_from_four(monos, world, phys=("P0", "P1")):
    """contact two-rho part of ordered four-rho monomials under full symmetrisation.
       X1X2X3X4 -> sum_{(ij)(kl)} S(1/2[Xi,Xj] 1/2[Xk,Xl]) + sum_{(ijk) l} S(Z_ijk X_l),
       Z_ijk = 1/4[[Xi,Xj],Xk] + 1/12 [Xi,[Xj,Xk]] + 1/12 [Xj,[Xi,Xk]]  (i<j<k)"""
    f = world.f
    M = np.zeros((world.n, world.n), complex)
    for m in monos:
        pts = m.pts(); assert len(pts) == 4 and len(set(pts)) == 4, m.src
        for (i, j), (k, l) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
            for pm, pn in (phys, phys[::-1]):
                a = {pts[i]: pm, pts[j]: pm, pts[k]: pn, pts[l]: pn}
                T = np.moveaxis(ev(m, world, a), (i, j, k, l), (0, 1, 2, 3))
                X = 0.25 * (1j) ** 2 * np.einsum("abc,kld,abkl->cd", f, f, T)
                M += X if pm == phys[0] else X.T
        for l in range(4):
            i, j, k = [q for q in range(4) if q != l]
            for pm, pn in (phys, phys[::-1]):
                a = {pts[i]: pm, pts[j]: pm, pts[k]: pm, pts[l]: pn}
                T = np.moveaxis(ev(m, world, a), (i, j, k, l), (0, 1, 2, 3))
                Z = 0.25 * np.einsum("abc,ckd,abkl->dl", f, f, T) \
                    + (1 / 12) * np.einsum("bkc,acd,abkl->dl", f, f, T) \
                    + (1 / 12) * np.einsum("akc,bcd,abkl->dl", f, f, T)
                X = (1j) ** 2 * Z
                M += X if pm == phys[0] else X.T
    return M

# ---------- file access ----------
def read(fn):
    lines = open(UP + fn).read().split("\n")
    return lines

def align_blocks(fn, keep_comments=False):
    """list of (first line number (1-based), text) for each align block not commented out"""
    lines = read(fn); out = []; cur = None
    for n, l in enumerate(lines, 1):
        if l.lstrip().startswith("%") and not keep_comments: continue
        if r"\begin{align}" in l and cur is None:
            cur = [n, l.split(r"\begin{align}", 1)[1]]
            if r"\end{align}" in cur[1]:
                cur[1] = cur[1].split(r"\end{align}")[0]; out.append(tuple(cur)); cur = None
            continue
        if cur is not None:
            if r"\end{align}" in l:
                cur[1] += "\n" + l.split(r"\end{align}")[0]; out.append(tuple(cur)); cur = None
            else:
                cur[1] += "\n" + l
    return out

def list_terms(blk):
    """monomial strings of a list block:  ... \\Bigg[ term \\notag\\\\ term ... \\Bigg]"""
    a = blk.index(r"\Bigg[") + len(r"\Bigg[")
    b = blk.rindex(r"\Bigg]") if r"\Bigg]" in blk[a:] else len(blk)
    body = blk[a:b]
    parts = re.split(r"\\(?:notag|nonumber)\s*\\\\", body)
    return [p for p in (q.strip() for q in parts) if p.replace("&", "").strip()]

def list_monos(blk):
    res = []
    for t in list_terms(blk): res += parse_mono(t)
    return res

def header_factor(blk):
    """numeric prefactor of a list block in units of  g^4/(16 pi^5) / (2pi)^3 / k^+ * log ; returns (factor, has_i)"""
    a = blk.index("=") + 1
    b = blk.index(r"\int")
    h = blk[a:b].replace(r"^{\prime}", "'")
    h = re.sub(r"\{\\color\{blue\}.*", "", h, flags=re.S) if r"{\color{blue}" in h else h
    h = h.replace(r"\,", "").replace("&", "").replace(" ", "").replace("\n", "")
    fac = 1.0; has_i = False
    pats = [(r"\frac{1}{(2\pi)^{3}}", 1, False), (r"\frac{i}{(2\pi)^{3}}", 1, True), (r"\frac{g^{4}}{16\pi^{5}}", 1, False),
            (r"\frac{ig^{4}}{16\pi^{5}}", 1, True), (r"\frac{g^{4}}{8\pi^{3}}\frac{1}{2\pi^{2}}", 1, False),
            (r"\frac{ig^{4}}{8\pi^{5}}", 2, True), (r"\frac{g^{4}}{8\pi^{5}}", 2, False), (r"\frac{1}{k^{+}}", 1, False)]
    for p, v, ii in pats:
        if p in h:
            assert h.count(p) == 1, (p, h)
            h = h.replace(p, ""); fac *= v; has_i = has_i or ii
    if h.startswith("-"): fac = -fac; h = h[1:]
    if h.startswith("+"): h = h[1:]
    assert h == "", ("header leftover", h)
    return fac, has_i

def list_monos_h(blk, want_i=None):
    fac, has_i = header_factor(blk)
    if want_i is not None: assert has_i == want_i, "i mismatch in header"
    res = list_monos(blk)
    for m in res: m.coef *= fac
    return res
