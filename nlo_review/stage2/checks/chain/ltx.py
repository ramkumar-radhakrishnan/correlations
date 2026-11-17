# Minimal LaTeX -> numeric evaluator for the kernel/colour formulas of the NLO files.
# Value types: complex scalar | Vec (transverse 2-vector) | Op (list of terms)
#   term = (coef, tensors[(array, idxnames)], rhos[(idx, point)], tags)
import re, numpy as np
from scipy.linalg import expm
UP = "/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
FILES = dict(F="690ba62f-Final_terms.tex", L="56ca25e0-Limiting_Cases.tex", O="f10fa593-Other_terms.tex",
             G1="16e62d79-Group_I.tex", G2="24c35375-Group_II.tex", G3="96b5b029-Group_III.tex")
_cache = {}
def lines(f, a, b):
    if f not in _cache: _cache[f] = open(UP + FILES[f]).read().split("\n")
    return "\n".join(_cache[f][a - 1:b])

class Vec:
    def __init__(s, v): s.v = np.asarray(v, dtype=complex)
class Op:
    def __init__(s, terms): s.terms = terms
def is_s(v): return isinstance(v, (int, float, complex, np.number))
def promote(v):
    if isinstance(v, Op): return v
    if is_s(v): return Op([(complex(v), [], [], ())])
    raise TypeError("cannot promote Vec to Op")
def mul(a, b):
    if is_s(a) and is_s(b): return a * b
    if is_s(a) and isinstance(b, Vec): return Vec(a * b.v)
    if isinstance(a, Vec) and is_s(b): return Vec(b * a.v)
    if isinstance(a, Vec) or isinstance(b, Vec): raise TypeError("Vec times Vec/Op without \\cdot")
    a, b = promote(a), promote(b)
    return Op([(c1 * c2, t1 + t2, r1 + r2, g1 + g2) for (c1, t1, r1, g1) in a.terms for (c2, t2, r2, g2) in b.terms])
def add(a, b, sgn=1):
    if is_s(a) and is_s(b): return a + sgn * b
    if isinstance(a, Vec) and isinstance(b, Vec): return Vec(a.v + sgn * b.v)
    if isinstance(a, Vec) or isinstance(b, Vec): raise TypeError("Vec + non-Vec")
    a, b = promote(a), promote(b)
    return Op(a.terms + [(sgn * c, t, r, g) for (c, t, r, g) in b.terms])

def make_env(N=2, seed=0, pts=("x", "y", "z", "w", "x'", "y'", "w'")):
    rng = np.random.default_rng(seed)
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
    def adj():
        H = sum(rng.normal() * t[a] for a in range(n)); V = expm(1j * H)
        return np.real(np.array([[2 * np.trace(t[a] @ V @ t[b] @ V.conj().T) for b in range(n)] for a in range(n)]))
    env = dict(N=N, n=n, f=f, U={p: adj() for p in pts}, X={p: rng.normal(size=2) for p in pts},
               k=rng.normal(size=2), kp=float(rng.uniform(0.7, 1.9)), pp=float(rng.uniform(0.1, 0.6)),
               vee=1e300, Lam=1e-300, dperp=2.0, sym={}, g=1.3, logmask=None, logs_seen=set())
    env["qp"] = env["pp"]
    return env

SKIP = [r"\notag", r"\nonumber", r"\left", r"\right", r"\Bigg", r"\bigg", r"\Big", r"\big", r"\,", r"\;", r"\!",
        r"\qquad", r"\quad", r"\displaystyle", "&", "~"]
def normalize(s):
    s = s.replace(r"^{\prime}", "'").replace(r"^\prime", "'")
    s = re.sub(r"\\\\(\[[^\]]*\])?", " ", s)
    for k in SKIP: s = s.replace(k, " ")
    s = re.sub(r"\\mathrm\{arc\s*cos\}", r"\\arccos ", s)
    s = re.sub(r"\\mathrm\{arc\s*tan\}", r"\\arctan ", s)
    s = re.sub(r"ik\s*\(", r"ik\\cdot(", s)
    s = re.sub(r"\{\\color\{[a-z]+\}", "{", s)
    s = s.replace(r"\tfrac", r"\frac").replace(r"\rm", " ")
    s = s.replace(r"\mu^{2}_{\overline{MS}}", r"\musq ")
    s = re.sub(r"\\widetilde\{\\mathcal\s*\{?P\}?\}_\{gg\}\(\\xi\)", r"\\Pgg ", s)
    return s
TOK = re.compile(r"\\[a-zA-Z]+|\\[{}]|[0-9]+(?:\.[0-9]+)?|[a-zA-Z]'*|[{}()\[\]^_+\-,=.]")
TIDX = set("ijmkln")
class P:
    def __init__(s, text, env):
        s.env = env; s.src = normalize(text); s.t = TOK.findall(s.src); s.i = 0
        rest = TOK.sub("", s.src).strip()
        if rest: raise SyntaxError("untokenized: " + rest[:80])
    def peek(s, k=0): return s.t[s.i + k] if s.i + k < len(s.t) else None
    def nxt(s):
        tok = s.peek(); s.i += 1; return tok
    def expect(s, tok):
        got = s.nxt()
        if got != tok: raise SyntaxError(f"expected {tok} got {got} at token {s.i}: ...{' '.join(s.t[max(0,s.i-12):s.i+6])}")
    def raw_group(s):
        s.expect("{"); d = 1; out = []
        while d:
            tok = s.nxt()
            if tok is None: raise SyntaxError("unbalanced {")
            if tok == "{": d += 1
            elif tok == "}": d -= 1
            if d: out.append(tok)
        return out
    def group(s):
        s.expect("{"); v = s.expr(); s.expect("}"); return v
    # ---------- grammar
    def expr(s):
        sgn = 1
        if s.peek() in ("+", "-"): sgn = 1 if s.nxt() == "+" else -1
        v = s.product()
        if sgn < 0: v = mul(-1, v)
        while s.peek() in ("+", "-"):
            sg = 1 if s.nxt() == "+" else -1
            v = add(v, s.product(), sg)
        return v
    ENDS = ("+", "-", ")", "]", "}", r"\}", ",", None, "=", ".")
    def product(s):
        acc = 1.0; pend = None; awaiting = False; n = 0
        while True:
            tok = s.peek()
            if tok == r"\times": s.nxt(); continue
            if tok == r"\cdot":
                s.nxt()
                if pend is None: raise SyntaxError("\\cdot without vector on the left")
                awaiting = True; continue
            if tok in s.ENDS: break
            F = s.factor(); n += 1
            if isinstance(F, Vec):
                if pend is None: pend = F
                elif awaiting: acc = mul(acc, complex((pend.v * F.v).sum())); pend = None; awaiting = False
                else: raise SyntaxError("two vectors juxtaposed without \\cdot near: " + " ".join(s.t[max(0, s.i - 10):s.i]))
            elif awaiting and is_s(F): pend = Vec(pend.v * F)
            elif awaiting: raise SyntaxError("\\cdot followed by non-vector")
            else: acc = mul(acc, F)
        if n == 0: raise SyntaxError("empty product near: " + " ".join(s.t[max(0, s.i - 10):s.i + 5]))
        if awaiting: raise SyntaxError("dangling \\cdot")
        if pend is not None: return mul(acc, pend)
        return acc
    def idxs(s):
        return s.raw_group()
    def factor(s):
        v = s.atom()
        while s.peek() == "^":
            s.nxt(); g = s.raw_group()
            if isinstance(v, Vec):
                if len(g) == 1 and g[0].isdigit():
                    p = int(g[0]); assert p % 2 == 0; v = complex((v.v * v.v).sum()) ** (p // 2)
                elif len(g) == 1 and g[0][0] in TIDX: v = Op([(1.0, [(v.v, ("t" + g[0],))], [], ())])
                else: raise SyntaxError("bad vector superscript " + str(g))
            elif is_s(v):
                if len(g) == 1 and g[0].isdigit(): v = v ** int(g[0])
                else: raise SyntaxError("bad scalar superscript " + str(g))
            else: raise SyntaxError("superscript on operator")
        return v
    def bracket(s, close):
        i0 = s.i; v = s.expr(); content = "".join(s.t[i0:s.i]); s.expect(close)
        if s.peek() == "_" and s.peek(1) == "{" and s.peek(2) == "+" and s.peek(3) == "}":
            s.i += 4
            var = "q" if "q^{+}" in content else "p"
            pole = "k" if ("k^{+}-" + var + "^{+}") in content else "0"
            if not is_s(v): raise SyntaxError("plus bracket with non-scalar")
            return Op([(complex(v), [], [], (f"plus[{var}@{pole}]",))])
        return v
    def logarg(s):
        tok = s.peek()
        if tok == "{": return s.group()
        if tok in ("(", "["): s.nxt(); return s.bracket(")" if tok == "(" else "]")
        if tok == r"\frac": return s.atom()
        raise SyntaxError("log argument")
    def atom(s):
        e = s.env; tok = s.nxt()
        if tok is None: raise SyntaxError("unexpected end")
        if re.fullmatch(r"[0-9.]+", tok): return float(tok)
        if tok == r"\frac":
            a = s.group(); b = s.group()
            if not is_s(b): raise SyntaxError("non-scalar denominator")
            return mul(a, 1.0 / b)
        if tok == "(": return s.bracket(")")
        if tok == "[": return s.bracket("]")
        if tok == "{": s.i -= 1; return s.group()
        if tok == r"\{":
            a = s.expr()
            if s.peek() == ",":
                s.nxt(); b = s.expr(); s.expect(r"\}")
                return add(mul(a, b), mul(b, a))
            s.expect(r"\}"); return a
        if tok == r"\pi": return np.pi
        if tok == r"\vee": return e["vee"]
        if tok == r"\Lambda": return e["Lam"]
        if tok == r"\epsilon": return e["sym"]["eps"]
        if tok == r"\musq": return e["sym"]["musq"]
        if tok == r"\xi": return e["sym"]["xi"]
        if tok == r"\gamma":
            if s.peek() == "_": g = (s.nxt(), s.raw_group()); assert g == ("_", ["E"])
            return float(np.euler_gamma)
        if tok == r"\beta":
            g = (s.nxt(), s.raw_group()); assert g == ("_", ["0"]); return e["sym"]["beta0"]
        if tok == r"\Pgg": return Op([(1.0, [], [], ("Pgg",))])
        if tok == r"\mathcal":
            g = s.raw_group()
            if g == ["K"]:
                s.expect("("); v = s.expr(); s.expect(")"); assert isinstance(v, Vec)
                return Vec(v.v / (v.v * v.v).sum())
            if g == ["A"]:
                assert (s.nxt(), s.raw_group()) == ("_", ["0"]); return e["sym"]["A0"]
            raise SyntaxError("mathcal " + str(g))
        if tok == r"\log":
            p = 1
            if s.peek() == "^": s.nxt(); p = int(s.raw_group()[0])
            i0 = s.i; a = s.logarg(); name = "".join(s.t[i0:s.i])
            if not is_s(a): raise SyntaxError("log of non-scalar")
            e["logs_seen"].add(name)
            if abs(complex(a).imag) > 1e-12 * abs(a): raise SyntaxError("log of complex")
            val = float(np.log(abs(complex(a).real)))
            if e["logmask"] is not None: val = 1.0 if name == e["logmask"] else 0.0
            return val ** p
        if tok == r"\sqrt":
            a = s.group(); return complex(np.sqrt(complex(a)))
        if tok in (r"\arccos", r"\arctan"):
            a = s.atom(); a = complex(a); assert abs(a.imag) < 1e-12
            return float(np.arccos(a.real)) if tok == r"\arccos" else float(np.arctan(a.real))
        if tok == "e" and s.peek() == "^":
            s.nxt(); a = s.group(); return complex(np.exp(complex(a)))
        if tok == "U":
            s.expect("^"); ix = s.idxs(); assert len(ix) == 2, ix
            s.expect("("); pt = s.nxt(); s.expect(")")
            return Op([(1.0, [(e["U"][pt], tuple(ix))], [], ())])
        if tok == "f" and s.peek() == "^":
            s.nxt(); ix = s.idxs(); assert len(ix) == 3, ix
            return Op([(1.0, [(e["f"], tuple(ix))], [], ())])
        if tok == r"\rho":
            s.expect("^"); ix = s.idxs(); assert len(ix) == 1, ix
            s.expect("("); pt = s.nxt(); s.expect(")")
            assert pt in e["X"], pt
            return Op([(1.0, [], [(ix[0], pt)], ())])
        if tok == r"\delta":
            assert s.nxt() in "^_"; ix = s.idxs(); assert len(ix) == 2 and all(c[0] in TIDX for c in ix), ix
            return Op([(1.0, [(np.eye(2), tuple("t" + c for c in ix))], [], ())])
        if tok == "N" and s.peek() == "_":
            s.nxt(); assert s.raw_group() == ["c"]; return float(e["N"])
        if tok == "d":
            if s.peek() == "_": s.nxt(); assert s.raw_group() == [r"\perp"]; return e["dperp"]
            if s.peek() in ("p", "q"): s.nxt(); s.expect("^"); assert s.raw_group() == ["+"]; return 1.0
            if s.peek() == r"\xi": s.nxt(); return 1.0
            raise SyntaxError("bare d")
        if tok in ("k", "p", "q") and s.peek() == "^" and s.peek(1) == "{" and s.peek(2) == "+":
            s.nxt(); g = s.raw_group(); val = e[{"k": "kp", "p": "pp", "q": "qp"}[tok]]
            if g == ["+"]: return val
            if len(g) == 2 and g[1].isdigit(): return val ** int(g[1])
            raise SyntaxError("bad plus-momentum power " + str(g))
        if tok == "k": return Vec(e["k"])
        if tok == "g": return e["g"]
        if tok == "i": return 1j
        if tok in ("X", "Z") and s.peek() == "_":
            s.nxt(); assert s.raw_group() == [r"\perp"]; return Vec(e["sym"][tok + "perp"])
        if tok in e["X"]: return Vec(e["X"][tok])
        if tok == r"\int":
            while s.peek() in ("_", "^"): s.nxt(); s.raw_group()
            return 1.0
        raise SyntaxError(f"unknown token {tok} near: " + " ".join(s.t[max(0, s.i - 10):s.i + 5]))

def parse(text, env, tokmut=None):
    p = P(text, env)
    if tokmut: p.t = tokmut(p.t)
    v = p.expr()
    if p.peek() == ".": p.nxt()
    if p.peek() is not None: raise SyntaxError("trailing tokens: " + " ".join(p.t[p.i:p.i + 12]))
    return v
def rhs(text):
    """drop the left-hand side: everything up to the first '=' """
    return text[text.index("=") + 1:]
LET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
def contract(op, tagmap=None):
    """-> dict (rho points in written order, tags) -> coefficient tensor (axes in rho order)"""
    out = {}
    for coef, tens, rhos, tags in promote(op).terms:
        names = {}; cnt = {}
        for arr, ix in tens:
            for c in ix: cnt[c] = cnt.get(c, 0) + 1
        for c, pt in rhos: cnt[c] = cnt.get(c, 0) + 1
        bad = {c: m for c, m in cnt.items() if m != 2}
        if bad: raise ValueError(f"index multiplicity != 2: {bad}; rhos={rhos}; tensors={[ix for _, ix in tens]}")
        for c in cnt: names[c] = LET[len(names)]
        sub = ",".join("".join(names[c] for c in ix) for arr, ix in tens) + "->" + "".join(names[c] for c, pt in rhos)
        val = coef * (np.einsum(sub, *[arr for arr, ix in tens], optimize="greedy") if tens else 1.0)
        tg = tuple(sorted((tagmap or {}).get(t, t) for t in tags))
        key = (tuple(pt for c, pt in rhos), tg)
        out[key] = out.get(key, 0) + val
    return out
def compare(A, B):
    """returns (max |A-B|, max |A|, max|B|, list of per-key info)"""
    info = []; md = 0; ma = 0; mb = 0
    for key in sorted(set(A) | set(B), key=str):
        a = A.get(key, 0); b = B.get(key, 0)
        d = float(np.max(np.abs(np.asarray(a) - np.asarray(b)))); na = float(np.max(np.abs(a))); nb = float(np.max(np.abs(b)))
        if not np.isfinite(d): d = float("inf")
        md = max(md, d); ma = max(ma, na); mb = max(mb, nb); info.append((key, d, na, nb))
    return md, ma, mb, info
def report(name, A, B, verbose=False):
    md, ma, mb, info = compare(A, B)
    rel = md / max(ma, mb, 1e-300) if np.isfinite(ma) and np.isfinite(mb) else float("inf")
    print(f"{name}: maxdiff={md:.3e} |A|={ma:.3e} |B|={mb:.3e} rel={rel:.2e} {'AGREE' if rel < 1e-9 else 'MISMATCH'}")
    if verbose or rel >= 1e-9:
        for key, d, na, nb in info: print(f"    {key}: diff={d:.3e} |A|={na:.3e} |B|={nb:.3e}")
    return rel
