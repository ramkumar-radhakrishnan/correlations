"""Symbolic canonicalization of kernel+colour structures (3-rho sector ledger)."""
import itertools

POSORD = {"x": 0, "y": 1, "z": 2, "u": 3, "w": 4, "wp": 5, "xp": 6, "yp": 7}

def canon_kernel(term):
    sign = 1.0
    pairs = []
    for kf in term.kern:
        if kf[0] == 'K':
            a, b, c, d = kf[1:5]
            if POSORD[a] > POSORD[b]:
                a, b = b, a; sign = -sign
            if POSORD[c] > POSORD[d]:
                c, d = d, c; sign = -sign
            pairs.append(tuple(sorted([(a, b), (c, d)])))
        elif kf[0] == 'inv2':
            a, b = sorted(kf[1:3], key=lambda v: POSORD[v])
            pairs.append((('inv2', a, b),))
        else:
            pairs.append(('fun', id(kf[1])))
    return tuple(sorted(pairs)), sign

def colour_canonical(term):
    fs = [tuple(c[1:]) for c in term.color if c[0] == 'f']
    Us = [tuple(c[1:]) for c in term.color if c[0] == 'U']
    ds = [tuple(c[1:]) for c in term.color if c[0] == 'delta']
    rhos = sorted(term.rhos, key=lambda r: POSORD[r[1]])
    idxs = set()
    for t in fs: idxs.update(t)
    for t in Us: idxs.update(t[:2])
    for t in ds: idxs.update(t)
    for i, _ in rhos: idxs.add(i)
    idxs = sorted(idxs)

    def sig(ix):
        s = []
        for (a, b, c) in fs:
            s.append(("f", (ix == a) + (ix == b) + (ix == c)))
        for (a, b, p) in sorted(Us, key=lambda u: (u[2],)):
            s.append(("U", p, "L" if ix == a else ("R" if ix == b else "-")))
        for (a, b) in ds:
            s.append(("d", (ix == a) + (ix == b)))
        for (i, p) in rhos:
            s.append(("r", p, ix == i))
        return tuple(s)

    groups = {}
    for ix in idxs:
        groups.setdefault(str(sig(ix)), []).append(ix)
    gkeys = sorted(groups)
    base = {}
    cnt = 0
    for gk in gkeys:
        base[gk] = cnt
        cnt += len(groups[gk])

    def build(m):
        sign = 1.0
        parts = []
        for (a, b, c) in fs:
            tt = [m[a], m[b], m[c]]
            s = 1
            for i in range(3):
                for j in range(2 - i):
                    if tt[j] > tt[j + 1]:
                        tt[j], tt[j + 1] = tt[j + 1], tt[j]
                        s = -s
            if tt[0] == tt[1] or tt[1] == tt[2]:
                return None, 0.0
            sign *= s
            parts.append(("f",) + tuple(tt))
        for (a, b, p) in Us:
            parts.append(("U", m[a], m[b], p))
        for (a, b) in ds:
            parts.append(("d",) + tuple(sorted((m[a], m[b]))))
        for (i, p) in rhos:
            parts.append(("r", m[i], p))
        return tuple(sorted(parts)), sign

    best_key, best_sign, zero = None, None, False
    gperm_lists = [list(itertools.permutations(groups[gk])) for gk in gkeys]
    for combo in itertools.product(*gperm_lists):
        m = {}
        for gk, perm in zip(gkeys, combo):
            for slot, ix in enumerate(perm):
                m[ix] = f"i{base[gk] + slot:02d}"
        key, sign = build(m)
        if key is None:
            return "ZERO", 0.0
        if best_key is None or key < best_key:
            best_key, best_sign = key, sign
        elif key == best_key and sign != best_sign:
            zero = True
    if zero:
        return "SELFZERO", 0.0
    return best_key, best_sign

def ledger_report(terms, title, tol=1e-12):
    led = {}
    for t in terms:
        kkey, ksign = canon_kernel(t)
        ckey, csign = colour_canonical(t)
        if ckey in ("ZERO", "SELFZERO"):
            continue
        key = (kkey, ckey, t.Nc_pow)
        c, tags = led.get(key, (0.0, []))
        led[key] = (c + t.coef * ksign * csign, tags + [t.tag])
    surv = {k: v for k, v in led.items() if abs(v[0]) > tol}
    print(f"=== {title}: {len(led)} distinct structures, {len(surv)} survive ===")
    return led, surv

def pretty(key, val):
    kkey, ckey, ncp = key
    ks = []
    for p in kkey:
        if len(p) == 2 and isinstance(p[0], tuple):
            ks.append(f"K({p[0][0]}-{p[0][1]}).K({p[1][0]}-{p[1][1]})")
        else:
            ks.append(str(p))
    cs = " ".join(("%s^{%s%s}(%s)" % (x[0], x[1], x[2], x[3]) if x[0] == "U" else
                   ("f^{%s%s%s}" % x[1:] if x[0] == "f" else
                    ("rho^{%s}(%s)" % (x[1], x[2]) if x[0] == "r" else str(x))))
                  for x in ckey)
    return f"{val[0]:+.6g} * Nc^{ncp} * {' '.join(ks)} * {cs}   <= {sorted(set(val[1]))}"

# ---------- full canonicalization incl. internal-position relabeling ----------
import itertools as _it
from engine import substitute_var as _subs

_INTERNAL = ["x", "y", "z", "u", "xp", "yp"]

def canon_full(term):
    """min over permutations of the internal position labels present in the term."""
    pres = [v for v in _INTERNAL if v in term.vars()]
    best = None
    zero = False
    for perm in _it.permutations(pres):
        t = term
        # apply substitution via temp names
        for v in pres:
            t = _subs(t, v, "TMP_" + v)
        for v, tgt in zip(pres, perm):
            t = _subs(t, "TMP_" + v, tgt)
        kkey, ks = canon_kernel(t)
        ckey, cs = colour_canonical(t)
        if ckey == "ZERO":
            return ("ZERO", None), 0.0
        if ckey == "SELFZERO":
            zero = True
            continue
        cand = (kkey, ckey)
        s = ks * cs
        if best is None or cand < best[0]:
            best = (cand, s)
        elif cand == best[0] and s != best[1]:
            zero = True
    if zero and best is None:
        return ("SELFZERO", None), 0.0
    if zero:
        return ("SELFZERO", None), 0.0
    return best

def ledger_full(terms, title, tol=1e-12):
    led = {}
    for t in terms:
        (key, s) = canon_full(t)
        if key[0] in ("ZERO", "SELFZERO"):
            continue
        k = (key, t.Nc_pow)
        c, tags = led.get(k, (0.0, []))
        led[k] = (c + t.coef * s, tags + [t.tag])
    surv = {k: v for k, v in led.items() if abs(v[0]) > tol}
    print(f"=== {title}: {len(led)} distinct structures, {len(surv)} survive ===")
    return led, surv

def pretty_full(key, val):
    (kkey, ckey), ncp = key
    ks = []
    for p in kkey:
        if len(p) == 2 and isinstance(p[0], tuple):
            ks.append(f"K({p[0][0]}-{p[0][1]}).K({p[1][0]}-{p[1][1]})")
        else:
            ks.append(str(p))
    cs = " ".join(("U^{%s%s}(%s)" % (x[1], x[2], x[3]) if x[0] == "U" else
                   ("f^{%s%s%s}" % x[1:] if x[0] == "f" else
                    ("rho^{%s}(%s)" % (x[1], x[2]) if x[0] == "r" else str(x))))
                  for x in ckey)
    return f"{val[0]:+.6g} * Nc^{ncp} * {' '.join(ks)} * {cs}   <= {sorted(set(val[1]))}"

# ---------- algebraic normalization: delta absorption, same-site U-orthogonality ----------
from engine import Term as _Term

def normalize_term(t):
    """Returns (Term, extra_dimA_power) with deltas absorbed, same-site orthogonal
    U-pairs reduced, f with repeated indices -> None (zero). dimA = Nc^2-1 factors
    from self-contracted deltas."""
    coef = t.coef
    color = list(t.color)
    rhos = list(t.rhos)
    dimA = 0
    changed = True
    while changed:
        changed = False
        # absorb deltas
        for n, c in enumerate(color):
            if c[0] == 'delta':
                i, j = c[1], c[2]
                color.pop(n)
                if i == j:
                    dimA += 1
                else:
                    def ren(ix):
                        return i if ix == j else ix
                    color = [(cc[0], ren(cc[1]), ren(cc[2]), cc[3]) if cc[0] == 'U'
                             else (('f', ren(cc[1]), ren(cc[2]), ren(cc[3])) if cc[0] == 'f'
                                   else ('delta', ren(cc[1]), ren(cc[2])))
                             for cc in color]
                    rhos = [(ren(ix), p) for ix, p in rhos]
                changed = True
                break
        if changed:
            continue
        # f with repeated index -> zero
        for c in color:
            if c[0] == 'f' and len({c[1], c[2], c[3]}) < 3:
                return None, 0
        # same-site U orthogonality: U^{ab}(s) U^{ac}(s) -> delta^{bc}; U^{ba}U^{ca} -> delta^{bc}
        nU = len(color)
        done = False
        for n1 in range(nU):
            if color[n1][0] != 'U':
                continue
            for n2 in range(n1 + 1, nU):
                if color[n2][0] != 'U':
                    continue
                _, a1, b1, p1 = color[n1]
                _, a2, b2, p2 = color[n2]
                if p1 != p2:
                    continue
                if a1 == a2 and _count_index(color, rhos, a1) == 2:
                    newc = [c for k, c in enumerate(color) if k not in (n1, n2)]
                    newc.append(('delta', b1, b2))
                    color = newc; done = True; break
                if b1 == b2 and _count_index(color, rhos, b1) == 2:
                    newc = [c for k, c in enumerate(color) if k not in (n1, n2)]
                    newc.append(('delta', a1, a2))
                    color = newc; done = True; break
            if done:
                break
        if done:
            changed = True
    t2 = _Term(coef, t.kern, color, rhos, phase=t.phase, Nc_pow=t.Nc_pow, tag=t.tag)
    return t2, dimA

def _count_index(color, rhos, ix):
    n = 0
    for c in color:
        if c[0] == 'U':
            n += (c[1] == ix) + (c[2] == ix)
        elif c[0] == 'f':
            n += (c[1] == ix) + (c[2] == ix) + (c[3] == ix)
        elif c[0] == 'delta':
            n += (c[1] == ix) + (c[2] == ix)
    for i, _ in rhos:
        n += (i == ix)
    return n

def ledger_full2(terms, title, tol=1e-12):
    led = {}
    for t0 in terms:
        t, dA = normalize_term(t0)
        if t is None:
            continue
        (key, s) = canon_full(t)
        if key[0] in ("ZERO", "SELFZERO"):
            continue
        k = (key, t.Nc_pow, dA)
        c, tags = led.get(k, (0.0, []))
        led[k] = (c + t.coef * s, tags + [t.tag])
    surv = {k: v for k, v in led.items() if abs(v[0]) > tol}
    print(f"=== {title}: {len(led)} distinct structures, {len(surv)} survive ===")
    return led, surv

def pretty_full2(key, val):
    (kkey, ckey), ncp, dA = key
    return pretty_full(((kkey, ckey), ncp), val) + (f"  [dimA^{dA}]" if dA else "")

def normalize_term2(t):
    """normalize_term + same-site f-invariance: f^{a'b'c'} U^{aa'}(s) U^{bb'}(s) -> f^{abc'} U^{c c'}(s)
    wait, correct identity: f^{a'b'c'} U^{aa'} U^{bb'} = f^{abc} U^{cc'} ... implemented as:
    if an f has two of its indices contracted with the PRIMED (second) slots of two U's at
    the same site: f^{p q r} U^{a p}(s) U^{b q}(s) = f^{a b c} U^{c r}(s); similarly for
    first slots: f^{p q r} U^{p a}(s) U^{q b}(s) = f^{a b c} U^{r c}(s)."""
    t, dA = normalize_term(t)
    if t is None:
        return None, 0
    changed = True
    while changed:
        changed = False
        color = list(t.color)
        fs = [(n, c) for n, c in enumerate(color) if c[0] == 'f']
        Us = [(n, c) for n, c in enumerate(color) if c[0] == 'U']
        for nf, f in fs:
            fid = list(f[1:])
            # find U pairs at same site contracted with two f indices on the same slot
            for slot in (1, 2):
                hits = []
                for nu, u in Us:
                    uin = u[slot]      # index on this slot
                    if uin in fid and _count_index(color, t.rhos, uin) == 2:
                        hits.append((nu, u, uin))
                # group by site
                bysite = {}
                for nu, u, uin in hits:
                    bysite.setdefault(u[3], []).append((nu, u, uin))
                done = False
                for site, lst in bysite.items():
                    if len(lst) >= 2:
                        (n1, u1, i1), (n2, u2, i2) = lst[0], lst[1]
                        if i1 == i2:
                            continue
                        # f^{..i1..i2..} U^{x i1}(s) U^{y i2}(s)  (slot=2)
                        #   = f with i1->x', i2->y' ... new index for remaining contraction
                        rest = [ix for ix in fid if ix not in (i1, i2)][0]
                        import itertools as _i2
                        # sign to bring f to order (i1, i2, rest)
                        perm = [fid.index(i1), fid.index(i2), fid.index(rest)]
                        sgn = 1
                        pp = perm[:]
                        for a_ in range(3):
                            for b_ in range(2 - a_):
                                if pp[b_] > pp[b_ + 1]:
                                    pp[b_], pp[b_ + 1] = pp[b_ + 1], pp[b_], 
                                    sgn = -sgn
                        other1 = u1[1] if slot == 2 else u1[2]
                        other2 = u2[1] if slot == 2 else u2[2]
                        newidx = "NF0"
                        k = 0
                        while _count_index(color, t.rhos, newidx) > 0:
                            k += 1
                            newidx = f"NF{k}"
                        newcolor = [c for n, c in enumerate(color) if n not in (nf, n1, n2)]
                        newcolor.append(('f', other1, other2, newidx))
                        if slot == 2:
                            newcolor.append(('U', newidx, rest, site))
                        else:
                            newcolor.append(('U', rest, newidx, site))
                        from engine import Term as _T
                        t = _T(sgn * t.coef, t.kern, newcolor, t.rhos, phase=t.phase,
                               Nc_pow=t.Nc_pow, tag=t.tag)
                        t, dA2 = normalize_term(t)
                        if t is None:
                            return None, 0
                        dA += dA2
                        changed = True
                        done = True
                        break
                if done:
                    break
            if changed:
                break
    return t, dA

def ledger_full3(terms, title, tol=1e-12, verbose_zero=False):
    led = {}
    for t0 in terms:
        t, dA = normalize_term2(t0)
        if t is None:
            continue
        (key, s) = canon_full(t)
        if key[0] in ("ZERO", "SELFZERO"):
            continue
        k = (key, t.Nc_pow, dA)
        c, tags = led.get(k, (0.0, []))
        led[k] = (c + t.coef * s, tags + [t.tag])
    surv = {k: v for k, v in led.items() if abs(v[0]) > tol}
    print(f"=== {title}: {len(led)} distinct structures, {len(surv)} survive ===")
    return led, surv
