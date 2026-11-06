"""Diff the (orientation-fixed) 4rho master's 4rho sector against the printed
collected 2-Wilson (line 72) and 4-Wilson (lines 76-78) sums, with positions FIXED
(no relabeling), to localize residual mis-transcriptions in the master equation.
"""
import itertools
from engine import Term, canonical_order
from symcolor import colour_canonical
from run_test1 import sym_hc
from texparse import parse_block
from blocks import block_body, F3, PI
import scan4rho as s

RHO2W = r"\rho^{b^{\prime}}(x^{\prime})\rho^{b}(u)\rho^{c}(x)\rho^{c}(y)"
RHO4W = r"\rho^{b^{\prime}}(x^{\prime})\rho^{d}(x)\rho^{c}(u)\rho^{e}(y)"

F4C = "/home/user/correlations/gluon-production-nlo/corrected/4rho_corrected.tex"

def grab(l0, l1):
    with open(F4C) as f:
        lines = f.read().split("\n")
    return "\n".join(lines[l0 - 1:l1])

def extract_brackets(text):
    """return list of \bigg[...\bigg] bodies (balanced on \bigg[ / \bigg])"""
    out = []
    i = 0
    while True:
        j = text.find(r"\bigg[", i)
        if j < 0:
            break
        k = text.find(r"\bigg]", j)
        out.append(text[j + len(r"\bigg["):k])
        i = k + 1
    return out

def parse_sum(bodies, rho_tex, coef, tag):
    terms = []
    for n, b in enumerate(bodies):
        ts = parse_block(b + rho_tex, extra_coef=coef, tag=f"{tag}{n}")
        # hack: append rho product to each term instead
        terms.extend(ts)
    return terms

def split_raw(body):
    from texparse import clean
    body = clean(body)
    terms, cur, depth = [], "", 0
    for c in body:
        if c in "{(":
            depth += 1
        elif c in ")}":
            depth -= 1
        if c in "+-" and depth == 0:
            if cur.strip():
                terms.append(cur)
            cur = c
        else:
            cur += c
    if cur.strip():
        terms.append(cur)
    return [t for t in terms if "U^" in t]

def parse_sum2(bodies, rho_tex, coef, tag):
    """append rho product per TERM (split first)"""
    from texparse import parse_term
    terms = []
    for n, b in enumerate(bodies):
        for i, tt in enumerate(split_raw(b)):
            terms.extend(parse_term(tt + " " + rho_tex, extra_coef=coef, tag=f"{tag}{n}#{i}"))
    return terms

def fixed_pos_ledger(terms, sector):
    """canonical order into CANON rho order; key = (kern, colour_canonical, rhos) with
    positions fixed."""
    led = {}
    for t0 in terms:
        for ct in canonical_order(t0):
            if len(ct.rhos) != sector:
                continue
            kkey = tuple(sorted([k for k in ct.kern]))
            # normalize rho dummies: rename rho indices in order of appearance
            ren = {}
            def R(ix):
                if ix not in ren:
                    ren[ix] = f"r{len(ren)}"
                return ren[ix]
            rhos = tuple((R(ix), v) for ix, v in ct.rhos)
            color = []
            for cf in ct.color:
                if cf[0] == 'U':
                    color.append(('U', ren.get(cf[1], cf[1]), ren.get(cf[2], cf[2]), cf[3]))
                elif cf[0] == 'f':
                    color.append(('f', *[ren.get(x, x) for x in cf[1:]]))
                else:
                    color.append(('delta', *[ren.get(x, x) for x in cf[1:]]))
            ckey, sgn = colour_canonical(Term(1.0, ct.kern, color, list(rhos), tag="k"))
            key = (kkey, ckey, rhos, ct.Nc_pow)
            led[key] = led.get(key, 0.0) + ct.coef * sgn
    return {k: v for k, v in led.items() if abs(v) > 1e-12}

def show_diff(A, B, nameA, nameB, limit=200):
    keys = set(A) | set(B)
    bad = [(k, A.get(k, 0.0), B.get(k, 0.0)) for k in keys
           if abs(A.get(k, 0.0) - B.get(k, 0.0)) > 1e-12]
    print(f"--- {nameA} vs {nameB}: {len(bad)} mismatched structures "
          f"(|A|={len(A)}, |B|={len(B)})")
    sc = 1 / (16 * PI**5)
    for k, a, b in sorted(bad, key=lambda x: -abs(x[1] - x[2]))[:limit]:
        kkey, ckey, rhos, ncp = k
        print(f"   A={a/sc:+.3f} B={b/sc:+.3f}  rhos={rhos}")
        print(f"      color={ckey}")
    return bad

if __name__ == "__main__":
    c = 1 / (16 * PI**5)
    # ---- printed collected sums ----
    b2w = extract_brackets(grab(71, 73))
    t2w = parse_sum2(b2w, RHO2W, c, "c2w")
    print("parsed 2W collected terms:", len(t2w))
    b4w = extract_brackets(grab(75, 79))
    print("4W brackets found:", len(b4w))
    t4w = parse_sum2(b4w, RHO4W, c, "c4w")
    print("parsed 4W collected terms:", len(t4w))
    for t in t2w + t4w:
        t.kern = list(s.KERN)

    # ---- master (orientation-fixed), 4rho sector, split by #Wilson lines ----
    fl = dict.fromkeys(s.FLIPS, False)
    for k in ['fB1a_z', 'fB1a_wp', 'fB1b_z', 'fB2a2_z', 'fB2a2_wp', 'fB2b1', 'fB3a', 'fB3b']:
        fl[k] = True
    for label, flags in (("FIXED", fl), ("ASPRINTED", dict.fromkeys(s.FLIPS, False))):
        m2, m4 = [], []
        for t, hc in s.build(flags):
            nw = len([cf for cf in t.color if cf[0] == 'U'])
            (m2 if nw == 2 else m4).append(t)
        # do the printed collected sums include the c.c. of groups A/C already?
        for withhc in (False, True):
            mm2 = m2 + ([sym_hc(t) for t in m2 if True] if withhc else [])
            # hc flag per term: A/C only
            mm2 = m2[:]
            mm4 = m4[:]
            if withhc:
                mm2 += [sym_hc(t) for t in m2 if t.tag.startswith(("4r_A", "4r_C"))]
                mm4 += [sym_hc(t) for t in m4 if t.tag.startswith(("4r_A", "4r_C"))]
            L2 = fixed_pos_ledger(mm2, 4)
            L4 = fixed_pos_ledger(mm4, 4)
            C2 = fixed_pos_ledger(t2w, 4)
            C4 = fixed_pos_ledger(t4w, 4)
            n2 = len(show_diff(L2, C2, f"{label} master2W(hc={withhc})", "printed2W", limit=0))
            n4 = len(show_diff(L4, C4, f"{label} master4W(hc={withhc})", "printed4W", limit=0))
            print(f"=== {label} hc={withhc}: 2W mismatches {n2}, 4W mismatches {n4}")
