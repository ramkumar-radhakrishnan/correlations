"""Scan Wilson-line orientation variants of the 4rho master equation (groups A,B)
so that (i) its 4rho sector cancels symbolically and (ii) its 3rho contact remainder
equals the collection's "From four rho's" blocks. Group C is held fixed (user-corrected).

Per-term canonical-ledger contributions are cached, so the 2^N combo scan is cheap.
"""
import itertools, math
from engine import Term, canonical_order
from symcolor import normalize_term2, canon_full
from run_test1 import sym_hc
from texparse import parse_block
from blocks import block_body, F3, PI

K = lambda a, b, c, d: ('K', a, b, c, d)
U = lambda i, j, p: ('U', i, j, p)
KERN = [K("xp", "wp", "u", "w"), K("x", "z", "y", "z")]

def dedup(color, rhos):
    out_r, extra, used = [], [], {}
    for n, (ix, p) in enumerate(rhos):
        if ix in used:
            nix = ix + f"_{n}"
            extra.append(('delta', ix, nix))
            out_r.append((nix, p))
        else:
            used[ix] = True
            out_r.append((ix, p))
    return color + extra, out_r

def prod(coef, br1_list, br2_list, tag):
    out = []
    for s1, c1, r1 in br1_list:
        for s2, c2, r2 in br2_list:
            color, rhos = dedup(c1 + c2, r1 + r2)
            out.append(Term(coef * s1 * s2, KERN, color, rhos, tag=tag))
    return out

def anti(i1, p1, i2, p2):
    return [[(i1, p1), (i2, p2)], [(i2, p2), (i1, p1)]]

FLIPS = ["fA3", "fA4", "fB1a_z", "fB1a_wp", "fB1b_z", "fB2a1", "fB2a2_z",
         "fB2a2_wp", "fB2b1", "fB2b2", "fB3a", "fB3b"]

def build(fl):
    """fl: dict flipname -> bool. Returns list of (Term, hc_flag)."""
    def Uo(i, j, p, flip):
        return U(j, i, p) if flip else U(i, j, p)
    T = []
    c = 1 / (16 * PI**5)
    BR1 = [(+1, [U("a", "bp", "xp")], [("bp", "xp")]),
           (-1, [U("a", "bp", "wp")], [("bp", "xp")])]
    A1b = [(+1, [U("a", "b", "u")], [("b", "u"), ("c", "x"), ("c", "y")]),
           (-1, [U("a", "b", "w"), U("c", "d", "x"), U("c", "e", "y")], [("d", "x"), ("e", "y"), ("b", "u")])]
    T += prod(-c, BR1, A1b, "4r_A1")
    A2b = [(+1, [U("a", "b", "w")], [("c", "x"), ("c", "y"), ("b", "u")]),
           (-1, [U("c", "d", "x"), U("c", "e", "y"), U("a", "b", "u")], [("d", "x"), ("e", "y"), ("b", "u")])]
    T += prod(+c, BR1, A2b, "4r_A2")
    A3b = []
    for rr in anti("d", "x", "c", "u"):
        A3b.append((+1, [U("a", "c", "w"), Uo("b", "d", "z", fl["fA3"]), U("b", "e", "y")], [("e", "y")] + rr))
    A3b.append((-1, [U("b", "c", "y"), U("b", "d", "x"), U("a", "e", "u")], [("c", "y"), ("d", "x"), ("e", "u")]))
    A3b.append((-1, [U("b", "c", "y"), U("a", "e", "u"), U("b", "d", "x")], [("c", "y"), ("e", "u"), ("d", "x")]))
    T += prod(-c, BR1, A3b, "4r_A3")
    A4b = []
    for pe in ("y", "z"):
        sgn = +1 if pe == "y" else -1
        fe = fl["fA4"] and pe == "z"
        A4b.append((sgn, [U("b", "c", "x"), U("a", "d", "u"), Uo("b", "e", pe, fe)], [("c", "x"), ("d", "u"), ("e", "y")]))
        A4b.append((sgn, [U("a", "d", "u"), U("b", "c", "x"), Uo("b", "e", pe, fe)], [("d", "u"), ("c", "x"), ("e", "y")]))
    T += prod(-c, BR1, A4b, "4r_A4")

    B1a = []
    for rr in anti("b", "x", "cp", "xp"):
        B1a.append((+1, [Uo("cp", "a", "wp", fl["fB1a_wp"])], rr))
    B1a.append((-1, [U("cp", "dp", "x"), U("a", "d", "xp"), Uo("b", "cp", "z", fl["fB1a_z"])], [("dp", "x"), ("d", "xp")]))
    B1a.append((-1, [U("a", "d", "xp"), U("cp", "dp", "x"), Uo("b", "cp", "z", fl["fB1a_z"])], [("d", "xp"), ("dp", "x")]))
    B1b = []
    for rr in anti("c", "u", "b", "y"):
        B1b.append((+1, [U("a", "c", "w")], rr))
    B1b.append((-1, [Uo("b", "c", "z", fl["fB1b_z"]), U("c", "e", "y"), U("a", "bp", "u")], [("e", "y"), ("bp", "u")]))
    B1b.append((-1, [Uo("b", "c", "z", fl["fB1b_z"]), U("a", "bp", "u"), U("c", "e", "y")], [("bp", "u"), ("e", "y")]))
    T += prod(+c * 0.5, B1a, B1b, "4r_B1")
    B2a = [(+1, [U("a", "bp", "xp"), U("cp", "ap", "x"), Uo("cp", "b", "z", fl["fB2a1"])], [("bp", "xp"), ("ap", "x")]),
           (-1, [Uo("b", "dp", "z", fl["fB2a2_z"]), U("dp", "ap", "x"), Uo("cp", "a", "wp", fl["fB2a2_wp"])], [("cp", "xp"), ("ap", "x")])]
    B2b = [(+1, [Uo("b", "c", "z", fl["fB2b1"]), U("c", "d", "y"), U("a", "e", "u")], [("d", "y"), ("e", "u")]),
           (-1, [U("a", "c", "w"), U("d", "e", "y"), Uo("d", "b", "z", fl["fB2b2"])], [("e", "y"), ("c", "u")])]
    T += prod(+c * 2, B2a, B2b, "4r_B2")
    B3a = [(+1, [U("cp", "dp", "x"), U("a", "ep", "xp"), Uo("b", "cp", "z", fl["fB3a"])], [("dp", "x"), ("ep", "xp")]),
           (-1, [U("a", "bp", "xp")], [("b", "x"), ("bp", "xp")])]
    B3b = [(+1, [Uo("b", "c", "z", fl["fB3b"]), U("a", "d", "u"), U("c", "e", "y")], [("d", "u"), ("e", "y")]),
           (-1, [U("a", "c", "u")], [("c", "u"), ("b", "y")])]
    T += prod(+c * 2, B3a, B3b, "4r_B3")

    C1a = [(+1, [U("a", "bp", "xp"), U("cp", "ap", "x"), U("cp", "b", "z")], [("bp", "xp"), ("ap", "x")]),
           (-1, [U("dp", "b", "z"), U("dp", "ap", "x"), U("a", "cp", "wp")], [("cp", "xp"), ("ap", "x")])]
    C1b = []
    for rr in anti("c", "u", "b", "y"):
        C1b.append((+1, [U("a", "c", "w")], rr))
    C1b.append((-1, [U("c", "b", "z"), U("c", "e", "y"), U("a", "ep", "u")], [("e", "y"), ("ep", "u")]))
    C1b.append((-1, [U("c", "b", "z"), U("a", "ep", "u"), U("c", "e", "y")], [("ep", "u"), ("e", "y")]))
    T += prod(+c, C1a, C1b, "4r_C1")
    C2a = [(+1, [U("cp", "dp", "x"), U("a", "ep", "xp"), U("cp", "b", "z")], [("dp", "x"), ("ep", "xp")]),
           (-1, [U("a", "bp", "xp")], [("b", "x"), ("bp", "xp")])]
    C2b = []
    for rr in anti("c", "u", "b", "y"):
        C2b.append((+1, [U("a", "c", "w")], rr))
    C2b.append((-1, [U("c", "b", "z"), U("c", "e", "y"), U("a", "ap", "u")], [("e", "y"), ("ap", "u")]))
    C2b.append((-1, [U("c", "b", "z"), U("a", "ap", "u"), U("c", "e", "y")], [("ap", "u"), ("e", "y")]))
    T += prod(+c, C2a, C2b, "4r_C2")
    C3a = C2a
    C3b = [(+1, [U("c", "b", "z"), U("c", "d", "y"), U("a", "e", "u")], [("d", "y"), ("e", "u")]),
           (-1, [U("a", "c", "w"), U("d", "e", "y"), U("d", "b", "z")], [("e", "y"), ("c", "u")])]
    T += prod(+c * 2, C3a, C3b, "4r_C3")

    out = []
    for t in T:
        hc = t.tag.startswith("4r_A") or t.tag.startswith("4r_C")
        out.append((t, hc))
    return out

# ---------- per-term cached ledger contributions ----------
_cache = {}
def sig(t):
    return (t.coef, tuple(t.kern), tuple(t.color), tuple(t.rhos), t.phase, t.Nc_pow)

def contribs(t):
    s = sig(t)
    if s in _cache:
        return _cache[s]
    out = []
    for ct in canonical_order(t):
        sec = len(ct.rhos)
        if sec not in (2, 3, 4):
            continue
        n, dA = normalize_term2(ct)
        if n is None:
            continue
        key, sgn = canon_full(n)
        if key[0] in ("ZERO", "SELFZERO"):
            continue
        out.append(((sec, key, n.Nc_pow, dA), n.coef * sgn))
    _cache[s] = out
    return out

def combo_ledger(fl, with_hc=True):
    led = {}
    for t, hc in build(fl):
        variants = [t] + ([sym_hc(t)] if (with_hc and hc) else [])
        for v in variants:
            for k, cv in contribs(v):
                led[k] = led.get(k, 0.0) + cv
    return led

def sector(led, sec):
    return {k: v for k, v in led.items() if k[0] == sec and abs(v) > 1e-12}

def led_of_terms(terms):
    led = {}
    for t in terms:
        for k, cv in contribs(t):
            led[k] = led.get(k, 0.0) + cv
    return led

def diff_norm(ledA, ledB):
    keys = set(ledA) | set(ledB)
    return sum(abs(ledA.get(k, 0.0) - ledB.get(k, 0.0)) for k in keys), len(
        [k for k in keys if abs(ledA.get(k, 0.0) - ledB.get(k, 0.0)) > 1e-12])

if __name__ == "__main__":
    # targets: collection's from-4rho 3rho blocks
    tgt_terms = (parse_block(block_body(F3, 182, 187), extra_coef=-1j / (16 * PI**5), tag="f4a")
                 + parse_block(block_body(F3, 188, 193), extra_coef=+1j / (16 * PI**5), tag="f4b"))
    T1 = sector(led_of_terms(tgt_terms), 3)
    T2 = sector(led_of_terms(tgt_terms + [sym_hc(t) for t in tgt_terms]), 3)
    print(f"target T1 (as written): {len(T1)} structures; T2 (+hc): {len(T2)}")

    best = []
    # warm the cache for all variants once (loop over single-flip builds)
    import time
    t0 = time.time()
    results = []
    for bits in itertools.product((False, True), repeat=len(FLIPS)):
        fl = dict(zip(FLIPS, bits))
        led = combo_ledger(fl, with_hc=True)
        s4 = sector(led, 4)
        if s4:
            continue
        s3 = sector(led, 3)
        for tname, tgt in (("T1", T1), ("T2", T2)):
            dn, nk = diff_norm(s3, tgt)
            results.append((nk, dn, tname, fl))
            if nk == 0:
                print("WINNER:", tname, {k: v for k, v in fl.items() if v})
    results.sort(key=lambda r: (r[0], r[1]))
    print(f"\nscan time {time.time()-t0:.1f}s; combos with 4rho sector cancelling: "
          f"{len(results)//2 if results else 0}")
    for nk, dn, tname, fl in results[:10]:
        print(f"  resid structs={nk:3d} norm={dn:.3e} vs {tname}  flips={ {k for k,v in fl.items() if v} }")
