"""Per-row comparison: MY group-file L_J 3rho terms  vs  user's collected blocks."""
import numpy as np
from engine import Term, canonical_order
from symcolor import ledger_full2 as ledger_full, pretty_full2 as pretty_full
from run_test1 import sym_hc
import sources_LJ
from fourrho import terms_4rho
from g3cancel import get_block_terms

def neg(t):
    return Term(-t.coef, t.kern, t.color, t.rhos, phase=t.phase, Nc_pow=t.Nc_pow, tag="NEG_" + t.tag)

def mine_by_tag(prefixes, with_hc):
    base = list(sources_LJ.SRC) + list(terms_4rho())
    out = []
    for t, hc in base:
        if not any(t.tag.startswith(p) for p in prefixes):
            continue
        variants = [t]
        if with_hc and hc:
            variants.append(sym_hc(t))
        for v in variants:
            for ct in canonical_order(v):
                if len(ct.rhos) == 3:
                    out.append(ct)
    return out

PAIRS = [
    ("GI_R2",   ["GI_R2_2a"],                 ["GI_R2_N1", "GI_R2_N2"]),
    ("GI_R3",   ["GI_R3_K2", "GI_R3_R2"],     ["GI_R3"]),
    ("GI_R4",   ["GI_R4a", "GI_R4b"],         ["GI_R4"]),
    ("GI_R5p1", ["GI_R5p1_e1", "GI_R5p1_e2"], ["GI_R5_1", "GI_R5_2"]),
    ("GII_R1",  ["GII_R1p1_A", "GII_R1p1_B"], ["GII_R1"]),
    ("GIII_R1", ["GIII_R1"],                  ["GIII_R1"]),
    ("GIII_R2", ["GIII_R2"],                  ["GIII_R2"]),
    ("GIII_R3", ["GIII_R3p1"],                ["GIII_R3"]),
    ("GIII_R5", ["GIII_R5"],                  ["GIII_R5"]),
    ("GIII_R6", ["GIII_R6"],                  ["GIII_R6"]),
    ("4rho",    ["4r_"],                      ["from4rho_a", "from4rho_b"]),
]

if __name__ == "__main__":
    theirs = get_block_terms()
    for name, mytags, theirtags in PAIRS:
        tt = [t for tag in theirtags for t in theirs[tag]]
        for with_hc in (False, True):
            mm = mine_by_tag(mytags, with_hc)
            diff = mm + [neg(t) for t in tt]
            led, surv = ledger_full(diff, f"{name} (mine{'+hc' if with_hc else ''} - theirs)")
            if not surv:
                print(f"   --> {name}: EXACT MATCH (hc={'yes' if with_hc else 'no'})")
                break
        else:
            # report how big the residue is relative to each side
            _, s_m = ledger_full(mm, name + " mine-only", tol=1e-12)
            _, s_t = ledger_full([neg(t) for t in tt], name + " theirs-only", tol=1e-12)
            print(f"   !! {name}: mismatch; mine {len(s_m)} structs, theirs {len(s_t)} structs, diff {len(surv)}")
