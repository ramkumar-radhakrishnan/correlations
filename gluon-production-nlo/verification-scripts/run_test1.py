"""TEST 1: Do MY sources (group files + 4rho) cancel in the 3-rho sector?
And do they agree with the user's collection ledger?"""
import numpy as np
from engine import Term, canonical_order, substitute_var, Config, accumulate, total_norm
from symcolor import ledger_report, pretty
import sources_LJ
from fourrho import terms_4rho

def sym_hc(t):
    if any(k[0] == 'fun' for k in t.kern):
        raise RuntimeError("fun kernel hc: " + t.tag)
    t2 = Term(t.coef.conjugate(), t.kern, t.color, list(reversed(t.rhos)),
              phase=(t.phase[1], t.phase[0]), Nc_pow=t.Nc_pow, tag=t.tag + "_hc")
    t2 = substitute_var(t2, "w", "TMPW")
    t2 = substitute_var(t2, "wp", "w")
    t2 = substitute_var(t2, "TMPW", "wp")
    return t2

def my_sources(sector3_only=True):
    """all L_J source terms with hc applied where flagged; canonical_order applied;
    returns list of terms (optionally only those whose final rho count == 3)."""
    base = list(sources_LJ.SRC) + list(terms_4rho())
    out = []
    for t, hc in base:
        variants = [t]
        if hc:
            if any(k[0] == 'fun' for k in t.kern):
                variants.append(None)  # handled numerically later; skip here
            else:
                variants.append(sym_hc(t))
        for v in variants:
            if v is None:
                continue
            for ct in canonical_order(v):
                if sector3_only and len(ct.rhos) != 3:
                    continue
                out.append(ct)
    return out

if __name__ == "__main__":
    terms3 = my_sources(sector3_only=True)
    print(f"my 3rho-sector terms (after reordering): {len(terms3)}")
    led, surv = ledger_report(terms3, "MY sources, 3rho sector")
    for k, v in sorted(surv.items(), key=lambda kv: -abs(kv[1][0]))[:40]:
        print("  ", pretty(k, v))
