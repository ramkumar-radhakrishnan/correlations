"""Check 1: 3-rho sector of the user's collected set (Evolution_three_rho, lines 180-255)."""
import math
import numpy as np
from engine import Term, Config, accumulate, compare, total_norm
from texparse import parse_block
from blocks import block_body, F3, PI

BLOCKS = [
    ("from4rho_a", 182, 187, -1j / (16 * PI**5)),
    ("from4rho_b", 188, 193, +1j / (16 * PI**5)),
    ("GI_R2_N1",   200, 203, +1j / (16 * PI**5)),
    ("GI_R2_N2",   207, 209, +1j / (16 * PI**5)),
    ("GI_R3",      211, 214, +1j / (16 * PI**5)),
    ("GI_R4",      216, 222, +1j / (16 * PI**5)),
    ("GI_R5_1",    224, 226, +1j / (16 * PI**5)),
    ("GI_R5_2",    227, 231, +1j / (16 * PI**5)),
    ("GII_R1",     233, 235, +1j / (16 * PI**5)),
    ("GIII_R1",    237, 239, +1j / (16 * PI**5)),
    ("GIII_R2",    241, 243, +1j / (16 * PI**5)),
    ("GIII_R3",    245, 247, +1j / (16 * PI**5)),
    ("GIII_R5",    249, 251, +1j / (16 * PI**5)),
    ("GIII_R6",    253, 255, +1j / (16 * PI**5)),
]

def get_block_terms():
    per = {}
    for tag, l0, l1, coef in BLOCKS:
        body = block_body(F3, l0, l1)
        per[tag] = parse_block(body, extra_coef=coef, tag=tag)
    return per

if __name__ == "__main__":
    per = get_block_terms()
    print("terms per block:", {k: len(v) for k, v in per.items()})
    allterms = [t for ts in per.values() for t in ts]
    cfg = Config(6, "su2", seed=11)
    # per-block scale
    print("--- per-block 3rho norms (scale reference) ---")
    for tag, ts in per.items():
        n = total_norm(accumulate(ts, cfg, sector_filter={3}))
        print(f"  {tag}: {n.get(3, 0):.4e}")
    print("--- total, as written ---")
    tot = accumulate(allterms, cfg)
    print("  norms:", {k: f"{v:.4e}" for k, v in total_norm(tot).items()})
    print("--- total + h.c. ---")
    tot2 = accumulate(allterms + [t.hc() for t in allterms], cfg)
    print("  norms:", {k: f"{v:.4e}" for k, v in total_norm(tot2).items()})
    print("--- anti-hermitian part: total - h.c. ---")
    tot3 = accumulate(allterms, cfg)
    tot3h = accumulate([t.hc() for t in allterms], cfg)
    for sec in sorted(set(tot3) | set(tot3h)):
        A = tot3.get(sec); B = tot3h.get(sec)
        if A is None: A = np.zeros_like(B)
        if B is None: B = np.zeros_like(A)
        print(f"  sector {sec}: |T|={np.abs(A).max():.4e}  |T+hc|={np.abs(A+B).max():.4e}  |T-hc|={np.abs(A-B).max():.4e}")
