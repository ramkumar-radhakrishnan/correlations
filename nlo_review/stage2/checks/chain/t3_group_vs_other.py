# Task 3 (two-rho sources): Group files (Lambda, vee kept) vs Other_terms (Lambda->0, vee kept).
# Model tested:  Group - Other = a*log(k+/Lambda) + b*log(vee/k+) + c   with a,b independent of Lambda,vee; report c (should be 0).
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from t1_special import addd, scale, evt
def fixG(name, t):
    if name == "GII RowI 2f": 
        assert t.count(r"U^{ac}(w)\,^{c'b}(z)") == 1 or t.count(r"U^{ac}(w)^{c'b}(z)") == 1, "typo pattern"
        t = t.replace(r"U^{ac}(w)\,^{c'b}(z)", r"U^{ac}(w)\,U^{c'b}(z)").replace(r"U^{ac}(w)^{c'b}(z)", r"U^{ac}(w)U^{c'b}(z)")
    return t
PAIRS = [("GII RowI 2a", "G2", 65, 65, 159, 160), ("GII RowI 2b", "G2", 68, 68, 163, 165), ("GII RowI 2c", "G2", 71, 71, 168, 169),
         ("GII RowI 2d", "G2", 74, 74, 172, 173), ("GII RowI 2e", "G2", 77, 77, 176, 177), ("GII RowI 2f", "G2", 80, 80, 180, 182),
         ("GI RowV 2a", "G1", 244, 244, 130, 130), ("GI RowV 2b rest", "G1", 250, 262, 133, 142)]
def D(env, name, g, a, b, oa, ob, Lam, vee):
    e = dict(env); e["Lam"] = Lam * env["kp"]; e["vee"] = vee * env["kp"]
    G = evt(fixG(name, rhs(lines(g, a, b)).rstrip().rstrip(".")), e); O = evt(rhs(lines("O", oa, ob)), e)
    keys = sorted(set(G) | set(O), key=str)
    return keys, np.concatenate([np.ravel(np.asarray(G.get(k, 0) - O.get(k, 0), dtype=complex)) for k in keys]), max(np.max(np.abs(v)) for v in O.values())
for N, seed in ((2, 21), (3, 22)):
    env = make_env(N, seed); env["pp"] = env["qp"] = 0.37 * env["kp"]; print(f"--- N={N} seed={seed}")
    for row in PAIRS:
        try:
            L1, L2, V1, V2 = 1e-9, 1e-12, 5.0, 11.0
            k, d11, sc = D(env, *row, L1, V1); _, d21, _ = D(env, *row, L2, V1); _, d12, _ = D(env, *row, L1, V2); _, d22, _ = D(env, *row, L2, V2)
            a = (d21 - d11) / (np.log(1 / L2) - np.log(1 / L1)); b = (d12 - d11) / (np.log(V2) - np.log(V1))
            c = d11 - a * np.log(1 / L1) - b * np.log(V1); chk = d22 - a * np.log(1 / L2) - b * np.log(V2) - c
            print(f"{row[0]:16s} {row[1]}:{row[2]}-{row[3]} vs O{row[4]}-{row[5]}: |Other|={sc:.2e} |a|={np.max(np.abs(a)):.2e} |b|={np.max(np.abs(b)):.2e} remainder |c|={np.max(np.abs(c)):.1e} 4th-point residual={np.max(np.abs(chk)):.1e}")
        except Exception as ex: print(row[0], "ERROR", type(ex).__name__, str(ex)[:300])
