# Task 1: rows of Final_terms (F) vs Limiting_Cases (L) that need a relabelling / a sum of sub-equations / a p+ integral
import sys; sys.path.insert(0, __file__.rsplit("/", 1)[0])
from ltx import *
from scipy.integrate import quad_vec
def addd(*ds):
    out = {}
    for d in ds:
        for k, v in d.items(): out[k] = out.get(k, 0) + v
    return out
def scale(d, c): return {k: c * v for k, v in d.items()}
def evt(t, env, tagmap=None): return contract(parse(t, env), tagmap)
def ev(f, a, b, env, tagmap=None): return evt(rhs(lines(f, a, b)), env, tagmap)

def row_GI_V_2a(env):
    # L147 = [int_0^k dp+ ...] + [int_{k+}^{0} dq+ ...]; second piece: q+ = k+ - p+, reversed limits give a factor -1
    t = rhs(lines("L", 147, 147)); mark = r"+ \frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}\, \int^{0}_{k^{+}}"
    assert t.count(mark) == 1; t1, t2 = t.split(mark); t2 = mark[1:] + t2
    A = ev("F", 256, 264, env)
    e2 = dict(env); e2["qp"] = env["kp"] - env["pp"]
    B = addd(evt(t1, env), scale(evt(t2, e2, {"plus[q@0]": "plus[p@k]"}), -1.0))
    return report("2rho GI  RowV 2a (F256-264 vs L147; L 2nd piece q+=k+-p+, int^0_k -> -int_0^k)", A, B)
F2 = [(281, 282), (285, 287), (290, 294), (297, 303), (306, 307)]
L2 = [(173, 173), (176, 177), (180, 180), (183, 184), (187, 188), (191, 193)]
def row_GII_I_second(env, verbose=True):
    A = addd(*[ev("F", a, b, env) for a, b in F2]); B = addd(*[ev("L", a, b, env) for a, b in L2])
    return report("2rho GII RowI 2nd part: sum F 2a..2e (281-307) vs sum L 2a..2f (173-193), d_perp=2", A, B)
def row_GII_IV(env):
    # L201-202 is written with xi in (0,1); F315 with p+: xi = k+/(p+ + k+), d xi = k+/(p+ + k+)^2 dp+
    A = ev("F", 315, 315, env)
    e2 = dict(env); e2["sym"] = dict(env["sym"]); k, p = env["kp"], env["pp"]; e2["sym"]["xi"] = k / (p + k)
    B = scale(ev("L", 201, 202, e2), k / (p + k) ** 2)
    return report("2rho GII RowIV (F315 vs L201-202; xi=k+/(p+ + k+), Jacobian k+/(p+ + k+)^2)", A, B)
def row_GIII_III_second(env):
    tF = rhs(lines("F", 324, 325)); mF = r"+\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\int_{w,x',w',z,x}e^{-ik\cdot(w'-w)}"
    assert tF.count(mF) == 1; F1, F2_ = tF.split(mF); F2_ = mF[1:] + F2_
    tL = rhs(lines("L", 222, 227)); L1 = rhs(lines("L", 222, 225)); L2_ = lines("L", 226, 227).strip()
    assert L2_.startswith("&  +"); L2_ = L2_[4:]
    r1 = report("2rho GIII RowIII 2nd part, p+-dependent terms (F324-325[first three] vs L222-225)", ev_t(F1, env), ev_t(L1, env))
    # last term: F has the p+ integral done (vee->inf), L has it as int_0^{vee-k+} dp+/(2pi) of a phase-free integrand
    A = ev_t(F2_, env); key = list(A)[0]
    def integrand(p):
        e2 = dict(env); e2["pp"] = p; d = ev_t(L2_, e2); assert list(d) == [key]; return d[key]
    val, err = quad_vec(integrand, 0, np.inf, epsabs=1e-16, epsrel=1e-11)
    r2 = report(f"2rho GIII RowIII 2nd part, last term (F325[4th], closed form  vs  numerical int_0^inf dp+ of L226-227; quad err {np.max(np.abs(err)):.1e})", A, {key: val})
    return max(r1, r2)
ev_t = evt
def row_GI_II_text():
    import difflib
    from tdiff_lc import norm
    a = norm(lines("F", 160, 168)); b = norm(lines("L", 44, 52))
    print("3rho GI  RowII N_2b (momentum-space formula, text comparison F160-168 vs L44-52):", "IDENTICAL" if a == b else "DIFFERENT")
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if tag != "equal": print(f"    [{tag}] F: ...{a[max(0,i1-25):i1]}[[{a[i1:i2]}]]{a[i2:i2+25]}...   L: [[{b[j1:j2]}]]")
if __name__ == "__main__":
    for N, seed in ((2, 1), (3, 2), (2, 7)):
        print(f"--- N={N} seed={seed}")
        for frac in (None, 0.07, 0.83):
            env = make_env(N, seed)
            if frac: env["pp"] = env["qp"] = frac * env["kp"]
            print(f"  p+/k+ = {env['pp']/env['kp']:.3f}")
            row_GI_V_2a(env)
            if frac is None: row_GII_I_second(env)
            for f2 in ([None] if frac is None else [frac, 2.7]):
                e = dict(env)
                if f2: e["pp"] = f2 * env["kp"]
                if f2 == 2.7: print("  p+/k+ = 2.7")
                row_GII_IV(e)
                if frac is None: row_GIII_III_second(e)
                else: report("2rho GIII RowIII 2nd part, p+-dependent terms", ev_t(rhs(lines("F", 324, 325)).split(r"+\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}")[0], e), ev("L", 222, 225, e))
    row_GI_II_text()
