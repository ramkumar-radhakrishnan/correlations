# Power counting of every (2rho) contact line at every coincidence that involves the merged charge point,
# with a smooth SU(3) Wilson-line field (colour factor included) and with the colour factor set to 1.
# Exponent alpha: |f| ~ lambda^alpha when the points of the cluster S approach each other at scale lambda.
# Absolutely integrable at S iff alpha > -d_S, d_S = 2(|S|-1) (minus 1 if the delta(s) constraint ties y' to z in S).
import sys, os, itertools, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lib, lines

LAMS = [1e-2, 1e-3, 1e-4, 1e-5]

def base_config(unit, rng):
    P = {}
    for p in unit.base:
        P[p] = rng.uniform(-1.5, 1.5, size=2)
    return P

def make_P(unit, base, S, lam, dirs, ubase, ucl):
    """positions with the cluster S scaled to size lam around the merged point (or the first point of S)."""
    P = {p: base[p].copy() for p in unit.points}
    c = base[unit.keep] if unit.keep in S else base[S[0]]
    for q in S:
        P[q] = c + lam * dirs[q]
    if unit.keep in S: P[unit.keep] = c.copy()
    else: P[S[0]] = c.copy()
    if unit.dist == 'delta':
        if "y'" in S and 'z' in S:
            P["y'"] = P['z'] + lam * ucl * lib.EPERP
        elif "y'" in S:                      # z on the line through y' perpendicular to k
            P['z'] = np.array([P["y'"][0], P["y'"][1] - ubase])
        else:
            P["y'"] = P['z'] + ubase * lib.EPERP
    for m in getattr(unit, 'extra_merge', []):
        P[m] = P[unit.keep].copy()
    return P

def evaluate(unit, col, P):
    Q = unit.full(P)
    for m in getattr(unit, 'extra_merge', []): Q[m] = Q[unit.keep]
    k = unit.kern(Q); c = unit.colour(col, Q)
    v1 = unit.coef * k * c; v0 = unit.coef * k
    if unit.kern2 is not None:
        k2 = unit.kern2(Q); c2 = unit.colour2(col, Q)
        v1 = v1 + unit.coef * k2 * c2; v0 = v0 + unit.coef * k2
    if unit.dist == 'pv':
        s = lib.s3(Q); v1 = v1 / s; v0 = v0 / s
    return v1, v0

def dim(unit, S):
    d = 2 * (len(S) - 1)
    if unit.dist == 'delta' and "y'" in S and 'z' in S: d -= 1
    return d

def scan_unit(unit, field, nconf=3, nor=16, seed=0, clusters=None):
    col = lib.Colour(field)
    rng = np.random.default_rng(seed)
    pts = [p for p in unit.points if p not in getattr(unit, 'extra_merge', [])]
    others = [p for p in pts if p != unit.keep]
    if unit.dist == 'delta': others = [p for p in others]  # y' allowed (constrained)
    if clusters is None:
        clusters = []
        for n in range(1, len(others) + 1):
            for comb in itertools.combinations(others, n):
                clusters.append([unit.keep] + list(comb))
    res = []
    bases = [base_config(unit, rng) for _ in range(nconf)]
    for S in clusters:
        mx1 = np.zeros(len(LAMS)); mx0 = np.zeros(len(LAMS)); mean2 = np.zeros(len(LAMS), complex)
        for base in bases:
            ubase = rng.choice([-1, 1]) * rng.uniform(0.5, 1.5)
            for o in range(nor):
                dirs = {}
                for q in S:
                    ang = rng.uniform(0, 2 * np.pi); r = rng.uniform(0.5, 1.5)
                    dirs[q] = r * np.array([np.cos(ang), np.sin(ang)])
                ucl = rng.choice([-1, 1]) * rng.uniform(0.5, 1.5)
                for i, lam in enumerate(LAMS):
                    P = make_P(unit, base, S, lam, dirs, ubase, ucl)
                    v1, v0 = evaluate(unit, col, P)
                    mx1[i] = max(mx1[i], abs(v1)); mx0[i] = max(mx0[i], abs(v0))
        a1 = np.log10(mx1[-1] / mx1[-2]) / np.log10(LAMS[-1] / LAMS[-2])
        a1b = np.log10(mx1[-2] / mx1[-3]) / np.log10(LAMS[-2] / LAMS[-3])
        a0 = np.log10(mx0[-1] / mx0[-2]) / np.log10(LAMS[-1] / LAMS[-2])
        res.append(dict(S=S, d=dim(unit, S), a_col=a1, a_col_prev=a1b, a_nocol=a0,
                        max_col=list(mx1), max_nocol=list(mx0)))
    return res

def circle_mean(unit, field, q, nconf=3, nang=96, seed=1):
    """angular mean of the line around the merged point for the pair cluster {keep, q} (exponent of the mean)."""
    col = lib.Colour(field); rng = np.random.default_rng(seed); out = []
    for _ in range(nconf):
        base = base_config(unit, rng); ubase = 0.9
        m = []
        for lam in LAMS:
            acc = 0
            for t in np.linspace(0, 2 * np.pi, nang, endpoint=False):
                dirs = {q: np.array([np.cos(t), np.sin(t)])}
                P = make_P(unit, base, [unit.keep, q], lam, dirs, ubase, 1.0)
                acc += evaluate(unit, col, P)[0]
            m.append(abs(acc / nang))
        out.append(np.log10(m[-1] / m[-2]) / np.log10(LAMS[-1] / LAMS[-2]))
    return out

if __name__ == '__main__':
    field = lib.SmoothField(seed=7)
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = []
    for iu, u in enumerate(lines.U):
        if which != 'all' and which != u.row: continue
        r = scan_unit(u, field, seed=iu)
        for e in r:
            flag = 'OK' if e['a_col'] > -e['d'] + 0.2 else ('MARGINAL' if e['a_col'] > -e['d'] - 0.2 else 'FAIL')
            flag0 = 'ok' if e['a_nocol'] > -e['d'] + 0.2 else ('marg' if e['a_nocol'] > -e['d'] - 0.2 else 'fail')
            print('%-38s S=%-28s d=%d  alpha(col)=%6.2f (prev %6.2f)  alpha(U=const-free kernel)=%6.2f  %s/%s' % (
                u.label(), ','.join(e['S']), e['d'], e['a_col'], e['a_col_prev'], e['a_nocol'], flag, flag0), flush=True)
            out.append(dict(unit=u.label(), dist=u.dist, **e))
    json.dump(out, open('scan_%s.json' % which.replace('.', '_'), 'w'), indent=1, default=float)


def all_clusters(unit):
    pts = [p for p in unit.points if p not in getattr(unit, 'extra_merge', [])]
    cl = []
    for n in range(2, len(pts) + 1):
        for comb in itertools.combinations(pts, n):
            S = list(comb)
            if unit.keep in S: S = [unit.keep] + [p for p in S if p != unit.keep]
            cl.append(S)
    return cl

def run_variant(tag, units, field, allc=False):
    out = []
    with open('scan_%s.out' % tag, 'w') as fo:
        for iu, u in enumerate(units):
            r = scan_unit(u, field, seed=iu, clusters=all_clusters(u) if allc else None)
            for e in r:
                flag = 'OK' if e['a_col'] > -e['d'] + 0.2 else ('MARGINAL' if e['a_col'] > -e['d'] - 0.2 else 'FAIL')
                flag0 = 'ok' if e['a_nocol'] > -e['d'] + 0.2 else ('marg' if e['a_nocol'] > -e['d'] - 0.2 else 'fail')
                fo.write('%-44s S=%-26s d=%d  alpha(col)=%6.2f (prev %6.2f)  alpha(kernel only)=%6.2f  %s/%s\n' % (
                    u.label() + (' [n local: %s->%s]' % (u.extra_merge[0], u.keep) if getattr(u, 'extra_merge', []) else ''),
                    ','.join(e['S']), e['d'], e['a_col'], e['a_col_prev'], e['a_nocol'], flag, flag0))
                out.append(dict(unit=u.label(), extra=getattr(u, 'extra_merge', []), dist=u.dist, **e))
    json.dump(out, open('scan_%s.json' % tag, 'w'), indent=1, default=float)
    return out
