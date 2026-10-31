r"""The BFKL terms of line 1 of Group III Rows III(1), III(2), V, VI and Group II Row IV, by exact identities only
(vee > k+, p+ > Lambda; no large-p+ step).  For each row, line 1 is rewritten as

  G3:  line 1 = A + Bterm + S + Gz
         A     : line 1 with [1/p+]_+  ->  [k+/(p+(p+ + k+))]_+ + 1/(p+ + k+)  and without the BFKL-carrying term
         Bterm : the BFKL-carrying term (-1/k+ delta_im delta_jk', or -p+/k+^2 [...] in III(2)) with G(y') -> G(y') - G(z)
         S     : + f(0) int_0^{vee-k+} dp+/(p+ + k+) = f(0) log(vee/k+)  =  minus the log(vee/k+) part of line 2
         Gz    : the BFKL-carrying term with G(y') -> G(z)
       G(y') = (x'-y')^{k'}/(x'-y')^2 x (the colour bracket with U(y')).
       Gz after the y' integral (Fourier identity, exact):  GzFT = Gz with e^{-ik(y'-w)} e^{-ik p+/k+ (y'-z)} K(y'-z)
       -> e^{-ik(w'-w)} k+/(p+ + k+) K(w'-z);  GzFT = (B)  (III(2): (B) + Rem, partial fractions).
  G2 IV: line 1 = A' + term0,  A' = line 1 with xi/(1-xi)_+ -> xi + xi^2/(1-xi)_+ ;  term0 = your BFKL.tex G2 IV line 1.
       H4 is the (1-xi)/xi term with -(d_perp/2) log((y'-y)^2 mu^2/4) and y' -> y in kernel and colour, after the y' integral
       (int d^2r e^{-iq.r} log r^2 = -4 pi/q^2, int_{k+/vee}^1 dxi (1-xi)/xi = log(vee/k+) - 1 + k+/vee): coefficients in
       the note; the Fourier identities are checked here numerically.

Pointwise in the transverse points (p+, xi integrals done numerically with ../rows_check/pe.py), two configurations.
Output: exact_check.json"""
import json
import os
import sys
import numpy as np
import rows as R
import extract as X
import texeval as TE

HERE = R.HERE
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'rows_check')))
import pe  # noqa: E402

pe.PLUS += [(r'\left[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\right]_{+}', r'\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}'),
            (r'\frac{\xi^{2}}{(1-\xi)_{+}}', r'\frac{\xi^{2}}{(1-\xi)}')]
UB, UO, EX = X.extract()

KPLUS = r'\left[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\right]_{+}'
OLDPLUS = r'\left[\frac{1}{p^{+}}\right]_{+}'
NEWPLUS = r'\left(' + KPLUS + r'+\frac{1}{(p^{+}+k^{+})}\right)'


def rep1(t, a, b):
    assert t.count(a) == 1, (t.count(a), a[:80])
    return t.replace(a, b)


# ---------------------------------------------------------------- the pieces of the four Group III rows
FR_FULL = r"\frac{(y'-z)^{m}\,(x'-y')^{k'}\,(y-w)^{i}\,(x-z)^{j}}{(y'-z)^{2}\,(x'-y')^{2}\,(y-w)^{2}\,(x-z)^{2}}"
FR_NOG = r"\frac{(y'-z)^{m}\,(y-w)^{i}\,(x-z)^{j}}{(y'-z)^{2}\,(y-w)^{2}\,(x-z)^{2}}"
FR_Z = r"\frac{(y'-z)^{m}\,(x'-z)^{k'}\,(y-w)^{i}\,(x-z)^{j}}{(y'-z)^{2}\,(x'-z)^{2}\,(y-w)^{2}\,(x-z)^{2}}"
BR3 = r"\left[\frac{\delta_{k'm}\,\delta_{ij}}{(p^{+}+k^{+})}-\left[\frac{1}{p^{+}}\right]_{+}\,\delta_{jm}\,\delta_{ik'}-\frac{1}{k^{+}}\, \delta_{im}\, \delta_{jk'}\right]"
BR3_A = r"\left[\frac{\delta_{k'm}\,\delta_{ij}}{(p^{+}+k^{+})}-" + NEWPLUS + r"\,\delta_{jm}\,\delta_{ik'}\right]"
BR3_B = r"\left[-\frac{1}{k^{+}}\, \delta_{im}\, \delta_{jk'}\right]"
COL = {'G3 III(1)': r"\left[U^{e'c'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]",
       'G3 V': r"\left[U^{e'c'}(y')\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]",
       'G3 VI': r"\left[U^{e'c'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]",
       'G3 III(2)': "\\left[f^{c'd'a}\\,U^{e'c'}(y')\\,U^{d'b}(z)\\,\\rho^{e'}(x')-f^{c'd'a}\\, \n    U^{bd'}(z)\\,U^{c'e'}(x')\\,\\rho^{e'}(x')\\right]"}
PH = r"e^{-ik\cdot (y'-w)}\, e^{-ik\cdot \frac{p^{+}}{k^{+}}\,(y'-z)}"
PH2 = r"e^{-ik\cdot (y'-w)}\,e^{-ik\cdot \frac{p^{+}}{k^{+}}\left(y'-z\right)}"
PH_FT = r"e^{-ik\cdot (w'-w)}\,\frac{k^{+}}{(p^{+}+k^{+})}"
# Row III part 2
FR2_FULL = r"\frac{(y'-z)^{m}\,(x'-y')^{k'}\,}{(y'-z)^{2}\,(x'-y')^{2}\,}"
FR2_NOG = r"\frac{(y'-z)^{m}}{(y'-z)^{2}}"
FR2_Z = r"\frac{(y'-z)^{m}\,(x'-z)^{k'}}{(y'-z)^{2}\,(x'-z)^{2}}"
B2_TERM = (r"-\frac{p^{+}}{k^{+\,2}}\,\bigg[\frac{(x-z)^{k'}\,(z-w)^{m}}{(z-w)^{2}} -\,\frac{(x-z)^{k'}\,(x-w)^{m}}{2\,(x-w)^{2}}\bigg]")
DEN2 = r"\frac{1}{\left[p^{+}\,(x-z)^{2}+k^{+}\,(x-w)^{2}\right]}"


def gfac(col, a="y'"):
    return r"\frac{(x'-%s)^{k'}}{(x'-%s)^{2}}\,%s" % (a, a, col.replace("(y')", '(%s)' % a))


def g3_texts(name):
    """original line 1, A, Bterm, Gz, GzFT as LaTeX (one term each)"""
    T = UO[name].terms[0]
    col = COL[name]
    colz = col.replace("(y')", '(z)')
    gdiff = r"\Bigg[" + gfac(col) + '-' + gfac(col, 'z') + r"\Bigg]"
    if name != 'G3 III(2)':
        A = rep1(T, BR3, BR3_A)
        Bt = rep1(rep1(rep1(T, BR3, BR3_B), FR_FULL, FR_NOG), col, gdiff)
        Gz = rep1(rep1(rep1(T, BR3, BR3_B), FR_FULL, FR_Z), col, colz)
        ph = PH
    else:
        A = rep1(rep1(T, B2_TERM, ''), OLDPLUS, NEWPLUS)
        # Bterm: keep only the -p+/k+^2 term inside the big bracket
        i, j = T.index(r'\notag \\ & \Bigg[') + len(r'\notag \\ & \Bigg['), T.index(r'\Bigg]\times')
        Bo = T[:i] + ' ' + B2_TERM + ' ' + T[j:]
        Bt = rep1(rep1(Bo, FR2_FULL, FR2_NOG), col, gdiff)
        Gz = rep1(rep1(Bo, FR2_FULL, FR2_Z), col, colz)
        ph = PH2
    GzFT = rep1(Gz, ph, PH_FT).replace("(y'-z)", "(w'-z)").replace(r"\int_{w,x',y',z,y,x}", r"\int_{w,x',w',z,y,x}") \
        .replace(r"\int_{w,y',z,y,x,x'}", r"\int_{w,w',z,y,x,x'}").replace(r"\int_{w,w',y',z,y,x}", r"\int_{w,w',z,y,x}")
    assert "y'" not in GzFT, GzFT[:300]
    return T, A, Bt, Gz, GzFT


REM2 = r"\frac{(x-w)^{2}}{(x-z)^{2}\,(p^{+}+k^{+})\,\left[p^{+}\,(x-z)^{2}+k^{+}\,(x-w)^{2}\right]}"


def iii2_rem(GzFT):
    """the non-logarithmic part of GzFT in Row III(2): -p+/k+^2 x k+/(p++k+) x 1/[p+A+k+B] = -1/(k+ A (p++k+)) + B/(A (p++k+)(p+A+k+B))"""
    t = rep1(GzFT, PH_FT, r"e^{-ik\cdot (w'-w)}")
    t = rep1(t, DEN2, '')
    return rep1(t, r"-\frac{p^{+}}{k^{+\,2}}", REM2)


# ---------------------------------------------------------------- Group II Row IV
BR4 = r"\left[\xi\,(1-\xi)+\frac{(1-\xi)}{\xi}+\frac{\xi}{(1-\xi)_{+}}\right]"
BR4_A = r"\left[\xi\,(1-\xi)+\frac{(1-\xi)}{\xi}+\xi+\frac{\xi^{2}}{(1-\xi)_{+}}\right]"


def g2iv_texts():
    T = UO['G2 IV'].terms[0]
    return T, rep1(T, BR4, BR4_A), UB['G2 IV'].prefix + ' (' + UB['G2 IV'].terms[0] + ') ' + UB['G2 IV'].suffix


def hb_texts():
    import bterms as BT
    return BT.hidden_B_text()


# ---------------------------------------------------------------- display forms (remaining terms)
BR3_DISP = (r"\left[\frac{\delta_{k'm}\,\delta_{ij}}{(p^{+}+k^{+})}-\old{\left[\frac{1}{p^{+}}\right]_{+}}\new{" + NEWPLUS
            + r"}\,\delta_{jm}\,\delta_{ik'}\old{-\frac{1}{k^{+}}\, \delta_{im}\, \delta_{jk'}}\right]")
BR4_DISP = r"\left[\xi\,(1-\xi)+\frac{(1-\xi)}{\xi}+\old{\frac{\xi}{(1-\xi)_{+}}}\new{\xi+\frac{\xi^{2}}{(1-\xi)_{+}}}\right]"
POLE4 = (r"\left[-\frac{d_{\perp}}{\epsilon}-d_{\perp}\,\gamma_{E}-\frac{d_{\perp}}{2}\,\log{\left(\frac{(y'-y)^{2}\,"
         r"\mu^{2}_{\overline{MS}}}{4}\right)}-1\right]")
LOG4 = r"\left[-\frac{d_{\perp}}{2}\,\log{\left(\frac{(y'-y)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right)}\right]"
K4 = r"\frac{(x-y)\cdot (x'-y')}{(x-y)^{2}\,(x'-y')^{2}}"
K4Y = r"\frac{(x-y)\cdot (x'-y)}{(x-y)^{2}\,(x'-y)^{2}}"
U4 = r"\left[U^{e'c}(y') - U^{e'c}(x')\right]"
LOGV_H4 = r'\log{\left[\frac{\vee}{k^{+}}\right]}'


def negate(t):
    t = t.strip()
    if t.startswith('&'):
        t = t[1:].strip()
    if t.startswith('+'):
        return '-' + t[1:]
    if t.startswith('-'):
        return '+' + t[1:]
    return '-' + t


def g3_display(name):
    """remaining line 1 of a Group III row: your line 1 with the changes marked, and the new (blue) terms"""
    T, A, Bt, Gz, GzFT = g3_texts(name)
    col = COL[name]
    if name != 'G3 III(2)':
        disp = rep1(T, BR3, BR3_DISP)
    else:
        disp = rep1(rep1(T, B2_TERM, r'\old{' + B2_TERM + '}'), OLDPLUS, r'\old{' + OLDPLUS + r'}\new{' + NEWPLUS + '}')
    bt = rep1(Bt, gfac(col) + '-' + gfac(col, 'z'), gfac(col) + r'\notag \\ & \qquad-' + gfac(col, 'z'))
    out = [(disp, None), (r'+\,' + bt, 'blue')]
    if name == 'G3 III(2)':
        out.append((r'+\,' + iii2_rem(GzFT).replace(r"\int_{w,w',z,y,x}", r"\int_{w,x',w',z,x}"), 'blue'))
    return out


def hsub_text():
    """the (1-xi)/xi term of line 1 of Group II Row IV with -(d_perp/2) log((y'-y)^2 mu^2/4) and y' -> y in kernel and colour"""
    T = UO['G2 IV'].terms[0]
    t = rep1(T, BR4, r"\left[\frac{(1-\xi)}{\xi}\right]")
    t = rep1(t, POLE4, LOG4)
    t = rep1(t, K4, K4Y)
    return rep1(t, U4, U4.replace("(y')", '(y)'))


def g2iv_display():
    T = UO['G2 IV'].terms[0]
    h4 = UB['G2 IV'].terms[1]
    assert h4.count(LOGV_H4) == 1
    h4rem = negate(h4).replace(LOGV_H4, r'\left(1-\frac{k^{+}}{\vee}\right)')
    return [(rep1(T, BR4, BR4_DISP), None), (r'-\,' + hsub_text(), 'blue'), (h4rem, 'blue')]


# ---------------------------------------------------------------- Fourier identities (radial integrals)
def fourier_checks():
    """the two 2-D Fourier transforms, as radial integrals with a regulator that is removed (a -> 0)"""
    from scipy import integrate, special
    out = {}
    # int d^2r e^{-iq.r} r^m/r^2 = -2 pi i q^m/q^2   <=>  q int_0^inf J_1(q r) dr = 1;  regulator e^{-a r}
    for q in (0.7, 1.9):
        vals = [integrate.quad(lambda r: np.exp(-a * r) * special.j1(q * r), 0, 200 / a, limit=8000)[0] * q
                for a in (1e-2, 3e-3, 1e-3)]
        out['int J1(qr) dr * q, q=%.1f (a -> 0)' % q] = [float(x) for x in vals]
    # int d^2r e^{-iq.r} log r^2 = -4 pi/q^2 (q != 0);  regulator e^{-a r^2}: value = 1 + O(a/q^2), extrapolated a -> 0
    for q in (0.7, 1.9):
        A = (4e-3, 2e-3, 1e-3)
        vals = []
        for a in A:
            f = lambda r: 2 * np.pi * r * special.j0(q * r) * np.log(r * r) * np.exp(-a * r * r)
            v = integrate.quad(f, 0, np.sqrt(60 / a), limit=8000)[0]
            vals.append(float(v / (-4 * np.pi / q**2)))
        ext = (A[1] * vals[2] - A[2] * vals[1]) / (A[1] - A[2])
        out['FT log r^2 / (-4 pi/q^2), q=%.1f (a -> 0)' % q] = vals + [float(ext)]
    return out


def run():
    out = {'fourier': fourier_checks()}
    HB = hb_texts()
    for seed in (1, 2):
        c = TE.cfg(seed)
        ctx = pe.Ctx(c, KP=1.0, V=20.0, L=1e-3, EPS=0.37)
        res = {}
        for name in ('G3 III(1)', 'G3 III(2)', 'G3 V', 'G3 VI'):
            T, A, Bt, Gz, GzFT = g3_texts(name)
            vT, vA, vB, vG = (pe.side_value(t, ctx) for t in (T, A, Bt, Gz))
            line2_log = pe.side_value(UO[name].prefix + ' (' + EX[name][1]['B'] + ') ' + UO[name].suffix, ctx)
            S = -line2_log
            vFT = pe.side_value(GzFT, ctx)
            vHB = pe.side_value(HB[name], ctx)
            rem = pe.side_value(iii2_rem(GzFT), ctx) if name == 'G3 III(2)' else 0.0
            res[name] = dict(
                identity=float(abs(vT - (vA + vB + S + vG)) / max(abs(vT), abs(vA), abs(vB), abs(vG))),
                GzFT_vs_B=float(abs(vFT - (vHB + rem)) / max(abs(vFT), abs(vHB))),
                sizes=[float(abs(x)) for x in (vT, vA, vB, S, vG, vFT, vHB, rem)])
        T, Ap, t0 = g2iv_texts()
        vT, vA, v0 = pe.side_value(T, ctx), pe.side_value(Ap, ctx), pe.side_value(t0, ctx)
        line2_log = pe.side_value(UO['G2 IV'].prefix + ' (' + EX['G2 IV'][1]['B'] + ') ' + UO['G2 IV'].suffix, ctx)
        res['G2 IV'] = dict(identity=float(abs(vT - (vA + v0)) / max(abs(vT), abs(vA), abs(v0))),
                            term0_plus_line2=float(abs(v0 + line2_log) / abs(v0)),
                            sizes=[float(abs(x)) for x in (vT, vA, v0)])
        out[str(seed)] = res
        for k, v in res.items():
            print(seed, k, {a: (b if a == 'sizes' else '%.1e' % b) for a, b in v.items() if a != 'sizes'})
    print(json.dumps(out['fourier'], indent=1))
    json.dump(out, open(os.path.join(HERE, 'exact_check.json'), 'w'), indent=1)


if __name__ == '__main__':
    if sys.argv[1:] == ['--fourier']:
        d = json.load(open(os.path.join(HERE, 'exact_check.json')))
        d['fourier'] = fourier_checks()
        print(json.dumps(d['fourier'], indent=1))
        json.dump(d, open(os.path.join(HERE, 'exact_check.json'), 'w'), indent=1)
    else:
        run()
