# Shared TeX snippets (paper phase convention e^{-ik(unprimed-primed)}); phases wrapped in \ph{} (red).
AUTHOR = '/home/user/correlations/nlo_review/final_explicit/Final_terms_author_2026-11-20.tex'
LINES = open(AUTHOR).read().split('\n')
def L(a, b):  # author's lines a..b (1-indexed, inclusive)
    return '\n'.join(LINES[a-1:b])
REPL = [  # author's phase -> new phase (red); every replacement is counted by the caller
    (r"e^{-ik\cdot(w'-w)}", r"e^{-ik\cdot\ph{(w-w')}}"),
    (r"e^{-ik(w'-w)}", r"e^{-ik\cdot\ph{(w-w')}}"),
    (r"e^{-ik\cdot(w^{\prime}-w)}", r"e^{-ik\cdot\ph{(w-w^{\prime})}}"),
    (r"e^{-ik\cdot(y'-y)\frac{1}{\xi}}", r"e^{-ik\cdot\ph{(y-y')}\frac{1}{\xi}}"),
    (r"e^{-ik\cdot(y'-y)\frac{p^{+}}{k^{+}}}", r"e^{-ik\cdot\ph{(y-y')}\frac{p^{+}}{k^{+}}}"),
    (r"e^{-ik\cdot(y'-y)}", r"e^{-ik\cdot\ph{(y-y')}}"),
]
def phase(s, expect=None):
    n = 0
    for a, b in REPL:
        n += s.count(a); s = s.replace(a, b)
    if expect is not None: assert n == expect, (n, expect, s[:200])
    return s, n
K = r'\mathcal{K}'
LA = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')"
PH_WW = r"e^{-ik\cdot\ph{(w-w')}}"
PH_ZW = r"e^{-ik\cdot\ph{(z-w')}}"
PH_WY = r"e^{-ik\cdot\ph{(w-y')}}"
PF = r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}"
# ---------------- kernels (row = PF * A * int kernel * colour) ----------------
R_I2b = r"\int d^{2}p\,d^{2}q\;e^{iq\cdot(w-y)}\,e^{ip\cdot(y-x)}\,\frac{(x'-w')^{i}}{(x'-w')^{2}\,\mathbf{q}^{2}}\,\frac{1}{2}\bigg[\frac{R^{i}(\mathbf{q},\mathbf{p})}{\mathbf{p}^{2}}-\frac{R^{i}(\mathbf{q},\mathbf{q}-\mathbf{p})}{(\mathbf{q}-\mathbf{p})^{2}}\bigg]"
# Group I Row III, structure with U(z): 24-term combined kernel (lines) + (l - E) part + E part (pair D)
KZ_24 = [
 r"\bigg[\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)}{2\,(z-w)^{2}}+\mathcal{K}(x'-w')\cdot\mathcal{K}(z-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)",
 r"-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-w)-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\bigg]\log\frac{(x-z)^{2}}{(x-w)^{2}}",
 r"+\frac{\arccos\left(\frac{(x-w)^{2}-(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(z-w)^{2}}}\right)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}\,(z-w)^{2}}\,\bigg[\mathcal{K}(x'-w')\cdot(z-w)\,\mathcal{K}(y-z)\cdot(x-w)",
 r"-\frac{(x-w)\cdot(x-z)}{(x-z)^{2}}\,\mathcal{K}(x'-w')\cdot(z-w)\,\mathcal{K}(y-z)\cdot(x-z)",
 r"+\mathcal{K}(x'-w')\cdot(x-z)\,\mathcal{K}(y-z)\cdot(z-w)-\frac{(x-w)\cdot(x-z)}{(x-w)^{2}}\,\mathcal{K}(x'-w')\cdot(x-w)\,\mathcal{K}(y-z)\cdot(z-w)\bigg]",
 r"-\frac{1}{2\,(z-w)^{2}}\log\frac{(z-w)^{2}}{(x-w)^{2}}\,\mathcal{K}(x'-w')\cdot(z-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)",
]
KZ_LE = r"\bigg[-\frac{1}{2}\log\frac{(z-w)^{2}}{(x-w)^{2}}-\mathcal{E}\big(k\cdot\ph{(z-w)}\big)\bigg]\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)"
KZ_E = r"\mathcal{E}\big(k\cdot\ph{(z-w)}\big)\,\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)"
KW = r"-\mathcal{E}\big(k\cdot\ph{(z-w)}\big)\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)"
KV1_J = r"\mathcal{J}(s)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(z-x)\,\mathcal{K}(y-z)\cdot\mathcal{K}(y'-x)"
KV1_E = r"-\mathcal{E}(s)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)\,\mathcal{K}(z-x)\cdot\mathcal{K}(y'-x)-e^{-is}\,\mathcal{E}(-s)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y'-x)\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)"
def KT(wa, wb):  # Group II Row I / III.II kernel; wa: point of the amplitude-side gluon, wb of the conjugate side
    xp = "x'-%s" % wb
    return (r"\bigg[\frac{\mathcal{K}(%s)\cdot\mathcal{K}(y-z)}{2\,(z-%s)^{2}}+\mathcal{K}(%s)\cdot\mathcal{K}(z-%s)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)" % (xp, wa, xp, wa),
            r"-\mathcal{K}(%s)\cdot\mathcal{K}(x-%s)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)-\mathcal{K}(%s)\cdot\mathcal{K}(x-%s)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-%s)\bigg]\log\frac{(x-z)^{2}}{(x-%s)^{2}}" % (xp, wa, xp, wa, wa, wa))
DIST = r"\Big[\pi\delta(s)-i\,\mathcal{P}\frac{1}{s}\Big]"
KYP_EH = r"-\widehat{\mathcal{E}}(s)\Big[\mathcal{K}(y'-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(y-w)\cdot\mathcal{K}(x-z)-\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)\Big]"
KYP_PSI = r"+\Psi(s)\,\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)"
KYP_D = r"+" + DIST + r"\mathcal{K}(y'-z)\cdot\mathcal{K}(y-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-y')"
KZ3_D = r"-" + DIST + r"\mathcal{K}(y'-z)\cdot\mathcal{K}(y-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)"
S_G3 = r"s\equiv k\cdot\ph{(z-y')}"
S_G1V = r"s\equiv k\cdot\ph{(x-z)}"

# ---- i-epsilon (out-state) prescription of Group I Row III: Delta N kernels (prefactor (1/(2pi)^3) eta g^4/(8 pi^4) 1/k+)
KDZ = [r"\bigg[\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)}{2\,(z-w)^{2}}+\mathcal{K}(x'-w')\cdot\mathcal{K}(z-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)",
       r"-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-w)-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\bigg]"]
KDW = [r"-\theta\big(\rho_{0}-|z-w|\big)\,\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}{2\,(z-w)^{2}}+\pi\,\delta^{(2)}(z-w)\,\log\frac{\rho_{0}}{|x-w|}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)"]
KDW_COUL = r"+\pi\,\delta^{(2)}(z-w)\,\log\frac{k^{+}}{\Lambda}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)"
PF_DN = r"\frac{1}{(2\pi)^{3}}\,\frac{\eta\,g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}"
