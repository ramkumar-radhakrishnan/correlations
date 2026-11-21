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
# ---- longitudinal-momentum integrals with plus prescriptions (author's notation; replaces the closed-form functions) ----
def EP(a):  # recoil phase e^{-ik.(a) p+/k+} (red: differs from the author's (b-a))
    return r"e^{-ik\cdot\ph{(%s)}\frac{p^{+}}{k^{+}}}" % a
PL1 = r"\Big[\frac{1}{p^{+}}\Big]_{+}"
PLK = r"\Big[\frac{k^{+}}{p^{+}}\Big]_{+}"
PLKQ = r"\Big[\frac{k^{+}}{k^{+}-p^{+}}\Big]_{+}"
PL3 = r"\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}"
DPK = r"\int_{0}^{k^{+}}\frac{dp^{+}}{2\pi}"     # Group I: upper limit k+
DPI = r"\int_{0}^{\infty}\frac{dp^{+}}{2\pi}"    # Group III: upper limit vee-k+ -> infinity
def SE(a, two_pi=False):   # soft integral of Group I (= former E(k.(a))); plain dp+ unless two_pi (then dp+/(2pi), prefactor x 2pi)
    return (DPK if two_pi else r"\int_{0}^{k^{+}}dp^{+}")+PL1+EP(a)
def SP(a, two_pi=False):   # soft integral of Group III (= former Psi(k.(a)))
    return (DPI if two_pi else r"\int_{0}^{\infty}dp^{+}")+PL3+EP(a)
KZ_LE = r"\bigg[-\frac{1}{2}\log\frac{(z-w)^{2}}{(x-w)^{2}}-"+SE('z-w')+r"\bigg]\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)"
KZ_E = SE('z-w')+r"\,\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)"
KW = r"-"+SE('z-w')+r"\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)"
# Group I Row V first part, inside pref (1/k+^2) int_0^k+ dp+/(2pi) int e^{-ik(z-w')} [ ... ]
KV1_J = EP('x-z')+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(z-x)\,\mathcal{K}(y-z)\cdot\mathcal{K}(y'-x)"
KV1_E1 = r"-"+PLK+EP('x-z')+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)\,\mathcal{K}(z-x)\cdot\mathcal{K}(y'-x)"
KV1_E2 = r"-"+PLK+r"e^{-ik\cdot\ph{(x-z)}}"+EP('z-x')+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y'-x)\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)"
PF_I5 = r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{16\pi^{4}}\,\frac{1}{k^{+\,2}}"
def KT(wa, wb):  # Group II Row I / III.II kernel; wa: point of the amplitude-side gluon, wb of the conjugate side
    xp = "x'-%s" % wb
    return (r"\bigg[\frac{\mathcal{K}(%s)\cdot\mathcal{K}(y-z)}{2\,(z-%s)^{2}}+\mathcal{K}(%s)\cdot\mathcal{K}(z-%s)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)" % (xp, wa, xp, wa),
            r"-\mathcal{K}(%s)\cdot\mathcal{K}(x-%s)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)-\mathcal{K}(%s)\cdot\mathcal{K}(x-%s)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-%s)\bigg]\log\frac{(x-z)^{2}}{(x-%s)^{2}}" % (xp, wa, xp, wa, wa, wa))
# Group III rows, inside pref (x 2pi) int_0^inf dp+/(2pi) int e^{-ik(w-y')} e^{-ik(z-y') p+/k+} [ ... ]
PH_G3 = EP("z-y'")
G3_REG = r"\frac{1}{p^{+}+k^{+}}\Big[\mathcal{K}(y'-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(y-w)\cdot\mathcal{K}(x-z)-\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)\Big]"
G3_SOFT = r"-"+PL3+r"\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)"
G3_INST = r"-\frac{1}{k^{+}}\,\mathcal{K}(y'-z)\cdot\mathcal{K}(y-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-y')"
G3_INSTZ = r"+\frac{1}{k^{+}}\,\mathcal{K}(y'-z)\cdot\mathcal{K}(y-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)"

# ---- i-epsilon (out-state) prescription of Group I Row III: Delta N kernels (prefactor (1/(2pi)^3) eta g^4/(8 pi^4) 1/k+)
KDZ = [r"\bigg[\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)}{2\,(z-w)^{2}}+\mathcal{K}(x'-w')\cdot\mathcal{K}(z-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)",
       r"-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-w)-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\bigg]"]
KDW = [r"-\theta\big(\rho_{0}-|z-w|\big)\,\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}{2\,(z-w)^{2}}+\pi\,\delta^{(2)}(z-w)\,\log\frac{\rho_{0}}{|x-w|}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)"]
KDW_COUL = r"+\pi\,\delta^{(2)}(z-w)\,\log\frac{k^{+}}{\Lambda}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)"
PF_DN = r"\frac{1}{(2\pi)^{3}}\,\frac{\eta\,g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}"
