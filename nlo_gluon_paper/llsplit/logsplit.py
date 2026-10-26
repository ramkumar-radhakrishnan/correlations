r"""How every logarithm of your final row expressions splits:

    log(...) = cL * log(k+/Lambda) + cV * log(vee/k+) + R,

R finite for Lambda -> 0 and vee -> infinity.  R is written with Lambda -> 0 where that is possible (your instruction);
for the logs whose remainder depends on Lambda only through log(1 - Lambda/k+) or log(1 + Lambda/k+), R = 0.

SPLIT maps the log as written (whitespace and \, removed) to (cL, cV, R, R_with_Lambda).  R_with_Lambda is the exact
remainder (used for the numerical identity check, R is its Lambda -> 0 limit).
HIDDEN: in Group III Rows III, V, VI and Group II Row IV the log(vee/k+) inside log(vee/k+ - 1) cancels against the one
hidden in the [1/p+]_+ (resp. (1-xi)_+) integral of the same row; there log(vee/k+ - 1) stays in the finite terms."""
import re

LKL = r"\log{\left(\frac{k^{+}}{\Lambda}\right)}"
LVK = r"\log{\left(\frac{\vee}{k^{+}}\right)}"
L1MA = r"\log{\left(1-\frac{k^{+}}{\vee}\right)}"
LVK1 = r"\log{\left(\frac{\vee}{k^{+}}-1\right)}"


def norm(s):
    return re.sub(r'\s+|\\,', '', s)


def lam_rem(P, Q):
    """log[((vee-k+) P + k+ Q)/(Lambda P + k+ Q)] = log(vee/k+) + log[((1-k+/vee) P + (k+/vee) Q)/((Lambda/k+) P + Q)]"""
    exact = (r"\log{\left[\frac{\left(1-\frac{k^{+}}{\vee}\right)%s+\frac{k^{+}}{\vee}\,%s}{\frac{\Lambda}{k^{+}}\,%s+%s}\right]}"
             % (P, Q, P, Q))
    lim = r"\log{\left[\frac{\left(1-\frac{k^{+}}{\vee}\right)%s+\frac{k^{+}}{\vee}\,%s}{%s}\right]}" % (P, Q, Q)
    return exact, lim


A, B = r"(x-z)^{2}", r"(x-w)^{2}"
Ap, Bp, Bw = r"(x'-z)^{2}", r"(x'-w')^{2}", r"(x-w')^{2}"

_raw = []


def add(s, cL, cV, R, Rx=None):
    _raw.append((s, cL, cV, R, R if Rx is None else Rx))


# ---- explicit large logarithms
add(LKL, 1, 0, None)
add(r"\log{\left[\frac{k^{+}}{\Lambda}\right]}", 1, 0, None)
add(r"\log\left(\frac{k^{+}}{\Lambda}\right)", 1, 0, None)
add(r"\log{\left[\frac{\vee}{k^{+}}\right]}", 0, 1, None)
add(r"\log\left(\frac{\vee}{k^{+}}\right)", 0, 1, None)
add(LVK, 0, 1, None)
# log(vee/k+ - 1) = log(vee/k+) + log(1 - k+/vee)
add(LVK1, 0, 1, L1MA)
# log((k+-Lambda)/Lambda) = log(k+/Lambda) + log(1 - Lambda/k+)
add(r"\log{\left(\frac{k^{+}-\Lambda}{\Lambda}\right)}", 1, 0, None, r"\log{\left(1-\frac{\Lambda}{k^{+}}\right)}")
add(r"\log{\left[\frac{k^{+}-\Lambda}{\Lambda}\right]}", 1, 0, None, r"\log{\left(1-\frac{\Lambda}{k^{+}}\right)}")
# log(Lambda/(k+-Lambda)) = -log(k+/Lambda) - log(1 - Lambda/k+)
add(r"\log{\left[\frac{\Lambda}{k^{+}-\Lambda}\right]}", -1, 0, None, r"\left(-\log{\left(1-\frac{\Lambda}{k^{+}}\right)}\right)")
# log(vee/(k+ + Lambda)) = log(vee/k+) - log(1 + Lambda/k+)
for s in (r"\log{\left[\frac{\vee}{\Lambda+k^{+}}\right]}", r"\log{\left[\frac{\vee}{k^{+}+\Lambda}\right]}",
          r"\log{\left[\frac{\vee}{\Lambda + k^{+}}\right]}"):
    add(s, 0, 1, None, r"\left(-\log{\left(1+\frac{\Lambda}{k^{+}}\right)}\right)")
# ---- logarithms with transverse distances (Group I Row III)
add(r"\log{\left[\frac{\left[\vee\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}{\left[(\Lambda + k^{+})\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}\right]}",
    0, 1, r"\log{\left[\frac{(x-z)^{2}-\frac{k^{+}}{\vee}\,(x-w)^{2}}{(x-z)^{2}-(x-w)^{2}}\right]}",
    r"\log{\left[\frac{(x-z)^{2}-\frac{k^{+}}{\vee}\,(x-w)^{2}}{\left(1+\frac{\Lambda}{k^{+}}\right)(x-z)^{2}-(x-w)^{2}}\right]}")
add(r"\log{\left[\frac{\left[(k^{+}-\Lambda)\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}{\left[\Lambda\, (x-z)^{2} - k^{+}\, (x-w)^{2}\right]}\right]}",
    0, 0, r"\log{\left[\frac{(x-w)^{2}-(x-z)^{2}}{(x-w)^{2}}\right]}",
    r"\log{\left[\frac{(k^{+}-\Lambda)\,(x-z)^{2}-k^{+}\,(x-w)^{2}}{\Lambda\,(x-z)^{2}-k^{+}\,(x-w)^{2}}\right]}")
add(r"\log{\left[\frac{\left[k^{+}\,(x-w)^{2}-\vee\, (x-z)^{2}\right]}{\left[k^{+}\,(x-w)^{2}-\Lambda\, (x-z)^{2}\right]}\right]}",
    0, 1, r"\log{\left[\frac{\frac{k^{+}}{\vee}\,(x-w)^{2}-(x-z)^{2}}{(x-w)^{2}}\right]}",
    r"\log{\left[\frac{\frac{k^{+}}{\vee}\,(x-w)^{2}-(x-z)^{2}}{(x-w)^{2}-\frac{\Lambda}{k^{+}}\,(x-z)^{2}}\right]}")
add(r"\log{\left[\frac{(x-w)^{2}\,k^{+}-(x-z)^{2}\,\vee}{(x-w)^{2}\,k^{+} - (x-z)^{2}\,(k^{+}+\Lambda)}\right]}",
    0, 1, r"\log{\left[\frac{(x-w)^{2}\,\frac{k^{+}}{\vee}-(x-z)^{2}}{(x-w)^{2}-(x-z)^{2}}\right]}",
    r"\log{\left[\frac{(x-w)^{2}\,\frac{k^{+}}{\vee}-(x-z)^{2}}{(x-w)^{2}-(x-z)^{2}\left(1+\frac{\Lambda}{k^{+}}\right)}\right]}")
add(r"\log{\left[\frac{(x-w)^{2}\,k^{+}-(x-z)^{2}\,(k^{+}-\Lambda)}{(x-w)^{2}\,k^{+}-(x-z)^{2}\,\Lambda}\right]}",
    0, 0, r"\log{\left[\frac{(x-w)^{2}-(x-z)^{2}}{(x-w)^{2}}\right]}",
    r"\log{\left[\frac{(x-w)^{2}\,k^{+}-(x-z)^{2}\,(k^{+}-\Lambda)}{(x-w)^{2}\,k^{+}-(x-z)^{2}\,\Lambda}\right]}")
add(r"\log{\left[\frac{\left[k^{+}\,(x-w)-(k^{+}-\Lambda)\,(x-z)\right]^{2}}{\left[k^{+}\,(x-w)-\Lambda\,(x-z)\right]^{2}}\right]}",
    0, 0, r"\log{\left[\frac{(z-w)^{2}}{(x-w)^{2}}\right]}",
    r"\log{\left[\frac{\left[k^{+}\,(x-w)-(k^{+}-\Lambda)\,(x-z)\right]^{2}}{\left[k^{+}\,(x-w)-\Lambda\,(x-z)\right]^{2}}\right]}")
# ---- log[((vee-k+) P + k+ Q)/(Lambda P + k+ Q)] (Groups I Row IV, II Row I, III Rows I, II)
for P, Q, forms in (
        (A, B, [r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\,(x-w)^{2}}\right]}",
                r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{\Lambda\, (x-z)^{2}+k^{+}\,(x-w)^{2}}\right]}",
                r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{k^{+}\,(x-w)^{2}+\Lambda\, (x-z)^{2}}\right]}",
                r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2} + k^{+}\, (x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\, (x-w)^{2}}\right]}",
                r"\log{\left[\frac{(\vee - k^{+})\, (x-z)^{2} + k^{+}\, (x-w)^{2}}{\Lambda\, (x-z)^{2} + k^{+}\, (x-w)^{2}}\right]}"]),
        (A, Bw, [r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w')^{2}}{\Lambda\, (x-z)^{2} + k^{+}\,(x-w')^{2}}\right]}"]),
        (Ap, Bp, [r"\log{\left[\frac{(\vee-k^{+})\,(x'-z)^{2}+k^{+}\,(x'-w')^{2}}{\Lambda\,(x'-z)^{2}+k^{+}\,(x'-w')^{2}}\right]}",
                  r"\log{\left[\frac{(\vee-k^{+})\,(x'-z)^{2} + k^{+}\, (x'-w')^{2}}{\Lambda\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}\right]}",
                  r"\log{\left[\frac{(\vee - k^{+})\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}{\Lambda\, (x'-z)^{2} + k^{+}\, (x'-w')^{2}}\right]}",
                  r"\log{\left[\frac{k^{+}\, (x'-w')^{2} + (\vee - k^{+})\, (x'-z)^{2}}{k^{+}\, (x'-w')^{2} + \Lambda\, (x'-z)^{2}}\right]}"])):
    ex, lim = lam_rem(P, Q)
    for s in forms:
        add(s, 0, 1, lim, ex)
# ---- Group II Row I 2a: vee cancels between numerator and denominator; the Lambda log has a finite limit
add(r"\log{\left[\frac{(\vee-k^{+})\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{k^{+}\, (x'-w')^{2} + (\vee-k^{+})\, (x'-z)^{2}}\right]}", 0, 0,
    r"\log{\left[\frac{\left(1-\frac{k^{+}}{\vee}\right)(x-z)^{2}+\frac{k^{+}}{\vee}\,(x-w)^{2}}{\frac{k^{+}}{\vee}\, (x'-w')^{2} + \left(1-\frac{k^{+}}{\vee}\right) (x'-z)^{2}}\right]}")
add(r"\log{\left[\frac{\Lambda\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{k^{+}\, (x'-w')^{2} + \Lambda\, (x'-z)^{2}}\right]}", 0, 0,
    r"\log{\left[\frac{(x-w)^{2}}{(x'-w')^{2}}\right]}",
    r"\log{\left[\frac{\Lambda\,(x-z)^{2}+k^{+}\,(x-w)^{2}}{k^{+}\, (x'-w')^{2} + \Lambda\, (x'-z)^{2}}\right]}")

SPLIT = {}
for s, cL, cV, R, Rx in _raw:
    SPLIT[norm(s)] = (cL, cV, R, Rx)
# logs that need no splitting (finite, no Lambda, no vee or only through k+/vee)
FINITE = [norm(s) for s in (
    r"\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}", r"\log{\left[1-\frac{k^{+}}{\vee}\right]}", L1MA,
    r"\log{\left(1-\frac{k^{+}}{\vee}\right)}",
    r"\log{\left[\frac{(x-w)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right]}",
    r"\log{\left(\frac{(y'-y)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right)}",
    r"\log\frac{(y-z)^{2}\mu^{2}_{\overline{MS}}}{4}",
    r"\log\Big(1-\frac{\Lambda}{k^{+}}\Big)")]
