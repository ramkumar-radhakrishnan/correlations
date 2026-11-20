"""Split Group I row V (6.9)-(6.11) and Group III (8.1)-(8.5) into kernels grouped by colour factor.
Each entry: (label, colour_group, coefficient, measure, integration variable, python kernel).
coefficient is relative to (1/(2pi)^3) g^4/pi^4, without colour constants (f^{c'd'a} etc.)."""
import re, texeval as te
from pconv import to_py, clean

B = te.blocks('../Final_terms_original.tex')
def norm(b):
    b = b.replace(r'{\color{red},}', ',').replace(r'w^{\prime}', "w'").replace(r'w^\prime', "w'")
    return b
POS = re.compile(r"\\int_\{(?:[xyzw]'?,?\s*)+\}")
def after_pos(seg, n=0):
    ms = list(POS.finditer(seg)); return seg[ms[n].end():]
def cut(seg, marker):
    i = seg.find(marker); assert i >= 0, marker; return seg[:i]

T = []
# ---------------- (6.9)
b = norm(B[18])
s1 = cut(after_pos(b, 0), r"\times\, \left[U^{ab'}(x')")
s2 = cut(after_pos(b, 1), r"\notag \\ & \times\, \left[U^{ab'}(x')")
T += [('6.9a', 'C69', 1j/16, '1/K^2', 'p', to_py(s1)), ('6.9b', 'C69', -1j/16, '1/K', 'q', to_py(s2))]
# ---------------- (6.10)
b = norm(B[19])
s1 = cut(after_pos(b, 0), r"\Bigg[U^{ab'}(x')")
s2 = cut(after_pos(b, 1), r"\Bigg[U^{ab'}(x')")
T += [('6.10a', 'C610', 1/8, '1/K', 'p', to_py(s1)), ('6.10b', 'C610', -1/8, '1/K', 'q', to_py(s2))]   # int_{k+}^{0} -> -int_0^{k+}
# ---------------- (6.11)
b = norm(B[20])
s1 = cut(after_pos(b, 0), r"\times\Big[U^{ab'}(x')")
T += [('6.11', 'C611', 1/8, '1/K', 'p', to_py(s1))]

# ---------------- Group III, measure int_0^{vee-k+} dp/(2pi)
A_ = r"\frac{(x'-y')^{k'}}{(x'-y')^{2}}"
B_ = r"\frac{(x'-z)^{k'}}{(x'-z)^{2}}"
def two_part(seg, endmark_first, label, g1, g2, coef, meas):
    out = []
    i = seg.find(A_); base = seg[:i]
    base = base[:base.rfind('[')]     # drop the opening bracket in front of A_
    out.append((label + 'b1', g1, coef, meas, 'p', to_py(base + A_)))
    out.append((label + 'b2', g2, -coef, meas, 'p', to_py(base + B_)))
    return out
# (8.2)
b = norm(B[31])
s1 = cut(after_pos(b, 0), r"\left[U^{c'e'}(y')")
T += [('8.2a', 'X1', 1j/8, '1/K', 'p', to_py(s1))]
T += two_part(after_pos(b, 1), None, '8.2', 'X1', 'X2', 1j/8, '1/K')
# (8.3)
b = norm(B[32])
s1 = cut(after_pos(b, 0), r"\left[f^{c'd'a}")
T += [('8.3a', 'Y1', 1/4, '1', 'p', to_py(s1))]
seg2 = after_pos(b, 1); seg2 = seg2[:seg2.find(r"\left[f^{cbd}")]
T += two_part(seg2, None, '8.3', 'Y1', 'Y2', 1/4, '1')
s3 = cut(after_pos(b, 2), r"\left[f^{c'd'a}").replace("w'", "y'")      # relabel w' -> y'
T += [('8.3c', 'Y2', 1/4, '1', 'p', to_py(s3))]
# (8.4)
b = norm(B[33])
s1 = cut(after_pos(b, 0), r"\left[U^{c'e'}(y')")
T += [('8.4a', 'W1', 1j/4, '1/K', 'p', to_py(s1))]
T += two_part(after_pos(b, 1), None, '8.4', 'W1', 'W2', 1j/4, '1/K')
# (8.5)
b = norm(B[34])
s1 = cut(after_pos(b, 0), r"\left[U^{c'e'}(y')")
T += [('8.5a', 'Z1', 1j/4, '1/K', 'p', to_py(s1))]
T += two_part(after_pos(b, 1), None, '8.5', 'Z1', 'Z2', 1j/4, '1/K')

if __name__ == '__main__':
    for lab, g, c, meas, var, k in T:
        print(f'[{lab:6s}] colour {g:5s} coef {c}  meas {meas} d{var}:\n    {k}\n')
