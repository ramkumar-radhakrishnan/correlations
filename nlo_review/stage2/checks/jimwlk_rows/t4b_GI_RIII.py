import sys; sys.path.insert(0,".")
from t4_brackets import *
W=World(41); W2=World(43)
lst=[]
for l in b3(207,211): lst += list_monos_h(B3[l], True)
orig = BR["GI  RIII"]
# variant: second part with U^{bc}(w) instead of U^{bc}(z) in its first term (as in the fraction form, 3rd term of line 60)
BR["GI  RIII"] = (True, [orig[1][0],
    (-2.0, [("x'-w'", "y-w"), ("x-z", "w-z")], r"f^{abd}" + FB + r"[U^{de}(x)U^{bc}(w)\rho^{e}(x)\rho^{c}(y)-U^{de}(x)U^{bc}(y)\rho^{e}(x)\rho^{c}(y)]")])
print("GI RIII with U^{bc}(w) in N_2: ordered max|diff| = %.2e / %.2e" % (odiff(build("GI  RIII"), lst, W)[0], odiff(build("GI  RIII"), lst, W2)[0]))
# alpha-part and beta-part separately
for tag, ent in (("alpha/alpha'", [orig[1][0]]), ("beta/beta' (U(w))", [BR["GI  RIII"][1][1]]), ("beta/beta' as typed (U(z))", [orig[1][1]])):
    BR["tmp"] = (True, ent)
    sel = [m for m in lst if m.tag.startswith("alpha" if tag.startswith("alpha") else "beta")]
    print("   %-28s: %.2e" % (tag, odiff(build("tmp"), sel, W)[0]))
