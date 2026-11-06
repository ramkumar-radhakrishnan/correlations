"""Block extraction from the uploaded tex files."""
import re, math

UP = "/root/.claude/uploads/8b1e3540-8e36-5e91-9944-396d8d8beed8/"
F3 = UP + "b51ab81a-Evolution_three_rho_corrected.tex"
F2 = UP + "bc6b127f-Evolution_two_rho_corrected.tex"
FJ = UP + "dee634f9-Evolution_JIMWLK__terms_corrected.tex"
FB = UP + "bf616470-BFKL_corrected.tex"
FO = UP + "2a0b5ba5-Other_terms_corrected.tex"
F4 = UP + "f037673c-4rho_corrected.tex"

def read_lines(path):
    with open(path) as fh:
        return fh.readlines()

def block_body(path, l0, l1, marker=r"\Bigg["):
    """Return text of lines [l0..l1] (1-based inclusive), body after first marker
    following 'e^{-ik', trailing \Bigg] / \end{align} stripped."""
    txt = "".join(read_lines(path)[l0 - 1:l1])
    txt = txt.replace("\\-", "-")
    i = txt.find("e^{-ik")
    j = txt.find(marker, i if i >= 0 else 0)
    assert j >= 0, (path, l0, "no marker")
    body = txt[j + len(marker):]
    body = body.replace("\\end{align}", " ")
    # strip trailing \Bigg]
    k = body.rfind("\\Bigg]")
    if k >= 0:
        body = body[:k] + body[k + 6:]
    return body

PI = math.pi
