#!/usr/bin/env python3
"""Typeset the draft error review as a PDF.

    pip install reportlab
    python3 review/pdf/build_pdf.py [output.pdf]
"""

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "tools"))

from reportlab.platypus import PageBreak, Spacer          # noqa: E402
import mdpdf                                              # noqa: E402

SRC = os.path.join(REPO, "review", "draft-review-omega-g2.md")
OUT_DEFAULT = os.path.join(REPO, "review", "Omega_g2_draft_corrections.pdf")

TITLE = ("Corrections to the <i>Omega at g²</i> Draft")

SUBTITLES = [
    "Located errors, internal inconsistencies and conceptual gaps in",
    "<i>Soft Gluon Wave Function and Evolution Operator in the CGC "
    "at Next-to-Leading Order</i>",
    "R. Radhakrishnan &nbsp;·&nbsp; draft reviewed: 57 pp., prepared for JHEP",
]

BLURBS = [
    "<b>Method.</b> Every item gives the equation number and the specific object. "
    "Section A lists errors verified by explicit re-derivation from the paper's own "
    "matrix elements (D.1)–(D.7) and its own conventions. Section B lists internal "
    "contradictions between two places in the paper, and says which one I believe is "
    "right. Section C is conceptual; section D is typography and prose. What was "
    "checked and found <i>correct</i> is listed first, so the verified ground is explicit.",

    "<b>The one that breaks a stated result.</b> Eq. (6.4) prints the two "
    "A-terms with opposite signs, making the bracket anti-Hermitian — so it cannot "
    "cancel the Hermitian H<sub>g</sub>, and Eq. (6.5) does not follow. Two independent "
    "arguments fix the sign as +A(p²/2p⁺)a<super>†</super>.",

    "<b>A coefficient of Ω is wrong.</b> Eq. (4.55) is missing the factor i from the "
    "colour-charge commutator and the factor 1/√(k⁺p⁺). The error propagates into "
    "Eq. (5.4), where the two terms inside one bracket then differ dimensionally.",

    "<b>Also located.</b> Eq. (4.36) omits the instantaneous piece that Eq. (4.35) "
    "promises and Eq. (5.3) actually uses; a stray momentum k and a swapped p↔q "
    "assignment in Eq. (5.3); colour factor and ρ argument in Eq. (4.86); g⁴ for g² in "
    "Eq. (4.94); unbalanced brackets in Eq. (2.23); a (2π)⁹ measure mismatch between "
    "Eqs. (4.97) and (4.99); and two duplicate bibliography entries.",
]

FOOTNOTE = ("Prepared with Claude Code. Section A items were verified by re-derivation; "
            "section B items are internal contradictions requiring an authorial decision.")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else OUT_DEFAULT
    mdpdf.register_fonts()
    S = mdpdf.make_styles()

    doc = mdpdf.Doc(out,
                    header="Corrections to the Omega at g² draft",
                    footer="arXiv:2607.18373 — draft review")
    avail = doc.width

    with open(SRC, encoding="utf-8") as fh:
        md = fh.read()

    marks = []
    body = mdpdf.render(md, S, avail, heading_sink=marks)

    # contents: level-2 sections with their level-3 items nested underneath
    entries, cur = [], None
    for lvl, txt in marks:
        if lvl == 1:
            continue
        if lvl == 2:
            cur = (txt, [])
            entries.append(cur)
        elif lvl == 3 and cur is not None:
            cur[1].append(txt)

    story = mdpdf.cover(S, TITLE, SUBTITLES, BLURBS, FOOTNOTE, top=26)
    story.append(PageBreak())
    story += mdpdf.contents(S, entries, "Index of findings")
    story.append(PageBreak())
    story += body
    story.append(Spacer(1, 2))

    doc.build(story)
    print("wrote %s (%.1f KB)" % (out, os.path.getsize(out) / 1024.0))


if __name__ == "__main__":
    main()
