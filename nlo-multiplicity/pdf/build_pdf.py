#!/usr/bin/env python3
"""Typeset the nine-part analysis (plus the JHEP addendum) as a single PDF.

    pip install reportlab
    python3 nlo-multiplicity/pdf/build_pdf.py [output.pdf]
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
REPO = os.path.dirname(PROJ)
sys.path.insert(0, os.path.join(REPO, "tools"))

from reportlab.platypus import PageBreak                  # noqa: E402
import mdpdf                                              # noqa: E402

OUT_DEFAULT = os.path.join(PROJ, "NLO_CGC_multiplicity_analysis.pdf")

ORDER = [
    ("00-scope-and-caveats.md", "Scope, notation, and honest caveats"),
    ("01-formalism.md", "Part I — The formalism"),
    ("02-moments.md", "Part II — Moments of the multiplicity operator"),
    ("03-cumulants.md", "Part III — Connected correlations and cumulants"),
    ("04-correlations.md", "Part IV — Correlation functions"),
    ("05-models.md", "Part V — Model calculations"),
    ("06-numerics.md", "Part VI — Numerical studies"),
    ("07-literature.md", "Part VII — Literature survey"),
    ("08-publication.md", "Part VIII — Publication potential"),
    ("09-roadmap.md", "Part IX — Research roadmap"),
    ("10-jhep-assessment.md", "Addendum — Will the paper itself make it into JHEP?"),
]

TITLE = ("Gluon Multiplicity Moments and<br/>Particle-Number Fluctuations "
         "at Order g⁴")

SUBTITLES = [
    "A formalism, feasibility and publication analysis of the programme announced in",
    "<i>Soft Gluon Wave Function and Evolution Operator in the CGC "
    "at Next-to-Leading Order</i>",
    '<a href="https://arxiv.org/abs/2607.18373" color="#17559c">arXiv:2607.18373</a>'
    " &nbsp;·&nbsp; R. Radhakrishnan",
]

BLURBS = [
    "<b>Three structural results.</b> (1) Because Ω is written fully normal-ordered, "
    "its C, D and E structures annihilate the soft vacuum, so Ω|0⟩ is exactly a "
    "truncated displaced squeezed state through O(g²) — and every moment follows "
    "from a Fredholm determinant. "
    "(2) Because ⟨aa⟩<sub>c</sub> is O(g²) while ⟨a<super>†</super>a⟩<sub>c</sub> is "
    "O(g⁴), the leading non-Poissonian effect is displacement–squeezing interference, "
    "and C₂ is complete at O(g⁴) with the existing operator. "
    "(3) Bose enhancement peaks near-side; pair emission from the squeezing kernel B₂ "
    "peaks away-side — the only clean discriminator against classical colour-charge "
    "fluctuations, which are degenerate in both g and 1/N<sub>dof</sub>.",

    "<b>Two limits stated honestly.</b> C₃ and C₄ are <i>not</i> complete — "
    "the three-gluon component of Ω|0⟩ first appears at O(g³) and contributes "
    "at the same order. And the classical and quantum contributions cannot be "
    "separated by power counting.",

    "<b>Verdict.</b> One PRD paper, ~70% confidence, conditional on a two-week "
    "go/no-go test. Not JHEP as scoped. More than a note.",
]

FOOTNOTE = ("Prepared with Claude Code. Parts I–IX were first drafted without access "
            "to the paper and were subsequently reconciled against v2; the first section "
            "records exactly what changed.")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else OUT_DEFAULT
    mdpdf.register_fonts()
    S = mdpdf.make_styles()

    doc = mdpdf.Doc(out,
                    header="Gluon multiplicity moments from the NLO CGC evolution operator",
                    footer="arXiv:2607.18373 follow-up analysis")
    avail = doc.width

    bodies, entries = [], []
    for fn, label in ORDER:
        with open(os.path.join(PROJ, "docs", fn), encoding="utf-8") as fh:
            md = fh.read()
        marks = []
        bodies.append(mdpdf.render(md, S, avail, heading_sink=marks))
        entries.append((label, [t for lvl, t in marks if lvl == 2][:9]))

    story = mdpdf.cover(S, TITLE, SUBTITLES, BLURBS, FOOTNOTE)
    story.append(PageBreak())
    story += mdpdf.contents(S, entries)
    for flow in bodies:
        story.append(PageBreak())
        story.extend(flow)

    doc.build(story)
    print("wrote %s (%.1f KB)" % (out, os.path.getsize(out) / 1024.0))


if __name__ == "__main__":
    main()
