#!/usr/bin/env python3
"""Shared markdown -> PDF renderer used by the project build scripts.

Handles headings, pipe tables, fenced code, blockquotes, lists, inline spans
and links, with a cover page, an optional contents page, and running page
furniture.

Uses TrueType fonts throughout so that Greek letters, subscripts, arrows,
box-drawing characters and mathematical symbols render as real glyphs rather
than the black boxes ReportLab's built-in WinAnsi fonts would produce.
FreeSans is the body family (the only installed sans with a complete
regular/bold/oblique set); DejaVu Sans Mono is used for code (best
box-drawing and math coverage).

Consumers:  nlo-multiplicity/pdf/build_pdf.py, review/pdf/build_pdf.py
"""

import html
import os
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, HRFlowable, PageBreak, PageTemplate,
    Paragraph, Preformatted, Spacer, Table, TableStyle,
)

SANS_DIR = "/usr/share/fonts/truetype/freefont"   # complete regular/bold/oblique family
MONO_DIR = "/usr/share/fonts/truetype/dejavu"     # best box-drawing + math coverage
ACCENT = colors.HexColor("#1f3a5f")
RULE = colors.HexColor("#c8d2de")
CODE_BG = colors.HexColor("#f4f6f9")
MUTED = colors.HexColor("#5a6472")


def register_fonts():
    faces = {
        (SANS_DIR, "DJV"): "FreeSans.ttf",
        (SANS_DIR, "DJV-Bold"): "FreeSansBold.ttf",
        (SANS_DIR, "DJV-Oblique"): "FreeSansOblique.ttf",
        (SANS_DIR, "DJV-BoldOblique"): "FreeSansBoldOblique.ttf",
        (MONO_DIR, "DJVMono"): "DejaVuSansMono.ttf",
        (MONO_DIR, "DJVMono-Bold"): "DejaVuSansMono-Bold.ttf",
        (MONO_DIR, "DJVMono-Oblique"): "DejaVuSansMono-Oblique.ttf",
    }
    for (d, name), fn in faces.items():
        pdfmetrics.registerFont(TTFont(name, os.path.join(d, fn)))
    pdfmetrics.registerFontFamily(
        "DJV", normal="DJV", bold="DJV-Bold",
        italic="DJV-Oblique", boldItalic="DJV-BoldOblique")
    pdfmetrics.registerFontFamily(
        "DJVMono", normal="DJVMono", bold="DJVMono-Bold",
        italic="DJVMono-Oblique", boldItalic="DJVMono-Oblique")


def make_styles():
    ss = getSampleStyleSheet()
    S = {}
    S["body"] = ParagraphStyle(
        "body", parent=ss["Normal"], fontName="DJV", fontSize=9.1, leading=13.6,
        alignment=TA_JUSTIFY, spaceAfter=6.5, textColor=colors.HexColor("#12161c"))
    S["h1"] = ParagraphStyle(
        "h1", parent=S["body"], fontName="DJV-Bold", fontSize=17, leading=21,
        textColor=ACCENT, spaceBefore=4, spaceAfter=10, alignment=0)
    S["h2"] = ParagraphStyle(
        "h2", parent=S["body"], fontName="DJV-Bold", fontSize=12.6, leading=16.5,
        textColor=ACCENT, spaceBefore=15, spaceAfter=6, alignment=0)
    S["h3"] = ParagraphStyle(
        "h3", parent=S["body"], fontName="DJV-Bold", fontSize=10.4, leading=14,
        textColor=colors.HexColor("#2c4c74"), spaceBefore=11, spaceAfter=4, alignment=0)
    S["h4"] = ParagraphStyle(
        "h4", parent=S["body"], fontName="DJV-BoldOblique", fontSize=9.4, leading=13,
        textColor=colors.HexColor("#3d566f"), spaceBefore=9, spaceAfter=3, alignment=0)
    S["code"] = ParagraphStyle(
        "code", parent=ss["Code"], fontName="DJVMono", fontSize=7.1, leading=9.2,
        textColor=colors.HexColor("#1a2b3c"), leftIndent=0, spaceBefore=0, spaceAfter=0)
    S["bullet"] = ParagraphStyle(
        "bullet", parent=S["body"], leftIndent=13, bulletIndent=3, spaceAfter=3.5)
    S["bullet2"] = ParagraphStyle(
        "bullet2", parent=S["bullet"], leftIndent=27, bulletIndent=17)
    S["quote"] = ParagraphStyle(
        "quote", parent=S["body"], leftIndent=14, rightIndent=8,
        fontName="DJV-Oblique", textColor=colors.HexColor("#25405e"),
        borderPadding=0, spaceBefore=5, spaceAfter=7)
    S["tcell"] = ParagraphStyle(
        "tcell", parent=S["body"], fontSize=7.6, leading=10.2,
        alignment=0, spaceAfter=0)
    S["thead"] = ParagraphStyle(
        "thead", parent=S["tcell"], fontName="DJV-Bold", textColor=colors.white)
    S["title"] = ParagraphStyle(
        "title", parent=S["body"], fontName="DJV-Bold", fontSize=23, leading=29,
        alignment=TA_CENTER, textColor=ACCENT, spaceAfter=10)
    S["subtitle"] = ParagraphStyle(
        "subtitle", parent=S["body"], fontSize=11.5, leading=16,
        alignment=TA_CENTER, textColor=MUTED, spaceAfter=6)
    S["tiny"] = ParagraphStyle(
        "tiny", parent=S["body"], fontSize=8, leading=11.5,
        alignment=TA_CENTER, textColor=MUTED)
    S["toc"] = ParagraphStyle(
        "toc", parent=S["body"], fontSize=9.6, leading=15, alignment=0, spaceAfter=1)
    S["toc2"] = ParagraphStyle(
        "toc2", parent=S["toc"], fontSize=8.3, leading=12)
    S["coverbox"] = ParagraphStyle(
        "coverbox", parent=S["body"], fontSize=8.7, leading=13,
        leftIndent=10 * mm, rightIndent=10 * mm, spaceAfter=7)
    return S


# ---------------------------------------------------------------- inline spans

def inline(text):
    """Markdown inline formatting -> ReportLab mini-HTML."""
    out, i, n = [], 0, len(text)
    while i < n:
        ch = text[i]
        if ch == "`":                                   # code span
            j = text.find("`", i + 1)
            if j == -1:
                out.append(html.escape(ch)); i += 1; continue
            body = html.escape(text[i + 1:j])
            out.append('<font face="DJVMono" size="8.2" '
                       'color="#8a2f5f">%s</font>' % body)
            i = j + 1
            continue
        if text.startswith("**", i):                    # bold
            j = text.find("**", i + 2)
            if j == -1:
                out.append(html.escape(text[i:i + 2])); i += 2; continue
            out.append("<b>%s</b>" % inline(text[i + 2:j]))
            i = j + 2
            continue
        if ch == "[":                                   # link
            m = re.match(r"\[([^\]]*)\]\(([^)\s]+)[^)]*\)", text[i:])
            if m:
                label = inline(m.group(1))
                href = html.escape(m.group(2), quote=True)
                out.append('<a href="%s" color="#17559c">%s</a>' % (href, label))
                i += m.end()
                continue
        if ch == "*" and not text.startswith("**", i):  # italic
            j = text.find("*", i + 1)
            if j != -1 and j > i + 1:
                out.append("<i>%s</i>" % inline(text[i + 1:j]))
                i = j + 1
                continue
        out.append(html.escape(ch))
        i += 1
    return "".join(out)


# --------------------------------------------------------------- block parsing

def split_row(line):
    """Split a pipe-table row, honouring backtick code spans (which may
    themselves contain `|`, e.g. `|α|²`)."""
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    cells, buf, esc, incode = [], [], False, False
    for ch in line:
        if esc:
            buf.append(ch); esc = False
        elif ch == "\\":
            esc = True
        elif ch == "`":
            incode = not incode; buf.append(ch)
        elif ch == "|" and not incode:
            cells.append("".join(buf).strip()); buf = []
        else:
            buf.append(ch)
    cells.append("".join(buf).strip())
    return cells


def is_divider(line):
    return bool(re.fullmatch(r"\|?[\s:|-]*-[\s:|-]*\|?", line.strip())) and "-" in line


def build_table(rows, S, avail):
    header, body = rows[0], rows[1:]
    ncol = max(len(r) for r in rows)
    header += [""] * (ncol - len(header))
    data = [[Paragraph(inline(c), S["thead"]) for c in header]]
    for r in body:
        r = r + [""] * (ncol - len(r))
        data.append([Paragraph(inline(c), S["tcell"]) for c in r])

    weights = []
    for c in range(ncol):
        longest = max((len(rows[r][c]) if c < len(rows[r]) else 0)
                      for r in range(len(rows)))
        weights.append(max(6.0, min(float(longest), 60.0)))
    total = sum(weights)
    widths = [avail * w / total for w in weights]

    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    style = [
        ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.4, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]
    for i in range(1, len(data)):
        if i % 2 == 0:
            style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f7f9fb")))
    t.setStyle(TableStyle(style))
    return t


def code_block(lines, S, avail):
    txt = "\n".join(lines) if lines else " "
    fs = 7.1
    widest = max((pdfmetrics.stringWidth(l, "DJVMono", fs) for l in lines), default=0)
    inner = avail - 12
    if widest > inner and widest > 0:
        fs = max(4.6, fs * inner / widest)
    st = ParagraphStyle("cb", parent=S["code"], fontSize=fs, leading=fs * 1.30)
    pre = Preformatted(txt, st)
    t = Table([[pre]], colWidths=[avail], hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), CODE_BG),
        ("BOX", (0, 0), (-1, -1), 0.5, RULE),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ]))
    return t


def render(md, S, avail, heading_sink=None):
    """Markdown string -> list of flowables.  If heading_sink is a list, it
    collects (level, text) for every heading, for building a contents page."""
    flow = []
    lines = md.split("\n")
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith("```"):                       # fenced code
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i].rstrip())
                i += 1
            i += 1
            flow.append(Spacer(1, 3))
            flow.append(code_block(buf, S, avail))
            flow.append(Spacer(1, 7))
            continue

        if not stripped:
            i += 1
            continue

        if re.fullmatch(r"(---+|\*\*\*+|___+)", stripped):   # rule
            flow.append(Spacer(1, 4))
            flow.append(HRFlowable(width="100%", thickness=0.7, color=RULE,
                                   spaceBefore=2, spaceAfter=8))
            i += 1
            continue

        m = re.match(r"(#{1,6})\s+(.*)", stripped)           # heading
        if m:
            lvl = len(m.group(1))
            key = {1: "h1", 2: "h2", 3: "h3"}.get(lvl, "h4")
            txt = m.group(2).strip()
            if heading_sink is not None:
                heading_sink.append((lvl, re.sub(r"[*`]", "", txt)))
            flow.append(Paragraph(inline(txt), S[key]))
            i += 1
            continue

        if stripped.startswith("|") and i + 1 < n and is_divider(lines[i + 1]):
            rows = [split_row(stripped)]
            i += 2
            while i < n and lines[i].strip().startswith("|"):
                rows.append(split_row(lines[i]))
                i += 1
            flow.append(Spacer(1, 3))
            flow.append(build_table(rows, S, avail))
            flow.append(Spacer(1, 9))
            continue

        if stripped.startswith(">"):                          # blockquote
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip())
                i += 1
            para = " ".join(x for x in buf if x)
            flow.append(Table(
                [[Paragraph(inline(para), S["quote"])]], colWidths=[avail],
                hAlign="LEFT",
                style=TableStyle([
                    ("LINEBEFORE", (0, 0), (0, -1), 2.2, ACCENT),
                    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#f3f6fa")),
                    ("LEFTPADDING", (0, 0), (-1, -1), 9),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                    ("TOPPADDING", (0, 0), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 5)])))
            flow.append(Spacer(1, 7))
            continue

        m = re.match(r"^(\s*)([-*+]|\d+[.)])\s+(.*)", line)   # list item
        if m:
            indent = len(m.group(1))
            marker = m.group(2)
            body = [m.group(3)]
            i += 1
            while i < n:
                nxt = lines[i]
                if not nxt.strip():
                    break
                if re.match(r"^\s*([-*+]|\d+[.)])\s+", nxt) or \
                   re.match(r"\s*(#{1,6}\s|```|\||>)", nxt):
                    break
                body.append(nxt.strip())
                i += 1
            style = S["bullet2"] if indent >= 2 else S["bullet"]
            bullet = marker if re.match(r"\d", marker) else "•"
            flow.append(Paragraph(inline(" ".join(body)), style, bulletText=bullet))
            continue

        buf = [stripped]                                      # paragraph
        i += 1
        while i < n and lines[i].strip() and not re.match(
                r"\s*(#{1,6}\s|```|\||>|[-*+]\s|\d+[.)]\s|---)", lines[i]):
            buf.append(lines[i].strip())
            i += 1
        flow.append(Paragraph(inline(" ".join(buf)), S["body"]))
    return flow


# ------------------------------------------------------------- page furniture

class Doc(BaseDocTemplate):
    def __init__(self, path, header="", footer="", **kw):
        super().__init__(path, pagesize=A4,
                         leftMargin=19 * mm, rightMargin=17 * mm,
                         topMargin=18 * mm, bottomMargin=17 * mm, **kw)
        self.header_text = header
        self.footer_text = footer
        frame = Frame(self.leftMargin, self.bottomMargin,
                      self.width, self.height, id="main",
                      leftPadding=0, rightPadding=0,
                      topPadding=0, bottomPadding=0)
        self.addPageTemplates([
            PageTemplate(id="cover", frames=[frame]),
            PageTemplate(id="main", frames=[frame], onPage=self.decorate),
        ])

    def decorate(self, canv, doc):
        canv.saveState()
        canv.setFont("DJV", 7.2)
        canv.setFillColor(MUTED)
        canv.drawString(self.leftMargin, self.pagesize[1] - 12 * mm, self.header_text)
        canv.setStrokeColor(RULE)
        canv.setLineWidth(0.4)
        canv.line(self.leftMargin, self.pagesize[1] - 13.6 * mm,
                  self.pagesize[0] - self.rightMargin, self.pagesize[1] - 13.6 * mm)
        canv.line(self.leftMargin, self.bottomMargin - 4 * mm,
                  self.pagesize[0] - self.rightMargin, self.bottomMargin - 4 * mm)
        canv.drawRightString(self.pagesize[0] - self.rightMargin,
                             self.bottomMargin - 8.5 * mm, str(doc.page))
        canv.drawString(self.leftMargin, self.bottomMargin - 8.5 * mm, self.footer_text)
        canv.restoreState()


def cover(S, title, subtitles=(), blurbs=(), footnote=None, top=34 * mm):
    """Standard cover page flowables (no trailing PageBreak)."""
    flow = [Spacer(1, top), Paragraph(title, S["title"]), Spacer(1, 4 * mm)]
    for s in subtitles:
        flow.append(Paragraph(s, S["subtitle"]))
    if blurbs:
        flow += [Spacer(1, 12 * mm),
                 HRFlowable(width="55%", thickness=0.8, color=RULE, hAlign="CENTER"),
                 Spacer(1, 10 * mm)]
    for b in blurbs:
        flow.append(Paragraph(b, S["coverbox"]))
    if footnote:
        flow += [Spacer(1, 14 * mm), Paragraph(footnote, S["tiny"])]
    return flow


def contents(S, entries, heading="Contents"):
    """entries: list of (top_label, [sub_labels])."""
    flow = [Paragraph(heading, S["h1"]),
            HRFlowable(width="100%", thickness=0.7, color=RULE, spaceAfter=9)]
    for label, subs in entries:
        flow.append(Paragraph("<b>%s</b>" % html.escape(label), S["toc"]))
        for s in subs:
            flow.append(Paragraph(
                "&nbsp;&nbsp;&nbsp;&nbsp;<font color='#5a6472'>%s</font>"
                % html.escape(s), S["toc2"]))
        flow.append(Spacer(1, 4))
    return flow
