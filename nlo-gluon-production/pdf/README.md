# Building the PDF

`document.html` is the typeset source. KaTeX renders the math and headless Chromium
prints it to PDF — no LaTeX toolchain required.

```bash
npm install katex                       # provides node_modules/katex/dist/*
pip install playwright
python3 make_figure.py                  # -> figure_factor.svg
python3 build_pdf.py                    # -> ../NLO_gluon_production_DGLAP.pdf
```

`document.html` keeps a `<!--FIGURE-->` placeholder that `build_pdf.py` replaces with
`figure_factor.svg` at build time, so the HTML stays the single editable source.

`build_pdf.py` points at the pre-installed Chromium via `executable_path` and fails
loudly on any KaTeX parse error or unrendered `$...$` fragment. Note that a literal
`<` inside math must be written `&lt;` — otherwise the browser reads it as a tag and
the math silently fails to render.

The figure palette is the data-viz reference palette, slots 1–3, validated for
colour-vision deficiency; the curves carry direct labels rather than relying on
colour alone.
