# Building the PDF

`document.html` is the typeset source; KaTeX renders the math and headless Chromium
prints it to PDF (no LaTeX toolchain required).

```bash
npm install katex          # provides node_modules/katex/dist/*
pip install playwright
python3 build_pdf.py       # -> ../z_integration_UV_finite_regular.pdf
```

`build_pdf.py` fails loudly on any KaTeX parse error or unrendered `$...$` fragment.
