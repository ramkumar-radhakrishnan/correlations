# Building the PDF

`document.html` is the typeset source; KaTeX renders the math and headless
Chromium prints it to PDF (no LaTeX toolchain required).

```bash
npm install katex                       # provides node_modules/katex/dist/*
pip install playwright
python3 build_pdf.py                    # -> ../CSL_lattice_Hamiltonian_general.pdf
```

`build_pdf.py` points at the pre-installed Chromium via `executable_path` and
fails loudly on any KaTeX parse error or unrendered `$...$` fragment.
The content mirrors `../DERIVATION.md`.
