#!/bin/bash
# Builds main.pdf: copies the sources to a build directory (BUILD_DIR, default .build next to the sources),
# runs pdflatex [/bibtex/pdflatex/pdflatex with "full"], lists LaTeX errors and overfull boxes > 5pt,
# and copies main.pdf back.   Usage: bash scripts/build.sh [full]
D="$(cd "$(dirname "$0")/.." && pwd)"
B="${BUILD_DIR:-$D/.build}"
mkdir -p "$B/figures" "$B/rows"
cp "$D"/*.tex "$D"/*.bib "$D"/*.sty "$D"/*.bst "$D"/*.tikzstyles "$B"/
cp "$D"/figures/* "$B"/figures/
cp "$D"/rows/*.tex "$B"/rows/
cd "$B" || exit 1
pdflatex -interaction=nonstopmode main.tex > log1.txt 2>&1
if [ "$1" == "full" ]; then bibtex main > bib.txt 2>&1; pdflatex -interaction=nonstopmode main.tex > log2.txt 2>&1; pdflatex -interaction=nonstopmode main.tex > log3.txt 2>&1; fi
grep -n '^!' main.log | head -5
grep -c 'undefined' main.log | sed 's/^/undefined references or citations (lines in log): /'
python3 -I - "$B/main.log" <<'PY'
import re, sys
log = open(sys.argv[1], errors='ignore').read()
cur = None
for m in re.finditer(r'\((\./[\w\-./]+\.(?:tex|tikz))|Overfull \\hbox \(([\d.]+)pt too wide\) (?:in paragraph at lines|detected at line) (\d+)', log):
    if m.group(1):
        if 'figures' not in m.group(1) and 'tikzstyles' not in m.group(1): cur = m.group(1)
    elif float(m.group(2)) > 5: print('overfull', cur, m.group(3), m.group(2))
PY
grep 'Output written' main.log
cp main.pdf "$D/main.pdf"
