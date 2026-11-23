#!/bin/bash
B=/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/paperbuild/build_v2
mkdir -p $B/figures
cp /home/user/correlations/nlo_review/paper/jhep_v2/*.tex /home/user/correlations/nlo_review/paper/jhep_v2/*.bib /home/user/correlations/nlo_review/paper/jhep_v2/*.sty /home/user/correlations/nlo_review/paper/jhep_v2/*.bst /home/user/correlations/nlo_review/paper/jhep_v2/*.tikzstyles $B/
cp /home/user/correlations/nlo_review/paper/jhep_v2/figures/* $B/figures/
mkdir -p $B/rows; cp /home/user/correlations/nlo_review/paper/jhep_v2/rows/*.tex $B/rows/
cd $B
pdflatex -interaction=nonstopmode main.tex > log1.txt 2>&1
if [ "$1" == "full" ]; then bibtex main > bib.txt 2>&1; pdflatex -interaction=nonstopmode main.tex > log2.txt 2>&1; pdflatex -interaction=nonstopmode main.tex > log3.txt 2>&1; fi
grep -n '^!' main.log | head -5
python3 -I - <<'PY'
import re
log=open('/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/paperbuild/build_v2/main.log',errors='ignore').read()
cur=None
for m in re.finditer(r'\((\./[\w\-./]+\.(?:tex|tikz))|Overfull \\hbox \(([\d.]+)pt too wide\) (?:in paragraph at lines|detected at line) (\d+)',log):
    if m.group(1):
        if 'figures' not in m.group(1) and 'tikzstyles' not in m.group(1): cur=m.group(1)
    elif float(m.group(2))>2: print(cur,m.group(3),m.group(2))
PY
grep 'Output written' main.log
cp main.pdf /home/user/correlations/nlo_review/paper/jhep_v2/main.pdf
