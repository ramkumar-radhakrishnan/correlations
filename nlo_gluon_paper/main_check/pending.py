r"""Which corrections of ../groups23/main/Main_corrections.pdf (E1-E11, X1-X7, G1-G4, M6) are already in the new
Calculation_main.tex (input/) and which are still needed.  Runs every section of ../groups23/main/make_corrected_main.py
on the new file: a section whose 'as written' text is still found is still needed; one whose text is no longer found
was changed by you (then the corrected text is looked for).  Output: pending.json"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, '..', 'groups23', 'main', 'make_corrected_main.py')).read()
head = SRC[:SRC.index('# ------------------------------------------------------------------ header')]
head = head.replace("src = open(os.path.join(HERE, '..', 'input', 'Calculation_main.tex')).read()",
                    "src = open(NEWFILE).read()")
body = SRC[SRC.index('# ------------------------------------------------------------------ X1'):SRC.index("open(os.path.join(HERE, 'Calculation_main_corrected.tex'), 'w')")]
parts = re.split(r'(# -{60,} [^\n]*)\n', body)
ns = {'NEWFILE': os.path.join(HERE, 'input', 'Calculation_main.tex'), '__file__': os.path.join(HERE, '..', 'groups23', 'main', 'x.py')}
exec(head, ns)
out = {}
pre = ''
for k in range(1, len(parts), 2):
    label = parts[k].split(' ', 2)[2].strip() if False else parts[k].lstrip('# -').strip()
    code = parts[k + 1]
    if code.lstrip().startswith('def ') and 'fix_line' in code[:200]:
        pass
    s_before = ns['s']
    try:
        exec(code, ns)
        out[label] = 'still needed (the text as written is unchanged)'
    except AssertionError as e:
        ns['s'] = s_before
        out[label] = 'NOT matched: the text was changed (%s)' % str(e)[:160]
    print('%-55s %s' % (label, out[label][:120]))
json.dump(out, open(os.path.join(HERE, 'pending.json'), 'w'), indent=1)
