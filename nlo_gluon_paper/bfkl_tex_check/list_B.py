import os, sys, re
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'full_check'))
import ll  # noqa
import three_check as T3
import texeval as TE
c0 = T3.configs(1)
for lab, t, cc in T3.row_terms('B', True) + T3.hidden_B():
    try:
        d = T3.degree(t, c0)
    except Exception as e:
        d = 'ERR ' + str(e)[:80]
    code = TE.translate(t)
    ph = set(re.findall(r'KD\("([^"]*)","([^"]*)"\)', code))
    ch = sorted(set(re.findall(r'CR\("[^"]*","([^"]*)"\)', code)))
    print(lab, d, cc, ph, ch)
