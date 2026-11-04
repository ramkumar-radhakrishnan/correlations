r"""Every step of your new Group_I/II/III.tex (input/), with the step list of ../rows_check/check_all.py.
Usage: python3 step_check.py Group_I|Group_II|Group_III  -> step_<group>.json"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'rows_check')))
import pe  # noqa: E402

pe.LINES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n') for g in ('Group_I', 'Group_II', 'Group_III')}
import check_all as CA  # noqa: E402

if __name__ == '__main__':
    group = sys.argv[1]
    out = {}
    for name, kind, e1, e2, opt in CA.STEPS[group]:
        res = [CA.run_step(group, kind, e1, e2, opt, seed) for seed in (1, 2)]
        out[name] = res
        print('%-48s %s' % (name, res), flush=True)
    if group == 'Group_I':
        out['Row V 2b: UV term'] = [CA.extra_group_I(s) for s in (1, 2)]
        print(out['Row V 2b: UV term'], flush=True)
    out['unmatched braces'] = sorted('%s %d-%d' % u for u in CA.UNMATCHED)
    json.dump(out, open(os.path.join(HERE, 'step_%s.json' % group), 'w'), indent=1)
