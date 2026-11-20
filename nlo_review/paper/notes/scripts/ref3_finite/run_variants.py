import sys, os, copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib, lines, scan
field = lib.SmoothField(seed=7)
tag = sys.argv[1]
if tag == 'allclusters':
    scan.run_variant('allclusters', lines.U, field, allc=True)
elif tag == 'nlocal':
    us = []
    for u in lines.U:
        v = copy.copy(u); v.extra_merge = [c for c in u.remaining_charges if c != u.keep]
        assert len(v.extra_merge) == 1
        us.append(v)
    scan.run_variant('nlocal', us, field)
elif tag == 'field2':
    scan.run_variant('field2', lines.U, lib.SmoothField(seed=11, nmodes=6, amp=2.0, qscale=2.0))
