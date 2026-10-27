"""Equation numbers of your files as LaTeX prints them: (section.n), one number per align row without \\notag/\\nonumber."""
import re
import pe


def numbers(group):
    L = pe.LINES[group]
    sec, n, out = 0, 0, {}
    i = 0
    while i < len(L):
        s = L[i]
        if re.match(r'\s*\\section\{', s):
            sec += 1
            n = 0
        if r'\begin{align}' in s and not s.lstrip().startswith('%'):
            j = i
            while r'\end{align}' not in L[j]:
                j += 1
            body = '\n'.join(L[i:j + 1])
            body = re.sub(r'(?<!\\)%[^\n]*', '', body)
            body = body[body.index(r'\begin{align}') + 13:body.rindex(r'\end{align}')]
            rows = re.split(r'\\\\', body)
            nums = []
            for r in rows:
                if not r.strip():
                    continue
                if re.search(r'\\(notag|nonumber)', r):
                    continue
                n += 1
                nums.append('%d.%d' % (sec, n))
            out[(i + 1, j + 1)] = nums
            i = j
        i += 1
    return out


def eq(group, a, b):
    nums = numbers(group).get((a, b), [])
    return '(' + ', '.join(nums) + ')' if nums else '(unnumbered)'


if __name__ == '__main__':
    for g in ('Group_I', 'Group_II', 'Group_III'):
        for k, v in numbers(g).items():
            print(g, k, v)
