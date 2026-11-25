"""Audit every screened and residual label against exact normalized cone widths."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median
import sys

from verify import Cube


def widths(n, data):
    r = n // 2
    assert len(data) == r*r
    matrix = [data[i:i+r] for i in range(0, len(data), r)]
    eye = [[int(i == j) for j in range(r)] for i in range(r)]
    primal = Cube(eye + list(zip(*matrix)))
    dual = Cube([[-a for a in row] for row in matrix] + eye)
    signs = []
    for i in range(n):
        a, b = primal.maximum(primal.forms[i])
        c, d = dual.maximum(dual.forms[i])
        assert a*c == 0
        signs.append('1' if 3*a >= b else '0' if 3*c >= d else '?')
    return ''.join(signs)


def audit(matrices, timings):
    truth, keys = {}, {}
    for line in Path(matrices).read_text().splitlines():
        n, family, seed, *fields = line.split()
        data = tuple(map(int, fields))
        if data not in truth:
            truth[data] = widths(int(n), data)
        keys[n, family, seed] = truth[data]
    seen, summaries = set(), defaultdict(list)
    with open(timings) as stream:
        for row in csv.DictReader(stream):
            n, family, seed, gauge, rep = (row[k] for k in ('N','family','seed','gauge','rep'))
            key = (n, family, seed, gauge, int(rep))
            assert key not in seen
            seen.add(key)
            expected = keys[n, family, seed]
            screen, batch, result = (row[k] for k in ('screen','batch_mask','screened_mask'))
            assert all(len(m) == len(expected) for m in (screen, batch, result))
            assert set(batch+result) <= set('01') and set(screen) <= set('01?')
            for correct, cert, old, new in zip(expected, screen, batch, result):
                assert cert == '?' or cert == new, (key, cert, new)
                assert correct == '?' or correct == old == new, (key, correct, old, new)
            summaries[int(n), family, gauge].append(
                (float(row['batch_ms']),float(row['screened_ms']),len(screen)-screen.count('?')))
    assert seen == {(*key, g, rep) for key in keys for g in ('shared','varied') for rep in range(5)}
    print(f'{len(truth)} distinct matrices, {sum(len(v)-v.count("?") for v in truth.values())} required signs, {len(seen)} audited observations')
    print('N family gauge batch_ms screened_ms ratio resolved_min resolved_max')
    for (n,f,g), records in sorted(summaries.items()):
        old,new,counts = zip(*records)
        a,b = median(old),median(new)
        print(n,f,g,f'{a:.3f}',f'{b:.3f}',f'{a/b:.2f}',min(counts),max(counts))


if __name__ == '__main__':
    audit(*sys.argv[1:])
