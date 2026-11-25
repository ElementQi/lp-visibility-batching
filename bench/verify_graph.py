"""Independent rational coordinate optima for graph-derived public LP cases."""

import csv
from fractions import Fraction as F
from math import lcm
from random import Random
import argparse, time
from graph_inputs import mps, svm, graph
from verify import Cube


class RationalCube:
    def __init__(self, forms):
        self.scales = [lcm(*(F(v).denominator for v in row)) for row in forms]
        rows = [[int(x * s) for x in row] for row, s in zip(forms, self.scales)]
        self.lp = Cube(rows)
        for row, scale in zip(self.lp.rows, self.scales):
            row[-1] = scale

    def maximum(self, i):
        a, b = self.lp.maximum(self.lp.forms[i])
        return F(a, b * self.scales[i])


def run(case, data, seed, count, table):
    rows, n = svm(data, count) if count else mps(data)
    rng = Random(seed)
    prices = [rng.randint(-4, 4) if seed else 0 for _ in range(n)]
    m = graph(rows, n, prices)
    r = len(m)
    d = len(m[0])
    eye = lambda k: [[int(i == j) for i in range(k)] for j in range(k)]
    p = RationalCube(eye(r) + list(zip(*m)))
    q = RationalCube([[-x for x in row] for row in m] + eye(d))
    with open(table) as f:
        allrows = [row for row in csv.DictReader(f) if row["case"] == case]
    assert allrows and all(x["status"] != "failure" for x in allrows)
    masks = {row["labels"] for row in allrows}
    assert all(len(s) == n for s in masks)
    visible = 0
    votes = 0
    for i in range(n):
        t = time.monotonic()
        pos = p.maximum(i)
        neg = q.maximum(i)
        assert 0 <= pos <= 1 and 0 <= neg <= 1 and pos * neg == 0
        req = "1" if pos >= F(1, 3) else "0" if neg >= F(1, 3) else "?"
        if req != "?":
            visible += 1
            assert all(mask[i] == req for mask in masks), (case, i, req)
        votes += len(masks)
        print(case, i, req, round(time.monotonic() - t, 3), flush=True)
    print("COMPLETE", case, n, visible, len(masks), votes, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    for key in ("case", "data", "table"):
        p.add_argument(key)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--count", type=int, default=0)
    a = p.parse_args()
    run(a.case, a.data, a.seed, a.count, a.table)
