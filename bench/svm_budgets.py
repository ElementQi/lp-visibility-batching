"""Rational certificates on both sides of the actual L1 soft-margin optimum."""

from fractions import Fraction as F
from pathlib import Path
from math import ceil, floor
from concurrent.futures import ProcessPoolExecutor
from graph_inputs import svm, graph
import argparse


def prepare(job, tag_prefix="heart"):
    m, source, solution_dir, out = job
    model, alpha = [
        list(map(F, row.split()))
        for row in (Path(solution_dir) / f"visibility_svm_solution_{m}")
        .read_text()
        .splitlines()
    ]
    rows, n = svm(source, m)
    d = len(model) - 1
    train = []
    for i, line in enumerate(Path(source).read_text().splitlines()[:m]):
        y, *features = line.split()
        v = [F(0)] * d
        for f in features:
            j, x = f.split(":")
            v[int(j) - 1] = F(x)
        train.append((int(y), v))
    grid = 1 << 30
    alpha = [max(F(0), min(F(1), F(round(a * grid), grid))) for a in alpha]
    mass = {s: sum(a for a, (y, _) in zip(alpha, train) if y == s) for s in (-1, 1)}
    assert min(mass.values()) > 0
    alpha = [a * min(mass.values()) / mass[y] for a, (y, _) in zip(alpha, train)]
    correlations = [
        sum(a * y * v[j] for a, (y, v) in zip(alpha, train)) for j in range(d)
    ]
    alpha = [a / max([F(1), *map(abs, correlations)]) for a in alpha]
    assert (
        all(0 <= a <= 1 for a in alpha)
        and sum(a * y for a, (y, _) in zip(alpha, train)) == 0
    )
    for j in range(d):
        assert abs(sum(a * y * v[j] for a, (y, v) in zip(alpha, train))) <= 1
    lower = sum(alpha)
    model = [F(round(x * grid), grid) for x in model]

    def primal(w):
        margin = [y * (w[-1] + sum(x * t for x, t in zip(v, w))) for y, v in train]
        slack = [max(0, 1 - h) for h in margin]
        upper = sum(slack) + sum(map(abs, w[:-1]))
        return upper, slack, [h + s - 1 for h, s in zip(margin, slack)]

    if primal([-x for x in model])[0] < primal(model)[0]:
        model = [-x for x in model]
    upper, slack, surplus = primal(model)
    assert 0 <= upper - lower < F(1, 10000)
    for name, budget in [
        ("below", F(floor(990 * lower), 1000)),
        ("above", F(ceil(1010 * upper), 1000)),
    ]:
        rows[-1][n - 1] = -budget
        if name == "above":
            witness = (
                [x for t in model for x in (max(t, 0), max(-t, 0))]
                + slack
                + surplus
                + [budget - upper, 1]
            )
            assert min(witness) >= 0 and all(
                sum(witness[j] * v for j, v in row.items()) == 0 for row in rows
            )
        else:
            coeff = [-a for a in alpha] + [1]
            witness = [
                sum(c * row.get(j, 0) for c, row in zip(coeff, rows)) for j in range(n)
            ]
            assert min(witness) >= 0 and witness[-1] == lower - budget > 0
        # Fixed radius declared from sample count, with both price signs tested.
        radius = 2 * (m + 1).bit_length() + 4
        for price in (-radius, 0, radius):
            prices = [0] * (n - 1) + [price]
            matrix, order = graph(rows, n, prices, return_order=True)
            scale = F(2) ** (price if name == "below" else -price)
            head = witness[-1] * scale
            ratio = head / max(head, max(witness[:-1]))
            tag = f"{tag_prefix}{m}_{name}_{price}"
            with (Path(out) / tag).open("w") as f:
                print(tag, len(matrix), len(matrix[0]), file=f)
                for row in matrix:
                    print(*(float(x).hex() for x in row), file=f)
            print(
                tag,
                n,
                order.index(n - 1),
                str(budget),
                float(lower),
                float(upper),
                float(ratio),
                int(ratio >= F(1, 3)),
                int(name == "above"),
                sep=",",
                flush=True,
            )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--tag-prefix", default="heart")
    p.add_argument("source")
    p.add_argument("solution_dir")
    p.add_argument("out")
    p.add_argument("sizes", nargs="+", type=int)
    a = p.parse_args()
    Path(a.out).mkdir(exist_ok=True)
    print(
        "case,n,target,budget,dual_lower,primal_upper,witness_head_fraction,mandatory,true_side",
        flush=True,
    )
    from functools import partial

    with ProcessPoolExecutor(4) as pool:
        list(
            pool.map(
                partial(prepare, tag_prefix=a.tag_prefix),
                [(m, a.source, a.solution_dir, a.out) for m in a.sizes],
            )
        )
