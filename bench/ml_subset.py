"""Prepare a fixed balanced 64-example budget-decision benchmark from scikit-learn."""
from fractions import Fraction
from pathlib import Path
from random import Random
import argparse


def subset(name, destination):
    from sklearn.datasets import load_breast_cancer, load_digits, load_iris, load_wine

    cases = {
        "cancer": (load_breast_cancer, (0, 1)),
        "digits": (load_digits, (4, 9)),
        "iris": (load_iris, (1, 2)),
        "wine": (load_wine, (0, 1)),
    }
    loader, classes = cases[name]
    data, labels = loader(return_X_y=True)
    rng = Random(20261125)
    selected = []
    for label in classes:
        ids = [i for i, y in enumerate(labels) if y == label]
        selected.extend(rng.sample(ids, 32))
    rng.shuffle(selected)
    rows = [[Fraction(str(v)) for v in data[i]] for i in selected]
    limits = [(min(col), max(col)) for col in zip(*rows)]
    with Path(destination).open("x") as out:
        for idx, row in zip(selected, rows):
            entries = []
            for j, (v, (lo, hi)) in enumerate(zip(row, limits), 1):
                w = round(256 * (2 * v - hi - lo) / (hi - lo)) if hi > lo else 0
                if w:
                    entries.append(f"{j}:{w/256:.8f}")
            print(-1 if labels[idx] == classes[0] else 1, *entries, file=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("name", choices=("cancer", "digits", "iris", "wine"))
    p.add_argument("destination")
    a = p.parse_args()
    subset(a.name, a.destination)
