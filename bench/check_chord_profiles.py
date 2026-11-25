"""Exact small representatives for all rational two-chord determinant sign profiles."""

from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from visibility.projective import _two_chord_table, _profile
from math import isqrt, gcd


# Negation supplies the complementary signs. Keep one orientation per polynomial.
FORMS = [v for v in product((-1, 0, 1), repeat=4) if next((c for c in v if c), 0) == 1]


def polynomials():
    out = {(0, 1, 0)}
    for a, b, c, d in FORMS:
        out.update(((a, b, 0), (c, d, 0)))
    for (a, b, c, d), (e, f, g, h) in combinations(FORMS, 2):
        out.add((a * g - e * c, a * h + b * g - e * d - f * c, b * h - f * d))
    return out


def brackets(coeff, bits):
    a, b, c = coeff
    if not c:
        return [(F(-a, b),) * 2] if b else []
    if c < 0:
        return brackets((-a, -b, -c), bits)
    disc = b * b - 4 * a * c
    if disc < 0:
        return []
    unit = 1 << bits
    v = isqrt(disc * unit * unit)
    hi = v + (v * v != disc * unit * unit)
    return (
        [(F(-b * unit - v, 2 * c * unit), F(-b * unit - v, 2 * c * unit))]
        if disc == 0
        else [
            (F(-b * unit - hi, 2 * c * unit), F(-b * unit - v, 2 * c * unit)),
            (F(-b * unit + v, 2 * c * unit), F(-b * unit + hi, 2 * c * unit)),
        ]
    )


def horizontal_samples():
    polys = set()
    for w in polynomials():
        g = gcd(*w)
        if not g:
            continue
        z = tuple(v // g for v in w)
        if next(v for v in reversed(z) if v) < 0:
            z = tuple(-v for v in z)
        polys.add(z)
    # Distinct irrational roots have distinct primitive irreducible quadratics.
    # A repeated rational root has the same singleton bracket after reduction.
    bits = 4
    while True:
        intervals = sorted(
            set(x for p in polys for x in brackets(p, bits) if x[1] >= 0)
        )
        if all(a[1] < b[0] for a, b in zip(intervals, intervals[1:])):
            break
        bits *= 2
    samples = {b for a, b in intervals if a == b and b > 0}
    for (_, a), (b, _) in zip(intervals, intervals[1:]):
        samples.add((a + b) / 2)
    samples.add(intervals[-1][1] + 1)
    return sorted(samples)


def key(x, y):
    x, y = F(x), F(y)
    m, n = x.numerator, x.denominator
    p, q = y.numerator, y.denominator
    vals = ((a * n + b * m) * q + (c * n + d * m) * p for a, b, c, d in FORMS)
    return tuple((v > 0) - (v < 0) for v in vals)


@lru_cache(maxsize=1)
def profiles():
    table = {}
    for x in horizontal_samples():
        cuts = sorted(
            {F(0)}
            | {
                -(a + b * x) / (c + d * x)
                for a, b, c, d in FORMS
                if c + d * x and -(a + b * x) / (c + d * x) > 0
            }
        )
        samples = (
            cuts[1:] + [(a + b) / 2 for a, b in zip(cuts, cuts[1:])] + [cuts[-1] + 1]
        )
        for y in samples:
            k = key(x, y)
            size = lambda z: max(
                max(a.numerator.bit_length(), a.denominator.bit_length()) for a in z
            )
            if k not in table or size((x, y)) < size(table[k]):
                table[k] = (x, y)
    # The exact arrangement above supplies the completeness certificate.
    # Pick the smallest square numerator/denominator cap covering that table.
    simple = {}
    for height in range(1, 32):
        grid = sorted(
            {F(i, j) for i in range(1, height + 1) for j in range(1, height + 1)}
        )
        for x in grid:
            for y in grid:
                k = key(x, y)
                assert k in table
                simple.setdefault(k, (x, y))
        if len(simple) == len(table):
            return simple
    return table


def representative(values):
    if len(values) > 2 or any(x <= 0 for x in values):
        raise ValueError("at most two positive chord magnitudes")
    if not values:
        return ()
    if len(values) == 1:
        x = values[0]
        return (F(1, 2) if x < 1 else F(2) if x > 1 else F(1),)
    return profiles()[key(*values)]


if __name__ == "__main__":
    from random import Random

    rng = Random(44731)
    points = [(F(1, 2**p), F(2**p + 1, 2**p)) for p in (10, 50, 100, 1000)]
    for _ in range(500):
        x = F(rng.randint(1, 500), rng.randint(1, 100))
        cands = [
            -(a + b * x) / (c + d * x)
            for a, b, c, d in FORMS
            if c + d * x and -(a + b * x) / (c + d * x) > 0
        ]
        points.append((x, rng.choice(cands)))
        points.append((x, F(rng.randint(1, 500), rng.randint(1, 150))))
    for x, y in points:
        assert key(x, y) == key(*representative([x, y]))
    table = profiles()
    published = _two_chord_table()
    assert set(published) == set(table) and all(
        _profile(x, y) == k for k, (x, y) in published.items()
    )
    assert len(table) == 201
    assert (
        max(max(v.numerator, v.denominator) for pair in table.values() for v in pair)
        == 8
    )
    print(
        "forms",
        len(FORMS),
        "x samples",
        len(horizontal_samples()),
        "rational sign profiles",
        len(table),
    )
    print(
        "maximum representative numerator/denominator bits",
        max(
            max(v.numerator.bit_length(), v.denominator.bit_length())
            for pair in table.values()
            for v in pair
        ),
    )
    print("independent exact evaluations", len(points))
