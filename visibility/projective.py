"""Positive spanning-forest normalization in a fixed column order."""

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product

from .lift import GraphLift


@dataclass(frozen=True)
class ForestForm(GraphLift):
    # Indices (row, column) in the free block, in row-major order.
    chords: tuple

    def rows(self):
        """Return the normalized row matrix in original coordinate order."""
        matrix = []
        for pivot, entries in zip(self.pivots, self.entries):
            row = [Fraction(0)] * len(self.scales)
            row[pivot] = Fraction(1)
            for j, v in zip(self.free, entries):
                row[j] = v
            matrix.append(row)
        return matrix


def forest_form(source):
    """Remove all positive diagonal gauges from an explicit PhysicalSource.

    The signed forest is fixed by the support and column order. Only chord
    magnitudes can differ across inequivalent inputs with that signed support.
    """
    n = source.dimension
    pivots = source.pivots
    fixed = set(pivots)
    free = tuple(j for j in range(n) if j not in fixed)
    edges = [[] for _ in range(n)]
    nonzeros = []
    for i, (p, row) in enumerate(zip(pivots, source._rows)):
        for k, j in enumerate(free):
            v = row[j]
            if v:
                # q_j/q_p = 1/|M_pj| on tree edges.
                edges[p].append((j, 1 / abs(v), (i, k)))
                edges[j].append((p, abs(v), (i, k)))
                nonzeros.append((i, k))
    scales = [None] * n
    tree = set()
    for root in range(n):
        if scales[root] is not None:
            continue
        scales[root] = Fraction(1)
        pending = deque([root])
        while pending:
            i = pending.popleft()
            for j, factor, edge in edges[i]:
                if scales[j] is None:
                    scales[j] = scales[i] * factor
                    tree.add(edge)
                    pending.append(j)
    matrix = tuple(
        tuple(row[j] * scales[j] / scales[i] for j in free)
        for row, i in zip(source._rows, pivots)
    )
    chords = tuple(e for e in nonzeros if e not in tree)
    return ForestForm(tuple(scales), pivots, free, matrix, chords)


def _profile(x, y):
    # Clear positive denominators once before the 40 comparisons.
    x, y = Fraction(x), Fraction(y)
    m, n = x.numerator, x.denominator
    p, q = y.numerator, y.denominator
    out = []
    for a, b, c, d in product((-1, 0, 1), repeat=4):
        if next((z for z in (a, b, c, d) if z), 0) != 1:
            continue
        v = (a * n + b * m) * q + (c * n + d * m) * p
        out.append((v > 0) - (v < 0))
    return tuple(out)


@lru_cache(maxsize=1)
def _two_chord_table():
    table = {}
    grid = sorted({Fraction(p, q) for p in range(1, 9) for q in range(1, 9)})
    for x in grid:
        for y in grid:
            table.setdefault(_profile(x, y), (x, y))
    return table


def small_chord_rows(form):
    """A four-bit row representation with every maximal-minor sign preserved.

    The forest form must have at most two chords. The returned row matrix
    represents the same oriented matroid in original coordinate order.
    """
    chords = form.chords
    if len(chords) > 2:
        raise ValueError("the finite table covers at most two chords")
    data = [form.entries[i][j] for i, j in chords]
    if len(chords) == 1:
        x = abs(data[0])
        values = (Fraction(1, 2) if x < 1 else Fraction(2) if x > 1 else Fraction(1),)
    elif len(chords) == 2:
        values = _two_chord_table()[_profile(abs(data[0]), abs(data[1]))]
    else:
        values = ()
    rows = form.rows()
    for (i, j), original, value in zip(chords, data, values):
        rows[i][form.free[j]] = value if original > 0 else -value
    return rows


def lift_elementary_target(original_rows, short_rows, witness, target):
    """Recover an original nonnegative row from a short elementary target support.

    Both row matrices must span the selected side in common coordinate order.
    The exact tests also allow callers to try general approximate charts.
    Returns None if the support fails to specify a compatible original row.
    """
    from .physical import _echelon

    n = len(witness)
    if not 0 <= target < n or min(witness) < 0 or witness[target] <= 0:
        raise ValueError("expected a nonnegative witness with positive target")
    if any(len(row) != n for row in (*original_rows, *short_rows)):
        raise ValueError("coordinate count mismatch")
    zeros = [j for j, x in enumerate(witness) if not x]
    base, _ = _echelon(short_rows, zeros)
    rank = len(_echelon(short_rows, range(n))[0])
    if len(base) + 1 != rank:
        return None
    piv, a = _echelon(original_rows, [*base, target, *range(n)])
    if piv != [*base, target]:
        return None
    row = a[len(base)]
    return tuple(row) if min(row) >= 0 else None
