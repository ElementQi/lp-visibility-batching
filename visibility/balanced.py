"""Deterministic complementary graph calibration by balanced row elimination."""

from fractions import Fraction
from functools import lru_cache

from .calibration import GraphChart, _normalized_power


def _ceil_log2(n):
    return (n - 1).bit_length()


def balanced_precision(n, rank, refresh=True):
    """Absolute working precision and signed-cycle exponent for the source."""
    if n < 1 or not 0 <= rank <= n:
        raise ValueError("expected 0 <= rank <= number of coordinates")
    e = _ceil_log2(10**8 * (n + 1) ** 3)
    if not refresh:
        depth = _ceil_log2(max(1, rank))
        e += depth * _ceil_log2(128 * (n + 1) ** 3)
    return e, e + _ceil_log2(2 * n * (e + 5))


def _stabilize(pivots, free, matrix):
    """Exchange on the short rational matrix. No original minor calls occur."""
    swaps = 0
    while True:
        candidate = next(
            (
                (i, j)
                for i, row in enumerate(matrix)
                for j, value in enumerate(row)
                if abs(value) > 2
            ),
            None,
        )
        if candidate is None:
            return swaps
        i, j = candidate
        pivot = matrix[i][j]
        row = [value / pivot for value in matrix[i]]
        row[j] = 1 / pivot
        for k, other in enumerate(matrix):
            if k == i:
                continue
            factor = other[j]
            if factor:
                for s, value in enumerate(row):
                    if s != j:
                        other[s] -= factor * value
                other[j] = -factor / pivot
        matrix[i] = row
        pivots[i], free[j] = free[j], pivots[i]
        swaps += 1


def balanced_chart(prices, rank, minor, *, refresh=True):
    """Build [I M] using the ordered prefix-minor callback in calibration.py.

    With refresh=True the source precision is O(log N) and the callback count
    O(r N log(r+1)). With refresh=False one cofactor list per leaf suffices,
    at O(log N log(r+1)) precision. Initialize a compressed source with the
    corresponding K from balanced_precision. Restart exceptions propagate.
    """
    n = len(prices)
    bits, _ = balanced_precision(n, rank, refresh)
    if rank in (0, n):
        return GraphChart(
            tuple(range(rank)),
            tuple(range(rank, n)),
            tuple(() for _ in range(rank)),
            Fraction(0),
            0,
            0,
        )
    prices = tuple(map(Fraction, prices))
    lookup = lru_cache(maxsize=None)(minor)
    grid = 1 << bits
    amplification = 128 * (n + 1) ** 3
    swaps = 0

    def rounded(x):
        return Fraction(round(x * grid), grid)

    def cofactor(prefix, candidates):
        scores = {}
        for col in candidates:
            value = lookup(tuple(sorted((*prefix, col))))
            if value is None:
                continue
            if value.sign not in (-1, 1):
                raise ValueError("invalid source determinant sign")
            parity = sum(j > col for j in prefix) % 2
            scores[col] = (
                value.proxy - prices[col],
                -value.sign if parity else value.sign,
            )
        if not scores:
            raise ValueError("vanishing cofactor row in an independent flag")
        top = max(scores, key=lambda col: (scores[col][0], -col))
        return top, {
            j: sign * _normalized_power(score - scores[top][0], bits)
            for j, (score, sign) in scores.items()
        }

    def build(prefix, active, count):
        nonlocal swaps
        if count == 1:
            pivot, row = cofactor(prefix, active)
            free = [j for j in active if j != pivot]
            return (
                [pivot],
                free,
                [[row.get(j, 0) * row[pivot] for j in free]],
                Fraction(9, grid),
            )
        left = count // 2
        a, rest, top, h_top = build(prefix, active, left)
        b, free, bottom, h_bottom = build(prefix + a, rest, count - left)
        pos = {col: j for j, col in enumerate(rest)}
        x = [pos[col] for col in b]
        y = [pos[col] for col in free]
        # Top rows become [I 0 Y-XZ] after subtracting rows [0 I Z].
        matrix = [
            [
                row[col] - sum(row[src] * z[j] for src, z in zip(x, bottom))
                for j, col in enumerate(y)
            ]
            for row in top
        ] + bottom
        basis = a + b
        swaps += _stabilize(basis, free, matrix)
        if refresh:
            for i, pivot in enumerate(basis):
                _, row = cofactor(prefix + basis[:i] + basis[i + 1 :], [pivot, *free])
                if abs(row.get(pivot, 0)) < Fraction(1, 4):
                    raise ValueError(
                        "minor replies violate the common-lift chart bound"
                    )
                matrix[i] = [rounded(row.get(j, 0) / row[pivot]) for j in free]
            error = Fraction(160, grid)
        else:
            matrix = [[rounded(value) for value in row] for row in matrix]
            error = amplification * max(h_top, h_bottom)
        return basis, free, matrix, error

    basis, free, matrix, error = build([], list(range(n)), rank)
    return GraphChart(
        tuple(basis),
        tuple(free),
        tuple(map(tuple, matrix)),
        error,
        swaps,
        lookup.cache_info().misses,
    )
