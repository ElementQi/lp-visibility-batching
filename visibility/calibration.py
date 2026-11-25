"""Complementary graph charts for visibility pricing.

A minor callback supplies exact signs and the published compressor's rational
log proxies. Its restart exception propagates to the enclosing price learner.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from typing import Callable, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Minor:
    """Sign of the canonically ordered prefix minor and its log proxy."""

    sign: int
    proxy: Fraction


@dataclass(frozen=True)
class GraphChart:
    """In pivot/free order the positive and polar charts are [I M], [-M^T I]."""

    pivots: Tuple[int, ...]
    free: Tuple[int, ...]
    entries: Tuple[Tuple[Fraction, ...], ...]
    entry_error: Fraction
    exchanges: int
    distinct_requests: int


def price_precision(n: int) -> Tuple[int, int]:
    """Return the dyadic output precision E and the common-lift type exponent K."""
    if n < 1:
        raise ValueError("expected at least one coordinate")
    e = (10**7 * (n + 1) ** 3 - 1).bit_length()
    k = e + (2 * n * (e + 5) - 1).bit_length()
    return e, k


@lru_cache(maxsize=16)
def _roots(precision: int):
    # Scaled lower bounds on 2^(1/2^j). Exact square roots are 1/2-Lipschitz
    # on [1,2], so each lower bound is within two units of the true value.
    count = precision + 4
    scale = 1 << (precision + (64 * count - 1).bit_length())
    value = 2 * scale
    roots = []
    for _ in range(count):
        value = isqrt(value * scale)
        roots.append(value)
    return scale, roots


def _normalized_power(q: Fraction, precision: int) -> Fraction:
    if q > 0:
        raise ValueError("a normalized log ratio must be nonpositive")
    if q < -precision:
        return Fraction(0)
    whole = q.numerator // q.denominator
    scale, roots = _roots(precision)
    count = len(roots)
    bits = int((q - whole) * (1 << count))
    value = scale
    for j, root in enumerate(roots):
        if bits & (1 << (count - 1 - j)):
            value = value * root // scale
    # Binary truncation of q costs <= 2^(1-count) before scaling by 2^whole.
    # Each root and product rounding costs at most 3/scale in relative error.
    # The product is <=2, so rounding once more to nearest grid is safe.
    return Fraction(
        round(Fraction(value << precision, scale << -whole)), 1 << precision
    )


def calibrated_chart(
    prices: Sequence[Fraction],
    rank: int,
    minor: Callable[[Tuple[int, ...]], Optional[Minor]],
) -> GraphChart:
    """Build one LP chart from the positive prefix rows at this price.

    minor(cols) refers to the first len(cols) positive rows and the listed
    sorted columns. None denotes an exactly zero determinant. All proxies must
    belong to the same compression epoch and satisfy the common-lift contract
    with the K returned by price_precision. The cache lasts one query.
    """
    n = len(prices)
    if n < 1 or not 0 <= rank <= n:
        raise ValueError("expected 0 <= rank <= number of prices")
    e, _ = price_precision(n)
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
    cache = {}

    def cofactor(prefix):
        scores = {}
        for j in range(n):
            if j in prefix:
                continue
            columns = (*prefix, j)
            key = tuple(sorted(columns))
            if key not in cache:
                cache[key] = minor(key)
            record = cache[key]
            if record is None:
                continue
            if record.sign not in (-1, 1):
                raise ValueError("a nonzero minor must have sign -1 or +1")
            # The sorting parity internal to the prefix is common to the row.
            # It cancels when normalizing by the pivot sign.
            sign = -record.sign if sum(i > j for i in prefix) % 2 else record.sign
            scores[j] = (Fraction(record.proxy) - prices[j], sign)
        if not scores:
            raise ValueError("a prefix cofactor row vanishes")
        winner = max(scores, key=lambda j: (scores[j][0], -j))
        return scores, winner

    pivots = []
    for _ in range(rank):
        _, j = cofactor(pivots)
        pivots.append(j)

    exchanges = 0
    while True:
        stable = []
        for k, i in enumerate(pivots):
            row, j = cofactor(pivots[:k] + pivots[k + 1 :])
            if i not in row:
                raise ValueError("the previously selected pivot minor vanished")
            if row[j][0] - row[i][0] > 1:
                pivots[k] = j
                exchanges += 1
                if exchanges > 2 * rank * (1 + rank.bit_length()):
                    raise ValueError(
                        "minor proxies exceed the common-lift exchange bound"
                    )
                break
            stable.append((row, j))
        else:
            break

    free = tuple(i for i in range(n) if i not in pivots)
    entries = []
    for i, (row, j) in zip(pivots, stable):
        top = row[j][0]
        pivot = _normalized_power(row[i][0] - top, e)
        entries.append(
            tuple(
                row[k][1] * row[i][1] * _normalized_power(row[k][0] - top, e) / pivot
                if k in row
                else Fraction(0)
                for k in free
            )
        )
    return GraphChart(
        tuple(pivots), free, tuple(entries), Fraction(96, 1 << e), exchanges, len(cache)
    )
