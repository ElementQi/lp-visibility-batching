"""Compile all one-coordinate price changes from one epoch's stable minor table."""

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .calibration import (
    GraphChart,
    _normalized_power,
    calibrated_chart,
    price_precision,
)


def _value(q, bits):
    if q <= 0:
        return _normalized_power(q, bits)
    if q > 1:
        raise ValueError("a stable chart coefficient must have log at most one")
    return 2 * _normalized_power(q - 1, bits + 1)


def _order_sign(cols):
    return (
        -1 if sum(a > b for i, a in enumerate(cols) for b in cols[i + 1 :]) % 2 else 1
    )


@dataclass(frozen=True)
class PriceLine:
    """Two charts and one bounded scalar evaluation at each point of the line."""

    center: Fraction
    peak: Fraction
    positive_direction: bool
    pieces: tuple
    precision: int
    constant: bool

    def chart(self, price):
        if self.constant:
            return self.pieces[0][0]
        offset = Fraction(price) - self.center
        u = self.peak + (offset if self.positive_direction else -offset)
        exchanged = u > 1
        template, along_row, position = self.pieces[exchanged]
        factor = (
            _normalized_power(-u, self.precision)
            if exchanged
            else _value(u, self.precision)
        )
        entries = tuple(
            tuple(
                x * factor if (i if along_row else j) == position else x
                for j, x in enumerate(row)
            )
            for i, row in enumerate(template.entries)
        )
        return GraphChart(
            template.pivots,
            template.free,
            entries,
            template.entry_error,
            int(exchanged),
            0,
        )


class AxisTable:
    """One-time source access for every axis through the supplied center price."""

    def __init__(self, prices, rank, minor):
        self.prices = tuple(map(Fraction, prices))
        cached = lru_cache(maxsize=None)(minor)
        self.base = calibrated_chart(self.prices, rank, cached)
        self._e = price_precision(len(prices))[0]
        p = self.base.pivots
        self._scores = []
        if self.base.free and p:
            base = cached(tuple(sorted(p)))
            if base is None:
                raise ValueError("the stable pivot minor vanished")
            for i, col in enumerate(p):
                row = []
                for j in self.base.free:
                    cols = (*p[:i], j, *p[i + 1 :])
                    value = cached(tuple(sorted(cols)))
                    if value is None:
                        row.append(None)
                    else:
                        sign = (
                            _order_sign(cols) * _order_sign(p) * value.sign * base.sign
                        )
                        score = (
                            value.proxy - base.proxy + self.prices[col] - self.prices[j]
                        )
                        if score > 1:
                            raise ValueError(
                                "minor table is inconsistent with the stable chart"
                            )
                        row.append((score, sign))
                self._scores.append(row)
        self.distinct_requests = cached.cache_info().misses

    def line(self, coordinate):
        if not 0 <= coordinate < len(self.prices):
            raise ValueError("invalid coordinate")
        piv = list(self.base.pivots)
        free = list(self.base.free)
        is_row = coordinate in piv
        idx = piv.index(coordinate) if is_row else free.index(coordinate)
        locations = (
            [(idx, j) for j in range(len(free))]
            if is_row
            else [(i, idx) for i in range(len(piv))]
        )
        active = [(i, j) for i, j in locations if self._scores[i][j] is not None]
        if not active:
            constant = GraphChart(
                tuple(piv), tuple(free), self.base.entries, self.base.entry_error, 0, 0
            )
            return PriceLine(
                self.prices[coordinate],
                Fraction(0),
                is_row,
                ((constant, is_row, idx),),
                self._e,
                True,
            )
        a, b = max(active, key=lambda ij: self._scores[ij[0]][ij[1]][0])
        peak, sign = self._scores[a][b]
        # The neutral template uses a signed unit at its largest affected entry.
        neutral = [
            [
                s[1] * _value(s[0], self._e) if s is not None else Fraction(0)
                for s in row
            ]
            for row in self._scores
        ]
        for i, j in active:
            t, h = self._scores[i][j]
            neutral[i][j] = h * _normalized_power(t - peak, self._e)
        error = Fraction(512, 1 << self._e)
        ordinary = GraphChart(
            tuple(piv), tuple(free), tuple(map(tuple, neutral)), error, 0, 0
        )
        swapped = [row[:] for row in neutral]
        for i, row in enumerate(neutral):
            for j, val in enumerate(row):
                if i == a and j == b:
                    swapped[i][j] = Fraction(sign)
                elif i == a:
                    swapped[i][j] = sign * val
                elif j == b:
                    swapped[i][j] = -sign * val
                else:
                    swapped[i][j] = val - sign * neutral[i][b] * neutral[a][j]
        piv[a], free[b] = free[b], piv[a]
        alternative = GraphChart(
            tuple(piv), tuple(free), tuple(map(tuple, swapped)), error, 1, 0
        )
        # After pivoting the inverse scalar appears in the orthogonal axis.
        return PriceLine(
            self.prices[coordinate],
            peak,
            is_row,
            ((ordinary, is_row, idx), (alternative, not is_row, b if is_row else a)),
            self._e,
            False,
        )
