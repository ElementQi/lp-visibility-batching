"""Recover a rational positive column scaling consistent with a returned graph."""

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class GraphLift:
    scales: tuple
    pivots: tuple
    free: tuple
    entries: tuple

    def original_vector(self, coefficients, *, polar=False):
        """Evaluate a graph row combination and undo its positive or polar scaling."""
        coefficients = tuple(map(Fraction, coefficients))
        piv, free = (self.free, self.pivots) if polar else (self.pivots, self.free)
        if len(coefficients) != len(piv):
            raise ValueError("one coefficient per pivot is required")
        values = [Fraction(0)] * len(self.scales)
        for i, a in zip(piv, coefficients):
            values[i] = a
        for j, col in enumerate(free):
            if polar:
                values[col] = -sum(a * v for a, v in zip(coefficients, self.entries[j]))
            else:
                values[col] = sum(
                    a * row[j] for a, row in zip(coefficients, self.entries)
                )
        return tuple(v * s if polar else v / s for v, s in zip(values, self.scales))


def realize_graph(chart, original):
    """Find exact scales for an unscaled graph block in the chart's pivot/free order.

    A supplied entry uses its signed interval. A real positive solution implies
    a rational solution with all coordinates products along simple paths. A
    ValueError identifies an inconsistent graph or incompatible entry bounds.
    """
    piv, free = chart.pivots, chart.free
    n = len(piv) + len(free)
    if n == 0 or set((*piv, *free)) != set(range(n)):
        raise ValueError("expected a partition of local coordinate indices")
    if len(original) != len(piv) or len(chart.entries) != len(piv):
        raise ValueError("row count differs from the supplied pivot set")
    error = Fraction(chart.entry_error)
    if error < 0:
        raise ValueError("negative error bound")
    edges = []
    source = []
    for i, row, approx in zip(piv, original, chart.entries):
        row = tuple(map(Fraction, row))
        source.append(row)
        if len(row) != len(free) or len(approx) != len(free):
            raise ValueError("column count differs from the supplied free set")
        for j, v, w in zip(free, row, approx):
            w = Fraction(w)
            if not v:
                if abs(w) > error:
                    raise ValueError("an exactly zero minor is outside its interval")
                continue
            oriented = w if v > 0 else -w
            lower = max(Fraction(0), oriented - error)
            upper = oriented + error
            if upper <= 0:
                raise ValueError("a nonzero minor has no compatible signed magnitude")
            edges.append((i, j, upper / abs(v)))
            if lower:
                edges.append((j, i, abs(v) / lower))
    # The all-ones vector represents a supersource with cost one to each vertex.
    # Synchronous updates cap every intermediate product length even on bad input.
    q = [Fraction(1)] * n
    for _ in range(n - 1):
        next_q = q.copy()
        for i, j, c in edges:
            next_q[j] = min(next_q[j], q[i] * c)
        if next_q == q:
            break
        q = next_q
    if any(q[j] > q[i] * c for i, j, c in edges):
        raise ValueError("incompatible multiplicative cycle")
    scaled = tuple(
        tuple(v * q[j] / q[i] for j, v in zip(free, row)) for i, row in zip(piv, source)
    )
    return GraphLift(tuple(q), piv, free, scaled)
