"""Pricing with actual rational minors and implicit reciprocal scale factors."""

from decimal import Decimal, localcontext
from fractions import Fraction
from math import log2

from .calibration import Minor, _normalized_power, calibrated_chart, price_precision
from .minor_view import MinorView
from .balanced import balanced_chart, balanced_precision
from .lift import realize_graph


def _log_proxy(x, precision):
    """Absolute 2^-precision approximation with a rationally checked fast proposal."""
    if x <= 0:
        raise ValueError("a logarithmic proxy needs a positive argument")
    a, b = x.numerator, x.denominator
    k = a.bit_length() - b.bit_length()
    below_power = a < (b << k) if k >= 0 else (a << -k) < b
    if below_power:
        k -= 1
    # Physical exponents are used only to normalize an existing input ratio.
    target = x / Fraction(2) ** (k + 1)
    count = 1 << (precision + 3)
    lo, hi = -count, 0
    if target == Fraction(1, 2):
        return Fraction(k)
    # Fast approximations only propose a logarithm. Acceptance uses outward
    # rational bounds on both endpoints; bisection handles doubtful proposals.
    if precision <= 45:
        guess = Fraction(round(log2(float(target)) * count), count)
    else:
        with localcontext() as ctx:
            ctx.prec = (precision + 24) * 30103 // 100000 + 8
            arg = Decimal(target.numerator) / Decimal(target.denominator)
            guess = Fraction(round(arg.ln() / Decimal(2).ln() * count), count)
    radius = Fraction(1, 1 << precision)
    power_error = radius / 256
    lower, upper = guess - radius, guess + radius
    if (
        lower <= -1
        or (
            lower <= 0
            and _normalized_power(lower, precision + 8) + power_error <= target
        )
    ) and (
        upper >= 0 or _normalized_power(upper, precision + 8) - power_error >= target
    ):
        return k + 1 + guess
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _normalized_power(Fraction(mid, count), precision + 8) <= target:
            lo = mid
        else:
            hi = mid
    return k + 1 + Fraction(lo, count)


def _echelon(rows, order):
    a = [list(map(Fraction, row)) for row in rows]
    pivots = []
    for j in order:
        k = len(pivots)
        pivot = next((i for i in range(k, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        divisor = a[k][j]
        a[k] = [v / divisor for v in a[k]]
        for i, row in enumerate(a):
            if i != k and row[j]:
                factor = row[j]
                a[i] = [v - factor * w for v, w in zip(row, a[k])]
        pivots.append(j)
    return pivots, a


def _det(rows):
    a = [list(row) for row in rows]
    value = Fraction(1)
    for i in range(len(a)):
        p = next((j for j in range(i, len(a)) if a[j][i]), None)
        if p is None:
            return Fraction(0)
        if p != i:
            a[i], a[p] = a[p], a[i]
            value = -value
        pivot = a[i][i]
        value *= pivot
        for j in range(i + 1, len(a)):
            factor = a[j][i] / pivot
            if factor:
                for t in range(i + 1, len(a)):
                    a[j][t] -= factor * a[i][t]
    return value


class PhysicalSource:
    """A reusable original minor source for explicit rational LP pricing.

    Rows span the selected positive space. A constraint matrix can be supplied
    directly and its kernel selected with polar=True in view(). Precision of
    source operations follows the explicit input's rational height.
    """

    def __init__(self, rows, *, log_precision=None):
        if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
            raise ValueError("expected a nonempty rectangular source")
        self.dimension = len(rows[0])
        pivots, reduced = _echelon(rows, range(self.dimension))
        self.rank = len(pivots)
        self.pivots = tuple(pivots)
        self._index = {col: i for i, col in enumerate(pivots)}
        self._rows = reduced[: self.rank]
        free = [i for i in range(self.dimension) if i not in self._index]
        self._dual = []
        for j in free:
            row = [Fraction(i == j) for i in range(self.dimension)]
            for i, col in enumerate(pivots):
                row[col] = -self._rows[i][j]
            self._dual.append(row)
        self._precision = max(
            price_precision(self.dimension)[0] + 3, log_precision or 0
        )
        self._cache = {}

    def minor(self, cols, *, log_precision=None):
        """Exact determinant sign and approximate log magnitude relative to the first basis."""
        if (
            tuple(sorted(set(cols))) != cols
            or len(cols) != self.rank
            or any(j < 0 or j >= self.dimension for j in cols)
        ):
            raise ValueError("expected sorted distinct parent basis indices")
        precision = max(self._precision, log_precision or 0)
        key = (cols, precision)
        if key in self._cache:
            return self._cache[key]
        selected = [
            (j, self._index[col]) for j, col in enumerate(cols) if col in self._index
        ]
        fixed = {i for _, i in selected}
        free = [col for col in cols if col not in self._index]
        small = [
            [row[col] for col in free]
            for i, row in enumerate(self._rows)
            if i not in fixed
        ]
        value = _det(small)
        if sum(i + j for i, j in selected) % 2:
            value = -value
        result = (
            None
            if not value
            else Minor(1 if value > 0 else -1, _log_proxy(abs(value), precision))
        )
        self._cache[key] = result
        return result

    def _indices(self, pinned):
        pinned = tuple(pinned)
        if len(set(pinned)) != len(pinned) or any(
            i < 0 or i >= self.dimension for i in pinned
        ):
            raise ValueError("invalid pinned original coordinates")
        excluded = set(pinned)
        remaining = tuple(i for i in range(self.dimension) if i not in excluded)
        return pinned, remaining

    def view(self, pinned=(), polar=False, *, log_precision=None):
        """Construct a chart's minor accessor in increasing surviving coordinate order."""
        pinned, remaining = self._indices(pinned)
        excluded = set(pinned)
        index = {i: j for j, i in enumerate(remaining)}
        rows = self._dual if polar else self._rows
        basis, _ = _echelon(rows, (*pinned, *remaining))
        jset = tuple(i for i in basis if i in excluded)
        piv = tuple(index[i] for i in basis if i not in excluded)
        return MinorView(
            self.dimension,
            self.rank,
            polar,
            jset,
            remaining,
            piv,
            lambda cols: self.minor(cols, log_precision=log_precision),
        )

    def chart(self, prices, pinned=(), polar=False):
        view = self.view(pinned, polar)
        if len(prices) != len(view.remaining):
            raise ValueError("expected one price per surviving coordinate")
        return calibrated_chart(prices, len(view.pivots), view)

    def balanced(self, prices, pinned=(), polar=False, *, refresh=True):
        """Choose enough physical-log precision for either balanced construction."""
        n = len(prices)
        bits, _ = balanced_precision(n, n, refresh)
        view = self.view(pinned, polar, log_precision=bits + 3)
        if n != len(view.remaining):
            raise ValueError("expected one price per surviving coordinate")
        return balanced_chart(prices, len(view.pivots), view, refresh=refresh)

    def realize(self, chart, pinned=(), polar=False):
        """Return an exact local diagonal and graph on the source's original rows."""
        pinned, remaining = self._indices(pinned)
        if len(remaining) != len(chart.pivots) + len(chart.free) or set(
            (*chart.pivots, *chart.free)
        ) != set(range(len(remaining))):
            raise ValueError("chart and restricted source have different coordinates")
        requested = tuple(remaining[j] for j in chart.pivots)
        basis, rows = _echelon(
            self._dual if polar else self._rows, (*pinned, *requested, *remaining)
        )
        start = sum(j in pinned for j in basis)
        if tuple(basis[start:]) != requested:
            raise ValueError("chart pivots do not span this restricted source")
        original = [
            [row[remaining[j]] for j in chart.free] for row in rows[start : len(basis)]
        ]
        return realize_graph(chart, original)
