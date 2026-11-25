"""Map a restricted primal or polar chart's prefix minors to one parent source."""

from typing import Callable, Optional, Sequence, Tuple

from .calibration import Minor


def _sorting_sign(columns):
    return (
        -1
        if sum(x > y for i, x in enumerate(columns) for y in columns[i + 1 :]) % 2
        else 1
    )


class MinorView:
    """A child chart in the gauge where its chosen pivot columns equal identity.

    Coordinates inside the callback are numbered in remaining order. pinned_basis
    consists of original coordinate indices of an independent basis for the full
    set of pinned columns on the selected side. pivots are local child indices.
    The client discards this view after the compressor reports a restart.
    """

    def __init__(
        self,
        dimension: int,
        rank: int,
        polar: bool,
        pinned_basis: Sequence[int],
        remaining: Sequence[int],
        pivots: Sequence[int],
        parent: Callable[[Tuple[int, ...]], Optional[Minor]],
    ):
        if dimension < 1 or not 0 <= rank <= dimension:
            raise ValueError("invalid parent dimensions")
        self.dimension = dimension
        self.polar = polar
        self.pinned = tuple(pinned_basis)
        self.remaining = tuple(remaining)
        self.pivots = tuple(pivots)
        self.parent = parent
        used = (*self.pinned, *self.remaining)
        if len(set(used)) != len(used) or any(x < 0 or x >= dimension for x in used):
            raise ValueError(
                "pinned and remaining columns must be disjoint parent indices"
            )
        if len(set(pivots)) != len(pivots) or any(
            x < 0 or x >= len(remaining) for x in pivots
        ):
            raise ValueError("invalid child pivot indices")
        if len(self.pinned) + len(self.pivots) != (dimension - rank if polar else rank):
            raise ValueError(
                "pinned basis and child pivots must complete a basis of this side"
            )
        self._base = self._completed(self.pivots)
        if self._base is None:
            raise ValueError("the supplied pivot completion is singular")

    def _completed(self, child_cols):
        cols = (*self.pinned, *(self.remaining[i] for i in child_cols))
        if len(set(cols)) != len(cols):
            return None
        sign = _sorting_sign(cols)
        if self.polar:
            if sum(cols) % 2:
                sign = -sign
            selected = set(cols)
            key = tuple(j for j in range(self.dimension) if j not in selected)
        else:
            key = tuple(sorted(cols))
        value = self.parent(key)
        if value is None:
            return None
        return Minor(sign * value.sign, value.proxy)

    def __call__(self, cols: Tuple[int, ...]) -> Optional[Minor]:
        if len(cols) > len(self.pivots) or any(
            j < 0 or j >= len(self.remaining) for j in cols
        ):
            raise ValueError("invalid prefix minor")
        value = self._completed((*cols, *self.pivots[len(cols) :]))
        if value is None:
            return None
        return Minor(value.sign * self._base.sign, value.proxy - self._base.proxy)
