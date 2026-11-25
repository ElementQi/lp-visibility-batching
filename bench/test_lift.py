"""Exact positive diagonal realization with clipping, disconnected blocks and bad cycles."""

from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys, unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from visibility.calibration import GraphChart
from visibility.lift import realize_graph


def fixture(piv, free, b, error):
    return GraphChart(tuple(piv), tuple(free), tuple(map(tuple, b)), F(error), 0, 0)


class LiftTest(unittest.TestCase):
    def test_extreme_scales_and_both_sides(self):
        rng = Random(88103)
        for seed in range(30):
            n = 8
            order = rng.sample(range(n), n)
            piv, free = order[:3], order[3:]
            true = [[F(rng.randint(-9, 9), 8) for _ in free] for _ in piv]
            true[0][0] = F(1, 1 << 150)
            true[1][1] = 0
            hidden = [F(2) ** rng.randrange(-400, 401) for _ in range(n)]
            original = [
                [v * hidden[i] / hidden[j] for v, j in zip(row, free)]
                for row, i in zip(true, piv)
            ]
            coarse = [[F(round(v * (1 << 16)), 1 << 16) for v in row] for row in true]
            result = realize_graph(fixture(piv, free, coarse, F(1, 1 << 16)), original)
            for a, b in zip(result.entries, coarse):
                self.assertTrue(all(abs(x - y) <= F(1, 1 << 16) for x, y in zip(a, b)))
            x = result.original_vector([F(3, 5), -2, 1])
            y = result.original_vector([2, -3, 1, F(1, 11), 0], polar=True)
            self.assertEqual(sum(a * b for a, b in zip(x, y)), 0)
            for j, c in enumerate(free):
                self.assertEqual(
                    x[c], sum(x[i] * row[j] for i, row in zip(piv, original))
                )

    def test_ties_zeros_and_negative_cycle(self):
        exact = fixture([0, 1], [2, 3], [[1, 1], [1, 2]], 0)
        with self.assertRaisesRegex(ValueError, "cycle"):
            realize_graph(exact, [[1, 1], [1, 1]])
        with self.assertRaisesRegex(ValueError, "zero minor"):
            realize_graph(exact, [[0, 1], [1, 2]])
        with self.assertRaisesRegex(ValueError, "signed magnitude"):
            realize_graph(exact, [[-1, 1], [1, 2]])
        lift = realize_graph(fixture([0, 1], [2], [[0], [0]], 0), [[0], [0]])
        self.assertEqual(lift.scales, (1, 1, 1))
        self.assertEqual(lift.original_vector([2, 3]), (2, 3, 0))

    def test_source_views_and_original_equations(self):
        from visibility.physical import PhysicalSource

        a = [[1, 0, 2, -3, 1], [0, 1, 5, 4, -1]]
        source = PhysicalSource(a)
        for pinned in ((), (1,), (0, 2)):
            live = [j for j in range(5) if j not in pinned]
            for polar in (False, True):
                for refresh in (False, True):
                    chart = source.balanced(
                        [10**99] + [0] * (len(live) - 1), pinned, polar, refresh=refresh
                    )
                    lift = source.realize(chart, pinned, polar)
                    x = [0] * 5
                    values = lift.original_vector(range(1, len(chart.pivots) + 1))
                    for j, v in zip(live, values):
                        x[j] = v
                    opposite = (
                        a
                        if polar
                        else [[-2, -5, 1, 0, 0], [3, -4, 0, 1, 0], [-1, 1, 0, 0, 1]]
                    )
                    for row in opposite:
                        self.assertEqual(sum(v * w for v, w in zip(x, row)), 0)
        with self.assertRaisesRegex(ValueError, "coordinates"):
            source.realize(source.balanced([0] * 5), (2,))
