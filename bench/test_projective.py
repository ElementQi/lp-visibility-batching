"""Forest normal forms under arbitrary row and positive column gauges."""

from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from visibility.physical import PhysicalSource
from visibility.projective import forest_form


class ProjectiveTest(unittest.TestCase):
    def test_scaled_equivalence_and_records(self):
        rng = Random(277991)
        for seed in range(60):
            r, n = 1 + seed % 5, 7 + seed % 4
            rows = [
                [
                    F(int(i == j))
                    if j < r
                    else F(rng.choice([0, 0, -3, 2, 4]), rng.randint(1, 7))
                    for j in range(n)
                ]
                for i in range(r)
            ]
            prices = [
                F(2) ** rng.randint(-150, 150)
                * F(rng.randint(1, 15), rng.randint(1, 15))
                for _ in range(n)
            ]
            shifted = [[x * d for x, d in zip(row, prices)] for row in rows]
            for i in range(1, r):
                shifted[i] = [x + 3 * y for x, y in zip(shifted[i], shifted[i - 1])]
            shifts = PhysicalSource(shifted[::-1])
            a, b = forest_form(PhysicalSource(rows)), forest_form(shifts)
            self.assertEqual((a.rows(), a.chords), (b.rows(), b.chords))
            for i, row in enumerate(a.entries):
                for j, x in enumerate(row):
                    if (i, j) not in a.chords:
                        self.assertIn(x, (-1, 0, 1))
            for side, count, equations in (
                (False, r, shifts._dual),
                (True, n - r, shifts._rows),
            ):
                v = b.original_vector(range(1, count + 1), polar=side)
                self.assertTrue(
                    all(sum(x * y for x, y in zip(row, v)) == 0 for row in equations)
                )

    def test_loops_and_coloops(self):
        for a in ([[0, 0, 0]], [[1, 0], [0, 1]], [[1, 0, 0, 1], [0, 1, 0, 0]]):
            chart = forest_form(PhysicalSource(a))
            self.assertFalse(chart.chords)
            self.assertTrue(all(x > 0 for x in chart.scales))

    def test_short_chords_and_original_circuit(self):
        from itertools import combinations
        from visibility.physical import _det
        from visibility.projective import small_chord_rows, lift_elementary_target

        for x, y in [(F(1, 2**300), F(2**400)), (F(3, 5), F(8, 5)), (F(1), F(1))]:
            a = [[1, 0, 2, 4, 3], [0, 1, 1, x, y]]
            source = PhysicalSource(a)
            form = forest_form(source)
            self.assertEqual(len(form.chords), 2)
            b = small_chord_rows(form)

            def sign(v):
                return (v > 0) - (v < 0)

            for cols in combinations(range(5), 2):

                def minor(m):
                    return _det([[F(row[j]) for j in cols] for row in m])

                self.assertEqual(sign(minor(a)), sign(minor(b)))
            self.assertEqual(lift_elementary_target(a, b, b[0], 0), tuple(a[0]))
            self.assertIsNone(
                lift_elementary_target(a, b, [v + w for v, w in zip(*b)], 0)
            )
        with self.assertRaises(ValueError):
            small_chord_rows(
                forest_form(PhysicalSource([[1, 0, 1, 1, 1, 1], [0, 1, 1, 1, 1, 1]]))
            )
