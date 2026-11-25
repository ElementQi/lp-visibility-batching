"""Exact invariants and small separating families for the two query interfaces."""

import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from visibility.axis import AxisTable
from visibility.balanced import balanced_chart, balanced_precision
from visibility.calibration import Minor
from visibility.physical import PhysicalSource, _det, _echelon, _log_proxy


class BalancedTest(unittest.TestCase):
    def test_growth_and_cofactor_count(self):
        r, d = 12, 4
        rows = [
            [F(i == j) if j <= i else F(-1) for j in range(r)] + [F(1)] * d
            for i in range(r)
        ]
        for refresh in (True, False):
            n = r + d
            bits, _ = balanced_precision(n, r, refresh)
            queried = set()

            def minor(cols):
                self.assertNotIn(cols, queried)
                queried.add(cols)
                sub = [[row[j] for j in cols] for row in rows[: len(cols)]]
                v = _det(sub)
                return (
                    Minor(1 if v > 0 else -1, _log_proxy(abs(v), bits + 4))
                    if v
                    else None
                )

            chart = balanced_chart([0] * n, r, minor, refresh=refresh)
            self.assertGreater(chart.exchanges, 0)
            self.assertLessEqual(
                len(queried), r * n * (1 + ((r - 1).bit_length() if refresh else 0))
            )
            self.assertEqual(chart.distinct_requests, len(queried))
            piv, physical = _echelon(rows, chart.pivots + chart.free)
            self.assertEqual(tuple(piv), chart.pivots)
            for a, b in zip(physical, chart.entries):
                for j, value in zip(chart.free, b):
                    self.assertLessEqual(abs(a[j] - value), chart.entry_error)
                    self.assertLessEqual(abs(a[j]), 3)
            self.assertLess(384 * n**3 * chart.entry_error, 1)

    def test_contraction_polar_and_price_lines(self):
        a = [[1, 0, 1, 1], [0, 1, 1, -1]]
        src = PhysicalSource([*a, a[1]], log_precision=100)
        for polar, pinned in (
            (False, ()),
            (True, (0,)),
            (False, (2,)),
            (False, (0, 1)),
            (True, ()),
        ):
            view = src.view(pinned, polar)
            r = len(view.pivots)
            n = len(view.remaining)
            basis = src._dual if polar else src._rows
            piv, mat = _echelon(basis, (*pinned, *view.remaining))
            constrained = sum(j in pinned for j in piv)
            rows = [
                [row[j] for j in view.remaining] for row in mat[constrained : len(piv)]
            ]
            table = AxisTable([0] * n, r, view)
            for i in range(n):
                for shift in (-35, 0, 37):
                    chart = table.line(i).chart(shift)
                    scaled = [
                        [
                            F(x) * F(2) ** (-shift if j == i else 0)
                            for j, x in enumerate(row)
                        ]
                        for row in rows
                    ]
                    _, physical = _echelon(scaled, chart.pivots + chart.free)
                    for row, approx in zip(physical, chart.entries):
                        for col, est in zip(chart.free, approx):
                            self.assertLessEqual(abs(row[col] - est), chart.entry_error)
                far = table.line(i).chart(10**400)
                self.assertTrue(all(abs(v) <= 7 for row in far.entries for v in row))

    def test_sign_identical_triangular_witnesses(self):
        for m in (5, 9, 16):
            a = [F(1, 2 ** min(i, m - 2)) for i in range(m)]
            eye = [[int(i == j) for j in range(m)] for i in range(m)]

            def upper(rows):
                return [
                    [v - sum(s[j] for s in rows[i + 1 :]) for j, v in enumerate(row)]
                    for i, row in enumerate(rows)
                ]

            arrays = []
            for sign in (1, -1):
                block = [
                    [
                        sign * x * y
                        + (
                            2 * x * a[-1]
                            if j == m - 1
                            else -2 * a[-1] * y
                            if i == m - 1
                            else 0
                        )
                        for j, y in enumerate(a)
                    ]
                    for i, x in enumerate(a)
                ]
                plus = [[*e, *v] for e, v in zip(eye, block)]
                minus = [[*(-x for x in v), *e] for v, e in zip(zip(*block), eye)]
                witness = (plus if sign == 1 else minus)[0]
                self.assertEqual(witness[0], 1)
                self.assertEqual(max(witness), 1)
                self.assertGreaterEqual(min(witness), 0)
                arrays.append(upper(plus) + upper(minus))
            diff = []
            for left, right in zip(*arrays):
                for x, y in zip(left, right):
                    self.assertEqual((x > 0) - (x < 0), (y > 0) - (y < 0))
                    diff.append(abs(x - y))
            self.assertEqual(max(diff), 2 * a[-1])

    def test_terminals_and_restart(self):
        def stop(_):
            raise RuntimeError("restart")

        with self.assertRaisesRegex(RuntimeError, "restart"):
            balanced_chart([0, 0, 0], 2, stop)
        self.assertEqual(balanced_chart([10], 1, stop).pivots, (0,))
        self.assertEqual(balanced_chart([10], 0, stop).free, (0,))
        with self.assertRaises(ValueError):
            balanced_chart([0, 0], 1, lambda _: None)
