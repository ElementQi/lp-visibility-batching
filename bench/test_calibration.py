"""Exact integration checks for graph selection, chart rounding, and visibility."""

import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from visibility.calibration import (
    Minor,
    _normalized_power,
    calibrated_chart,
    price_precision,
)
from visibility.reflection import shared_rows, verify_rows


def determinant(rows):
    if not rows:
        return 1
    return sum(
        (-1) ** k * x * determinant([row[:k] + row[k + 1 :] for row in rows[1:]])
        for k, x in enumerate(rows[0])
        if x
    )


class CalibrationTest(unittest.TestCase):
    def test_dyadic_powers(self):
        for e in (4, 8, 16, 32):
            for denominator in (1, 3, 8, 31):
                for numerator in range(-e * denominator - 1, 1, 3):
                    q = F(numerator, denominator)
                    x = _normalized_power(q, e)
                    error = F(1, 2**e)
                    # Exact endpoint tests after raising to a positive integer.
                    actual = F(2) ** q.numerator
                    self.assertLessEqual(max(F(0), x - error) ** q.denominator, actual)
                    self.assertGreaterEqual((x + error) ** q.denominator, actual)

    def test_hidden_gauge_and_exchange(self):
        # Every nonzero minor of this matrix has power-of-two magnitude except
        # some unused prefix minors. The callback computes a dyadic log proxy
        # exactly when requested at power-of-two minors.
        a = [[1, -1, -1, 1], [0, 1, -1, 1], [0, 0, 1, 1]]
        observed = []

        def minor(cols):
            observed.append(cols)
            x = determinant([[row[j] for j in cols] for row in a[: len(cols)]])
            if not x:
                return None
            m = abs(x)
            if m & (m - 1):
                raise AssertionError((cols, x))
            # Multiplying column j by 2^(3*j) is cancelled by price 3*j.
            return Minor(1 if x > 0 else -1, F(m.bit_length() - 1 + 3 * sum(cols)))

        chart = calibrated_chart([0, 3, 6, 9], 3, minor)
        self.assertEqual(chart.exchanges, 1)
        self.assertEqual(chart.distinct_requests, len(set(observed)))
        self.assertEqual(len(observed), len(set(observed)))
        for k, row in enumerate(a):
            for col, head in enumerate(chart.free):
                error = abs(
                    sum(
                        row[j] * coeff[col]
                        for j, coeff in zip(chart.pivots, chart.entries)
                    )
                    - row[head]
                )
                self.assertLessEqual(
                    error, chart.entry_error * sum(abs(row[j]) for j in chart.pivots)
                )

    def test_restart_and_terminal_ranks(self):
        class Restart(Exception):
            pass

        def restart(_):
            raise Restart()

        with self.assertRaises(Restart):
            calibrated_chart([0, 1], 1, restart)
        self.assertEqual(calibrated_chart([0, 1], 0, restart).free, (0, 1))
        self.assertEqual(calibrated_chart([0, 1], 2, restart).pivots, (0, 1))
        for n in (1, 2, 128, 4096):
            e, k = price_precision(n)
            self.assertGreaterEqual(2**k, 2 * n * (e + 5) * 2**e)
            self.assertLess(6 * n * n * F(96, 2**e), F(1, 64 * n))

    def test_shared_rows(self):
        for t in (F(-3), F(-1, 3), F(0), F(1, 3), F(3)):
            basis = [[1, t, 0, 0], [0, 0, -t, 1]]
            answer, rows, _ = shared_rows(basis)
            self.assertEqual(verify_rows(basis, 0, rows), answer)
            self.assertNotIn("?", answer)
            if t > 0:
                self.assertEqual(answer, "++")
            if t < 0:
                self.assertEqual(answer, "--")
            if t == 0:
                self.assertEqual(answer, "+-")


if __name__ == "__main__":
    unittest.main()
