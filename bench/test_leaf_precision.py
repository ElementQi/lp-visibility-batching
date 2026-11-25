"""A balanced leaf transcript at reduced precision with opposite mandatory signs."""

from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from visibility.balanced import balanced_chart
from visibility.calibration import Minor
from visibility.physical import PhysicalSource, _echelon


def staircase(height):
    """Return sparse T and a with T a = a[-1] e_last."""
    if not height:
        return [{0: 1}], [F(1)]
    tail, v = staircase(height - 1)
    half = len(tail)
    width = max(1, half // 2)
    left = [{i: 1, **{half + j: -1 for j in range(width)}} for i in range(half)]
    right = [{half + j: x for j, x in row.items()} for row in tail]
    return left + right, [F(1)] * half + [x / width for x in v]


def sign(x):
    return (x > 0) - (x < 0)


class LeafPrecisionTest(unittest.TestCase):
    def test_exponent_and_signs(self):
        for height in range(4, 10):
            t, a = staircase(height)
            r = len(t)
            loss = (height - 2) * (height - 1) // 2
            self.assertEqual(a[-1], F(1, 2**loss))
            center = [F(0)] * (r - 1) + [-2 * a[-1]]
            images = []
            for choice in (1, -1):
                v = [x + choice * y for x, y in zip(center, a)]
                images.append([sum(x * v[j] for j, x in row.items()) for row in t])
                if choice > 0:
                    self.assertEqual(v[0], 1)
                else:
                    self.assertTrue(all(-1 <= y <= 0 for y in v))
            for i, (x, y) in enumerate(zip(*images)):
                self.assertEqual(x - y, 2 * a[-1] if i == r - 1 else 0)
                self.assertEqual(sign(x), sign(y))
                self.assertLessEqual(max(abs(x), abs(y)), 3 * a[-1])

    def test_actual_source_and_adaptive_merges(self):
        # Height four already distinguishes the one-bit leaf data. Higher
        # heights above test the exact quasipolynomial growth law.
        t, a = staircase(4)
        r = len(a)
        price = 2 * r
        charts = []
        for choice in (1, -1):
            v = [choice * y - (2 * a[-1] if j == r - 1 else 0) for j, y in enumerate(a)]
            rows = [
                [int(j == i) for j in range(r)]
                + [row.get(j, 0) for j in range(r)]
                + [sum(x * v[j] for j, x in row.items())]
                for i, row in enumerate(t)
            ]
            source = PhysicalSource(rows)
            # Keep the upstream source precise. Only lower the precision used
            # to clip normalized leaf data in this deliberately weakened run.
            with patch("visibility.balanced.balanced_precision", return_value=(1, 20)):
                chart = balanced_chart(
                    [price] * r + [0] * (r + 1), r, source.view(), refresh=False
                )
            charts.append(chart)
            self.assertEqual(chart.exchanges, 0)
            self.assertEqual(chart.pivots, tuple(range(r, 2 * r)))
            effective = [
                [F(x, 2**price) if j < r else x for j, x in enumerate(row)]
                for row in rows
            ]
            _, exact = _echelon(effective, (*chart.pivots, *chart.free))
            if choice > 0:
                positive = exact[0]
                self.assertEqual(positive[r], 1)
                self.assertEqual(max(positive), 1)
                self.assertGreaterEqual(min(positive), 0)
            else:
                polar = [0] * r + [-x for x in v] + [1]
                self.assertEqual(polar[r], 1)
                self.assertEqual(max(polar), 1)
                self.assertGreaterEqual(min(polar), 0)
                self.assertTrue(
                    all(
                        sum(x * y for x, y in zip(row, polar)) == 0 for row in effective
                    )
                )
        self.assertEqual(charts[0], charts[1])
        self.assertTrue(all(v == 0 for row in charts[0].entries for v in row))

    def test_common_log_transcript_away_from_clipping(self):
        # Both true terminal coefficients are close to -1/2. A common signed
        # log reply on the leaf neighborhood is admissible for K=2.
        t, a = staircase(5)
        r = len(a)
        price = 2 * r
        traces = []
        for choice in (1, -1):
            v = [choice * x - (F(1, 2) if j == r - 1 else 0) for j, x in enumerate(a)]
            rows = [
                [int(i == j) for j in range(r)]
                + [row.get(j, 0) for j in range(r)]
                + [sum(x * v[j] for j, x in row.items())]
                for i, row in enumerate(t)
            ]
            source = PhysicalSource(rows)
            transcript = {}

            def common_minor(cols):
                physical = source.minor(cols)
                reply = physical
                if physical and abs(physical.proxy + 1) < F(1, 8):
                    reply = Minor(physical.sign, F(-1))
                transcript[cols] = reply
                return reply

            from visibility.minor_view import MinorView

            view = MinorView(
                2 * r + 1,
                r,
                False,
                (),
                tuple(range(2 * r + 1)),
                tuple(range(r)),
                common_minor,
            )
            with patch("visibility.balanced.balanced_precision", return_value=(4, 10)):
                chart = balanced_chart(
                    [price] * r + [0] * (r + 1), r, view, refresh=False
                )
            self.assertEqual(chart.pivots, tuple(range(r, 2 * r)))
            self.assertEqual(chart.exchanges, 0)
            self.assertEqual(
                [row[-1] for row in chart.entries], [0] * (r - 1) + [F(-1, 2)]
            )
            self.assertEqual(v[0], choice)
            traces.append(transcript)
        self.assertEqual(traces[0], traces[1])
