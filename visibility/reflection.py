"""Two exact shared visibility certificates from a calibrated rational paired basis."""

from fractions import Fraction

F = Fraction


def _solve(matrix, rhs):
    """Solve a positive definite rational system with multiple right-hand sides."""
    n = len(matrix)
    rows = [list(row) + list(right) for row, right in zip(matrix, rhs)]
    for i in range(n):
        pivot = rows[i][i]
        if pivot <= 0:
            raise ValueError("the supplied rows must be linearly independent")
        for j in range(i + 1, n):
            if rows[j][i]:
                ratio = rows[j][i] / pivot
                for k in range(i, len(rows[j])):
                    rows[j][k] -= ratio * rows[i][k]
    for i in reversed(range(n)):
        for k in range(n, len(rows[i])):
            rows[i][k] = (
                rows[i][k] - sum(rows[i][j] * rows[j][k] for j in range(i + 1, n))
            ) / rows[i][i]
    return [row[n:] for row in rows]


def _power_of_two_at_least(value):
    return 1 << (value - 1).bit_length()


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def verify_rows(basis, error, multipliers):
    """Each returned sign excludes a normalized opposite witness with head >= 1/3."""
    n = len(basis[0]) // 2
    outputs = []
    for side, weights in enumerate(multipliers):
        if len(weights) != len(basis):
            raise ValueError("one multiplier per basis row is required")
        radius = error * sum(map(abs, weights))
        lower = [
            _dot(weights, col) - radius
            for col in zip(*(row[side * n : (side + 1) * n] for row in basis))
        ]
        loss = -3 * sum(min(x, 0) for x in lower)
        outputs.append([x > loss for x in lower])
    return "".join("+" if a else "-" if b else "?" for a, b in zip(*outputs))


def shared_rows(basis, error=0):
    """Return signs, two coefficient vectors, and the number of fixed-point steps.

    basis has 2n columns arranged [positive | polar]. The exact paired product
    lies within error of every entry. Tight calibration as stated in the proof
    ensures a complete answer. With a larger bound a '?' is an explicit abstention.
    """
    basis = [list(map(F, row)) for row in basis]
    if not basis or not basis[0] or len(basis[0]) % 2:
        raise ValueError("expected a nonempty paired basis")
    m, n = len(basis), len(basis[0]) // 2
    if any(len(row) != 2 * n for row in basis) or error < 0:
        raise ValueError("inconsistent basis width or negative error")
    error = F(error)
    gram = [[_dot(a, b) for b in basis] for a in basis]
    weights = _solve(gram, basis)
    grid = _power_of_two_at_least(128 * n)
    scale = _power_of_two_at_least(65536 * n**3)
    denom = scale + n
    # The integer compression is shrunk so the implemented reflection is nonexpansive.
    matrix = [
        [
            int(
                scale
                * (2 * sum(basis[k][i] * weights[k][j] for k in range(m)) - int(i == j))
            )
            for j in range(n)
        ]
        for i in range(n)
    ]
    z = [grid] * n
    total = [0] * n
    total_image = [0] * n
    limit = 64 * n * n
    for step in range(1, limit + 1):
        image = [_dot(row, z) for row in matrix]
        for i in range(n):
            total[i] += z[i]
            total_image[i] += image[i]
            z[i] = max(grid, abs(image[i]) // denom)
        if step != limit and step > 8 and step & (step - 1):
            continue
        both = [
            [denom * a + sign * b for a, b in zip(total, total_image)]
            for sign in (1, -1)
        ]
        losses = [-3 * sum(min(value, 0) for value in row) for row in both]
        if step != limit and any(
            p <= losses[0] and q <= losses[1] for p, q in zip(*both)
        ):
            continue
        certs = [
            [_dot(row[side * n : (side + 1) * n], total) for row in weights]
            for side in range(2)
        ]
        mask = verify_rows(basis, error, certs)
        if "?" not in mask or step == limit:
            return mask, certs, step
    raise AssertionError("unreachable")
