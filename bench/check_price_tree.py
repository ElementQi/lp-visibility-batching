"""Rational certificates for reduced fixed schedules in the extremal-price learner."""

from fractions import Fraction as F
from math import ceil


def check(k, delta, gamma, beta, eta, reserve, decay, c, eps, sigma, r, exponent):
    # Exact rational log bounds from integrating the geometric series, 0<eta<1.
    hit = 1 - (1 - gamma) ** k
    low_plus = eta - eta**2 / 2 + eta**3 / 3 - eta**4 / 4
    low_minus = -eta - eta**2 / (2 * (1 - eta))
    assert hit * low_plus + (1 - hit) * low_minus >= eps
    ratio = beta / (1 + (beta - 1) * delta)
    C = ratio**k * max(F(1), 1 / (reserve * (delta - gamma)), ratio / decay)
    mass = 1 + reserve / (1 - decay)
    assert low_minus + sigma * r / C >= eps
    assert c.numerator == 1 and c.denominator & (c.denominator - 1) == 0
    log_fee = c.denominator.bit_length() * F(7, 10)  # ln(2/c) <= log2(2/c)*.7
    assert (eps / 2 - sigma - mass * c / (1 - c)) * r >= mass * log_fee
    assert beta**exponent > r * r * (1 + (beta - 1) * delta) ** (2 * exponent)
    root_slope = 2 * ceil(mass / (2 * sigma))
    assert root_slope * sigma >= mass and r % 2
    return C, mass, root_slope


if __name__ == "__main__":
    options = [
        (
            4,
            F(2, 11),
            F(3, 17),
            F(4),
            F(1, 16),
            F(4),
            F(1, 4),
            F(1, 2**14),
            F(3, 1000),
            F(7, 10000),
            205001,
            48,
        ),
        (
            5,
            F(5, 33),
            F(7, 48),
            F(5),
            F(1, 16),
            F(12),
            F(1, 3),
            F(1, 2**18),
            F(36, 10000),
            F(88, 100000),
            350001,
            39,
        ),
    ]
    for args in options:
        C, mass, slope = check(*args)
        print(
            "arity",
            args[0],
            "branch",
            args[-2],
            "exponent",
            args[-1],
            "transfer",
            C,
            "mass",
            mass,
            "even root slope",
            slope,
        )
