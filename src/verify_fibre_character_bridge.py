"""Standalone exact checker for the CW(105,36) fibre-character bridge.

This file deliberately imports neither the contraction enumerator nor the
CGW(35,12;6) classification.  Eisenstein integers are represented as pairs
``a + b*w`` with ``w^2 + w + 1 = 0``.
"""

from __future__ import annotations

from itertools import product
from typing import TypeAlias

Eisenstein: TypeAlias = tuple[int, int]

ZERO: Eisenstein = (0, 0)
ONE: Eisenstein = (1, 0)
OMEGA: Eisenstein = (0, 1)


def add(left: Eisenstein, right: Eisenstein) -> Eisenstein:
    return (left[0] + right[0], left[1] + right[1])


def neg(value: Eisenstein) -> Eisenstein:
    return (-value[0], -value[1])


def mul(left: Eisenstein, right: Eisenstein) -> Eisenstein:
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c - b * d)


def conjugate(value: Eisenstein) -> Eisenstein:
    a, b = value
    return (a - b, -b)


def power(value: Eisenstein, exponent: int) -> Eisenstein:
    result = ONE
    for _ in range(exponent):
        result = mul(result, value)
    return result


def scalar(value: int) -> Eisenstein:
    return (value, 0)


def character(column: tuple[int, int, int]) -> Eisenstein:
    omega2 = mul(OMEGA, OMEGA)
    return add(
        add(scalar(column[0]), mul(scalar(column[1]), OMEGA)),
        mul(scalar(column[2]), omega2),
    )


def contraction_row(positive_positions: tuple[int, int, int]) -> tuple[int, ...]:
    row = [0] * 35
    row[0] = -3
    for position in positive_positions:
        row[position] = 3
    return tuple(row)


EXPECTED_CONTRACTION_ROWS = (
    contraction_row((5, 10, 20)),
    contraction_row((15, 25, 30)),
)


def periodic_correlation(row: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(row[q] * row[(q + shift) % len(row)] for q in range(len(row)))
        for shift in range(len(row))
    )


def run_exact_checks() -> dict[str, object]:
    # The CRT map (q,r) -> 3q+35r is checked as a bijective group map.
    crt_images: set[int] = set()
    for q in range(35):
        for r in range(3):
            image = (3 * q + 35 * r) % 105
            assert (12 * image) % 35 == q
            assert (2 * image) % 3 == r
            crt_images.add(image)
    assert crt_images == set(range(105))
    for q, r, t, s in product(range(35), range(3), range(35), range(3)):
        left = (3 * ((q + t) % 35) + 35 * ((r + s) % 3)) % 105
        right = (3 * q + 35 * r + 3 * t + 35 * s) % 105
        assert left == right

    # Coefficient-by-coefficient proof of
    # sum_q D_q conjugate(D_{q+t}) = sum_s w^{-s} R(t,s).
    for r, u in product(range(3), repeat=2):
        left = mul(power(OMEGA, r), conjugate(power(OMEGA, u)))
        fibre_shift = (u - r) % 3
        right = power(OMEGA, (-fibre_shift) % 3)
        assert left == right

    legal_columns = tuple(product((-1, 0, 1), repeat=3))
    zero_sum_columns = tuple(column for column in legal_columns if sum(column) == 0)
    active_columns = tuple(
        column for column in zero_sum_columns if sum(value * value for value in column) == 2
    )
    assert len(legal_columns) == 27
    assert set(zero_sum_columns) == {(-1, 0, 1), (-1, 1, 0), (0, -1, 1),
                                     (0, 0, 0), (0, 1, -1), (1, -1, 0),
                                     (1, 0, -1)}
    assert len(active_columns) == 6
    assert tuple(column for column in legal_columns if sum(column) == 3) == ((1, 1, 1),)
    assert tuple(column for column in legal_columns if sum(column) == -3) == ((-1, -1, -1),)

    one_minus_omega = add(ONE, neg(OMEGA))
    zeta6 = add(ONE, OMEGA)
    states = (
        (1, -1, 0),
        (1, 0, -1),
        (0, 1, -1),
        (-1, 1, 0),
        (-1, 0, 1),
        (0, -1, 1),
    )
    assert set(states) == set(active_columns)
    for exponent, column in enumerate(states):
        assert character(column) == mul(one_minus_omega, power(zeta6, exponent))
    assert character((0, 0, 0)) == ZERO
    assert character((1, 1, 1)) == ZERO
    assert character((-1, -1, -1)) == ZERO

    norm = mul(one_minus_omega, conjugate(one_minus_omega))
    assert norm == scalar(3)

    active_counts: list[int] = []
    forced_counts: list[int] = []
    for row in EXPECTED_CONTRACTION_ROWS:
        assert sum(row) == 6
        assert sum(value * value for value in row) == 36
        assert periodic_correlation(row) == (36,) + (0,) * 34
        assert all(row[(4 * q) % 35] == row[q] for q in range(35))
        assert set(row) <= {-3, 0, 3}

        forced_count = sum(value != 0 for value in row)
        forced_weight = 3 * forced_count
        remaining_weight = 36 - forced_weight
        assert remaining_weight >= 0 and remaining_weight % 2 == 0
        active_count = remaining_weight // 2
        assert forced_count == 4
        assert active_count == 12
        forced_counts.append(forced_count)
        active_counts.append(active_count)

    # If AA#=36 on Z_35 x Z_3, then R(0,0)=36 and every other
    # R(t,s)=0.  The coefficient identity above therefore gives DD#=36.
    projected = []
    for t in range(35):
        value = ZERO
        for fibre_shift in range(3):
            correlation = 36 if (t, fibre_shift) == (0, 0) else 0
            value = add(
                value,
                mul(
                    scalar(correlation),
                    power(OMEGA, (-fibre_shift) % 3),
                ),
            )
        projected.append(value)
    assert tuple(projected) == (scalar(36),) + (ZERO,) * 34

    # Since D=(1-w)U and N(1-w)=3, DD#=36 is exactly UU#=12.
    divided = tuple((value[0] // 3, value[1] // 3) for value in projected)
    assert all(value[0] % 3 == value[1] % 3 == 0 for value in projected)
    assert divided == (scalar(12),) + (ZERO,) * 34

    return {
        "crt_points_checked": len(crt_images),
        "character_coefficients_checked": 9,
        "legal_columns_checked": len(legal_columns),
        "active_character_states_checked": len(states),
        "norm_one_minus_omega": norm[0],
        "contraction_rows_checked": len(EXPECTED_CONTRACTION_ROWS),
        "active_columns_by_row": active_counts,
        "forced_constant_columns_by_row": forced_counts,
        "contraction_correlations_checked": 35 * len(EXPECTED_CONTRACTION_ROWS),
        "projected_zero_shift": projected[0][0],
        "quotient_zero_shift": divided[0][0],
    }


def main() -> None:
    report = run_exact_checks()
    for key, value in report.items():
        if isinstance(value, list):
            rendered = ",".join(str(item) for item in value)
        else:
            rendered = str(value)
        print(f"{key}={rendered}")


if __name__ == "__main__":
    main()
