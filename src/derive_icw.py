"""Independent reproduction of the multiplier-4 ICW_3(35,36) reduction.

This module uses two separate finite enumerators:

1. recursive depth-first search with interval/gcd pruning;
2. meet-in-the-middle enumeration joined on the two scalar invariants.

They enumerate multiplier-orbit coefficient vectors in [-3,3] satisfying
sum b_i = 6 and sum b_i^2 = 36. A separate exact autocorrelation filter then
finds the invariant ICWs and classifies them under affine/sign equivalence.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from functools import reduce
from pathlib import Path
from typing import Iterable, Sequence

N = 35
MULTIPLIER = 4
COEFFICIENT_BOUND = 3
TARGET_SUM = 6
TARGET_WEIGHT = 36


@dataclass(frozen=True)
class ICWEnumeration:
    orbits: tuple[tuple[int, ...], ...]
    scalar_candidates_dfs: tuple[tuple[int, ...], ...]
    scalar_candidates_mitm: tuple[tuple[int, ...], ...]
    invariant_solutions: tuple[tuple[int, ...], ...]
    equivalence_classes: tuple[tuple[int, ...], ...]


def multiplier_orbits(n: int = N, multiplier: int = MULTIPLIER) -> tuple[tuple[int, ...], ...]:
    """Return sorted orbits of index multiplication by ``multiplier`` mod n."""
    if math.gcd(n, multiplier) != 1:
        raise ValueError("multiplier must be a unit modulo n")
    unseen = set(range(n))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit: list[int] = []
        current = start
        while current not in orbit:
            orbit.append(current)
            unseen.discard(current)
            current = (multiplier * current) % n
        orbits.append(tuple(orbit))
    return tuple(orbits)


def multiplier_order(n: int = N, multiplier: int = MULTIPLIER) -> int:
    if math.gcd(n, multiplier) != 1:
        raise ValueError("multiplier must be a unit modulo n")
    value = 1
    for order in range(1, n + 1):
        value = value * multiplier % n
        if value == 1:
            return order
    raise RuntimeError("failed to find multiplicative order")


def _reachable_linear_bounds(sizes: Sequence[int]) -> tuple[int, int]:
    total = COEFFICIENT_BOUND * sum(sizes)
    return -total, total


def _suffix_gcds(sizes: Sequence[int]) -> tuple[int, ...]:
    result = [0] * (len(sizes) + 1)
    for i in range(len(sizes) - 1, -1, -1):
        result[i] = math.gcd(sizes[i], result[i + 1])
    return tuple(result)


def enumerate_scalar_candidates_dfs(
    orbit_sizes: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    """DFS enumeration of orbit coefficients satisfying sum and norm only."""
    sizes = tuple(int(s) for s in orbit_sizes)
    suffix_sum = [0] * (len(sizes) + 1)
    for i in range(len(sizes) - 1, -1, -1):
        suffix_sum[i] = suffix_sum[i + 1] + sizes[i]
    suffix_gcd = _suffix_gcds(sizes)
    values = tuple(range(-COEFFICIENT_BOUND, COEFFICIENT_BOUND + 1))
    found: list[tuple[int, ...]] = []
    prefix: list[int] = []

    def recurse(index: int, linear: int, square: int) -> None:
        remaining_size = suffix_sum[index]
        max_linear = COEFFICIENT_BOUND * remaining_size
        if not (linear - max_linear <= TARGET_SUM <= linear + max_linear):
            return
        if square > TARGET_WEIGHT:
            return
        if square + COEFFICIENT_BOUND**2 * remaining_size < TARGET_WEIGHT:
            return
        gcd_remaining = suffix_gcd[index]
        if gcd_remaining and (TARGET_SUM - linear) % gcd_remaining:
            return

        if index == len(sizes):
            if linear == TARGET_SUM and square == TARGET_WEIGHT:
                found.append(tuple(prefix))
            return

        size = sizes[index]
        for coefficient in values:
            prefix.append(coefficient)
            recurse(
                index + 1,
                linear + size * coefficient,
                square + size * coefficient * coefficient,
            )
            prefix.pop()

    recurse(0, 0, 0)
    return tuple(sorted(found))


def _enumerate_half(sizes: Sequence[int]) -> list[tuple[tuple[int, ...], int, int]]:
    values = range(-COEFFICIENT_BOUND, COEFFICIENT_BOUND + 1)
    records: list[tuple[tuple[int, ...], int, int]] = []

    def recurse(index: int, coefficients: list[int], linear: int, square: int) -> None:
        if square > TARGET_WEIGHT:
            return
        if index == len(sizes):
            records.append((tuple(coefficients), linear, square))
            return
        size = sizes[index]
        for coefficient in values:
            coefficients.append(coefficient)
            recurse(
                index + 1,
                coefficients,
                linear + size * coefficient,
                square + size * coefficient * coefficient,
            )
            coefficients.pop()

    recurse(0, [], 0, 0)
    return records


def enumerate_scalar_candidates_mitm(
    orbit_sizes: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    """Meet-in-the-middle enumeration using independent half tables."""
    sizes = tuple(int(s) for s in orbit_sizes)
    split = len(sizes) // 2
    left_records = _enumerate_half(sizes[:split])
    right_records = _enumerate_half(sizes[split:])

    right_index: dict[tuple[int, int], list[tuple[int, ...]]] = defaultdict(list)
    for coefficients, linear, square in right_records:
        right_index[(linear, square)].append(coefficients)

    found: list[tuple[int, ...]] = []
    for left, linear, square in left_records:
        key = (TARGET_SUM - linear, TARGET_WEIGHT - square)
        for right in right_index.get(key, ()):
            found.append(left + right)
    return tuple(sorted(found))


def sequence_from_orbit_coefficients(
    coefficients: Sequence[int], orbits: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    if len(coefficients) != len(orbits):
        raise ValueError("one coefficient is required for every orbit")
    row = [0] * sum(len(orbit) for orbit in orbits)
    for coefficient, orbit in zip(coefficients, orbits, strict=True):
        for index in orbit:
            row[index] = int(coefficient)
    return tuple(row)


def integer_periodic_autocorrelation(values: Sequence[int]) -> tuple[int, ...]:
    """Local ICW correlation implementation, independent of verify_direct.py."""
    row = tuple(int(value) for value in values)
    n = len(row)
    return tuple(
        sum(row[index] * row[(index + shift) % n] for index in range(n))
        for shift in range(n)
    )


def is_icw(values: Sequence[int], weight: int = TARGET_WEIGHT) -> bool:
    correlations = integer_periodic_autocorrelation(values)
    return correlations[0] == weight and all(value == 0 for value in correlations[1:])


def units_modulo(n: int = N) -> tuple[int, ...]:
    return tuple(value for value in range(n) if math.gcd(value, n) == 1)


def affine_transform(
    values: Sequence[int], multiplier: int, shift: int, sign: int = 1
) -> tuple[int, ...]:
    """Map old index i to multiplier*i+shift and optionally negate values."""
    row = tuple(int(value) for value in values)
    n = len(row)
    if math.gcd(multiplier, n) != 1:
        raise ValueError("multiplier must be a unit modulo the sequence length")
    if sign not in (-1, 1):
        raise ValueError("sign must be -1 or 1")
    transformed = [0] * n
    for index, value in enumerate(row):
        transformed[(multiplier * index + shift) % n] = sign * value
    return tuple(transformed)


def canonical_affine_sign(values: Sequence[int]) -> tuple[int, ...]:
    row = tuple(int(value) for value in values)
    n = len(row)
    return min(
        affine_transform(row, multiplier, shift, sign)
        for multiplier in units_modulo(n)
        for shift in range(n)
        for sign in (-1, 1)
    )


def classify_affine_sign(
    rows: Iterable[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted({canonical_affine_sign(row) for row in rows}))


def normalized_b() -> tuple[int, ...]:
    row = [0] * N
    row[0] = -3
    for index in (5, 10, 20):
        row[index] = 3
    return tuple(row)


def enumerate_icw() -> ICWEnumeration:
    orbits = multiplier_orbits()
    sizes = tuple(len(orbit) for orbit in orbits)
    dfs = enumerate_scalar_candidates_dfs(sizes)
    mitm = enumerate_scalar_candidates_mitm(sizes)
    if dfs != mitm:
        raise AssertionError("DFS and meet-in-the-middle scalar enumerations disagree")

    invariant_solutions = tuple(
        sequence_from_orbit_coefficients(coefficients, orbits)
        for coefficients in dfs
        if is_icw(sequence_from_orbit_coefficients(coefficients, orbits))
    )
    classes = classify_affine_sign(invariant_solutions)
    return ICWEnumeration(
        orbits=orbits,
        scalar_candidates_dfs=dfs,
        scalar_candidates_mitm=mitm,
        invariant_solutions=tuple(sorted(invariant_solutions)),
        equivalence_classes=classes,
    )


def _nonzero_map(row: Sequence[int]) -> dict[str, int]:
    return {str(i): int(value) for i, value in enumerate(row) if value}


def enumeration_report() -> dict[str, object]:
    result = enumerate_icw()
    normalized = normalized_b()
    normalized_present = normalized in result.invariant_solutions
    alternate = affine_transform(normalized, 3, 0, 1)
    return {
        "parameters": {
            "n": N,
            "weight": TARGET_WEIGHT,
            "coefficient_bound": COEFFICIENT_BOUND,
            "target_sum": TARGET_SUM,
            "multiplier": MULTIPLIER,
            "multiplier_order": multiplier_order(),
        },
        "orbits": [list(orbit) for orbit in result.orbits],
        "orbit_sizes": [len(orbit) for orbit in result.orbits],
        "scalar_candidate_count_dfs": len(result.scalar_candidates_dfs),
        "scalar_candidate_count_mitm": len(result.scalar_candidates_mitm),
        "enumerators_agree_exactly": result.scalar_candidates_dfs
        == result.scalar_candidates_mitm,
        "invariant_solution_count": len(result.invariant_solutions),
        "invariant_solutions": [
            {
                "nonzero_entries": _nonzero_map(row),
                "autocorrelations": list(integer_periodic_autocorrelation(row)),
            }
            for row in result.invariant_solutions
        ],
        "affine_sign_equivalence_class_count": len(result.equivalence_classes),
        "normalized_candidate": {
            "polynomial": "-3 + 3 X^5 + 3 X^10 + 3 X^20 (mod X^35-1)",
            "nonzero_entries": _nonzero_map(normalized),
            "present_in_invariant_solutions": normalized_present,
            "autocorrelations": list(integer_periodic_autocorrelation(normalized)),
        },
        "explicit_equivalence": {
            "unit_multiplier": 3,
            "image_of_normalized_nonzero_entries": _nonzero_map(alternate),
            "explanation": "multiplication by 3 maps {5,10,20} to {15,30,25}",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="optional JSON report path")
    args = parser.parse_args()
    report = enumeration_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
