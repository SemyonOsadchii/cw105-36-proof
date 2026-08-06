"""Finite-field congruence obstruction for the CW(105,36) phase lift.

The exact phase equation is

    U(X) conjugate(U(X^-1)) = 12

in Z[omega][X]/(X^35-1).  Reducing modulo the Eisenstein primes
``1-omega`` and ``2`` gives zero-product equations over F_3 and F_4.
Because 3 and 2 are coprime to 35, the group algebras are semisimple and the
factor pairs under reciprocal conjugation give a small union of cyclic codes.

This module implements the complete proof computation independently of the
standalone C++ checker.  All field arithmetic is exact.
"""

from __future__ import annotations

import itertools
import math
from functools import lru_cache
from collections import Counter
from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np

from .phase_reduction import Q_MODULUS, SPECIAL_POSITIONS
from .symmetry import CONTRACTION_STABILIZER, canonical_support


# ---------------------------------------------------------------------------
# Generic low-degree polynomial helpers, coefficients low-to-high.


def _trim(poly: Sequence[int]) -> tuple[int, ...]:
    values = list(int(value) for value in poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values or [0])


def _poly_add(
    left: Sequence[int],
    right: Sequence[int],
    add: Callable[[int, int], int],
) -> tuple[int, ...]:
    length = max(len(left), len(right))
    return _trim(
        tuple(
            add(
                int(left[index]) if index < len(left) else 0,
                int(right[index]) if index < len(right) else 0,
            )
            for index in range(length)
        )
    )


def _poly_mul(
    left: Sequence[int],
    right: Sequence[int],
    add: Callable[[int, int], int],
    multiply: Callable[[int, int], int],
) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = add(result[i + j], multiply(int(a), int(b)))
    return _trim(result)


def _poly_divmod(
    dividend: Sequence[int],
    divisor: Sequence[int],
    add: Callable[[int, int], int],
    multiply: Callable[[int, int], int],
    inverse: Callable[[int], int],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    remainder = list(_trim(dividend))
    divisor_t = _trim(divisor)
    if divisor_t == (0,):
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [0] * max(1, len(remainder) - len(divisor_t) + 1)
    while len(remainder) >= len(divisor_t) and any(remainder):
        degree = len(remainder) - len(divisor_t)
        coefficient = multiply(remainder[-1], inverse(divisor_t[-1]))
        quotient[degree] = coefficient
        # In both F3 and characteristic two, subtraction is implemented by
        # adding the additive inverse.  The caller's add/multiply operations
        # below are paired with explicit negation where needed.
        for index, value in enumerate(divisor_t):
            product = multiply(coefficient, value)
            if add is _gf3_add:
                remainder[degree + index] = (remainder[degree + index] - product) % 3
            else:
                remainder[degree + index] = add(remainder[degree + index], product)
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return _trim(quotient), _trim(remainder)


def _poly_mod(
    dividend: Sequence[int],
    modulus: Sequence[int],
    add: Callable[[int, int], int],
    multiply: Callable[[int, int], int],
    inverse: Callable[[int], int],
) -> tuple[int, ...]:
    return _poly_divmod(dividend, modulus, add, multiply, inverse)[1]


def _poly_gcd(
    left: Sequence[int],
    right: Sequence[int],
    add: Callable[[int, int], int],
    multiply: Callable[[int, int], int],
    inverse: Callable[[int], int],
) -> tuple[int, ...]:
    a, b = _trim(left), _trim(right)
    while b != (0,):
        a, b = b, _poly_mod(a, b, add, multiply, inverse)
    if a == (0,):
        return a
    scale = inverse(a[-1])
    return _trim(tuple(multiply(scale, value) for value in a))


def _poly_pow_mod(
    base: Sequence[int],
    exponent: int,
    modulus: Sequence[int],
    add: Callable[[int, int], int],
    multiply: Callable[[int, int], int],
    inverse: Callable[[int], int],
) -> tuple[int, ...]:
    result: tuple[int, ...] = (1,)
    current = _poly_mod(base, modulus, add, multiply, inverse)
    power = int(exponent)
    while power:
        if power & 1:
            result = _poly_mod(
                _poly_mul(result, current, add, multiply),
                modulus,
                add,
                multiply,
                inverse,
            )
        current = _poly_mod(
            _poly_mul(current, current, add, multiply),
            modulus,
            add,
            multiply,
            inverse,
        )
        power >>= 1
    return result


def _prime_divisors(value: int) -> tuple[int, ...]:
    n = int(value)
    result: list[int] = []
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            result.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        result.append(n)
    return tuple(result)


def _is_irreducible(
    polynomial: Sequence[int],
    field_order: int,
    add: Callable[[int, int], int],
    multiply: Callable[[int, int], int],
    inverse: Callable[[int], int],
) -> bool:
    f = _trim(polynomial)
    degree = len(f) - 1
    if degree <= 0 or f[-1] != 1:
        return False
    x = (0, 1)
    x_mod = _poly_mod(x, f, add, multiply, inverse)
    if _poly_pow_mod(x, field_order**degree, f, add, multiply, inverse) != x_mod:
        return False
    for prime in _prime_divisors(degree):
        power = _poly_pow_mod(
            x, field_order ** (degree // prime), f, add, multiply, inverse
        )
        difference = _poly_add(power, x_mod, add)
        if add is _gf3_add:
            difference = _poly_add(
                power, tuple((-v) % 3 for v in x_mod), add
            )
        if len(_poly_gcd(f, difference, add, multiply, inverse)) > 1:
            return False
    return True


# ---------------------------------------------------------------------------
# F3 reduction modulo (1-omega).


def _gf3_add(left: int, right: int) -> int:
    return (left + right) % 3


def _gf3_mul(left: int, right: int) -> int:
    return (left * right) % 3


def _gf3_inv(value: int) -> int:
    reduced = value % 3
    if reduced == 1:
        return 1
    if reduced == 2:
        return 2
    raise ZeroDivisionError("zero has no inverse in F3")


F3_SELF_FACTORS: tuple[tuple[int, ...], ...] = (
    (2, 1),  # X-1
    (1, 1, 1, 1, 1),  # Phi_5
    (1, 1, 1, 1, 1, 1, 1),  # Phi_7
)
F3_RECIPROCAL_PAIR: tuple[tuple[int, ...], tuple[int, ...]] = (
    (1, 2, 2, 1, 2, 1, 0, 1, 2, 0, 1, 0, 1),
    (1, 0, 1, 0, 2, 1, 0, 1, 2, 1, 2, 2, 1),
)


def gf3_factorization_is_valid() -> bool:
    factors = (*F3_SELF_FACTORS, *F3_RECIPROCAL_PAIR)
    product: tuple[int, ...] = (1,)
    for factor in factors:
        product = _poly_mul(product, factor, _gf3_add, _gf3_mul)
    target = (2,) + (0,) * 34 + (1,)
    return (
        product == target
        and F3_RECIPROCAL_PAIR[0][::-1] == F3_RECIPROCAL_PAIR[1]
        and all(
            _is_irreducible(factor, 3, _gf3_add, _gf3_mul, _gf3_inv)
            for factor in factors
        )
    )


def gf3_code_generators() -> tuple[tuple[int, ...], tuple[int, ...]]:
    mandatory: tuple[int, ...] = (1,)
    for factor in F3_SELF_FACTORS:
        mandatory = _poly_mul(mandatory, factor, _gf3_add, _gf3_mul)
    generators = tuple(
        _poly_mul(mandatory, factor, _gf3_add, _gf3_mul)
        for factor in F3_RECIPROCAL_PAIR
    )
    if any(len(generator) != 24 for generator in generators):
        raise AssertionError("F3 generators must have degree 23")
    return generators  # type: ignore[return-value]


@dataclass(frozen=True)
class Mod3Candidate:
    code_index: int
    message: tuple[int, ...]
    sequence: tuple[int, ...]

    @property
    def support(self) -> tuple[int, ...]:
        return tuple(index for index, value in enumerate(self.sequence) if value)

    @property
    def support_mask(self) -> int:
        return sum(1 << index for index in self.support)


def _all_ternary_messages(length: int) -> np.ndarray:
    count = 3**length
    integers = np.arange(count, dtype=np.int64)
    messages = np.empty((count, length), dtype=np.int16)
    for column in range(length):
        messages[:, column] = integers % 3
        integers //= 3
    return messages


@lru_cache(maxsize=1)
def enumerate_mod3_candidates_numpy() -> tuple[tuple[Mod3Candidate, ...], ...]:
    """Enumerate both 3^12 cyclic codes with exact integer NumPy arithmetic."""
    messages = _all_ternary_messages(12)
    result: list[tuple[Mod3Candidate, ...]] = []
    for code_index, generator in enumerate(gf3_code_generators()):
        matrix = np.zeros((12, Q_MODULUS), dtype=np.int16)
        generator_array = np.asarray(generator, dtype=np.int16)
        for shift in range(12):
            matrix[shift, shift : shift + len(generator)] = generator_array
        codewords = (messages @ matrix) % 3
        weight_mask = np.count_nonzero(codewords, axis=1) == 12
        special_mask = np.all(codewords[:, sorted(SPECIAL_POSITIONS)] == 0, axis=1)
        selected = np.flatnonzero(weight_mask & special_mask)
        candidates = tuple(
            Mod3Candidate(
                code_index=code_index,
                message=tuple(int(value) for value in messages[index]),
                sequence=tuple(int(value) for value in codewords[index]),
            )
            for index in selected
        )
        result.append(candidates)
    return tuple(result)


def mod3_periodic_autocorrelation(sequence: Sequence[int]) -> tuple[int, ...]:
    row = tuple(int(value) % 3 for value in sequence)
    return tuple(
        sum(row[q] * row[(q + shift) % len(row)] for q in range(len(row))) % 3
        for shift in range(len(row))
    )


def phase_exponent_to_mod3_sign(exponent: int) -> int:
    """Image of zeta_6^exponent modulo 1-omega, encoded in F3."""
    return 1 if int(exponent) % 2 == 0 else 2


# ---------------------------------------------------------------------------
# F4 reduction modulo 2.  Elements are encoded as b0 + b1*alpha in bits,
# alpha^2+alpha+1=0: 0,1,alpha,alpha+1 -> 0,1,2,3.


def gf4_add(left: int, right: int) -> int:
    return int(left) ^ int(right)


def gf4_multiply(left: int, right: int) -> int:
    x, y = int(left), int(right)
    x0, x1 = x & 1, (x >> 1) & 1
    y0, y1 = y & 1, (y >> 1) & 1
    z0 = (x0 & y0) ^ (x1 & y1)
    z1 = (x0 & y1) ^ (x1 & y0) ^ (x1 & y1)
    return z0 | (z1 << 1)


def gf4_inverse(value: int) -> int:
    reduced = int(value)
    if reduced == 1:
        return 1
    if reduced == 2:
        return 3
    if reduced == 3:
        return 2
    raise ZeroDivisionError("zero has no inverse in F4")


def gf4_conjugate(value: int) -> int:
    return gf4_multiply(value, value)


F4_FACTORS: tuple[tuple[int, ...], ...] = (
    (1, 1),
    (1, 1, 0, 1),
    (1, 0, 1, 1),
    (1, 2, 1),
    (1, 3, 1),
    (1, 2, 1, 3, 3, 0, 1),
    (1, 3, 1, 2, 2, 0, 1),
    (1, 0, 2, 2, 1, 3, 1),
    (1, 0, 3, 3, 1, 2, 1),
)
F4_STAR_PAIRS: tuple[tuple[int, int], ...] = (
    (1, 2),
    (3, 4),
    (5, 7),
    (6, 8),
)


def gf4_star_polynomial(polynomial: Sequence[int]) -> tuple[int, ...]:
    raw = tuple(gf4_conjugate(value) for value in reversed(polynomial))
    scale = gf4_inverse(raw[-1])
    return _trim(tuple(gf4_multiply(scale, value) for value in raw))


def gf4_factorization_is_valid() -> bool:
    product: tuple[int, ...] = (1,)
    for factor in F4_FACTORS:
        product = _poly_mul(product, factor, gf4_add, gf4_multiply)
    target = (1,) + (0,) * 34 + (1,)
    star_map = tuple(F4_FACTORS.index(gf4_star_polynomial(f)) for f in F4_FACTORS)
    expected_star = (0, 2, 1, 4, 3, 7, 8, 5, 6)
    return (
        product == target
        and star_map == expected_star
        and all(
            _is_irreducible(factor, 4, gf4_add, gf4_multiply, gf4_inverse)
            for factor in F4_FACTORS
        )
    )


def gf4_code_generators() -> tuple[tuple[int, ...], ...]:
    generators: list[tuple[int, ...]] = []
    for choices in itertools.product((0, 1), repeat=len(F4_STAR_PAIRS)):
        generator = F4_FACTORS[0]
        for (left, right), choice in zip(F4_STAR_PAIRS, choices, strict=True):
            generator = _poly_mul(
                generator,
                F4_FACTORS[(left, right)[choice]],
                gf4_add,
                gf4_multiply,
            )
        if len(generator) != 19:
            raise AssertionError("F4 minimal generator must have degree 18")
        generators.append(generator)
    return tuple(generators)


def phase_exponent_to_gf4(exponent: int) -> int:
    """Image of zeta_6^exponent modulo 2, encoded in F4."""
    return (1, 3, 2)[int(exponent) % 3]


def _gf4_nullspace(matrix: Sequence[Sequence[int]], columns: int) -> tuple[tuple[int, ...], ...]:
    rows = [list(map(int, row)) for row in matrix]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        found = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column]),
            None,
        )
        if found is None:
            continue
        rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
        scale = gf4_inverse(rows[pivot_row][column])
        rows[pivot_row] = [gf4_multiply(scale, value) for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [
                gf4_add(value, gf4_multiply(factor, pivot_value))
                for value, pivot_value in zip(rows[row], rows[pivot_row], strict=True)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break

    free_columns = [column for column in range(columns) if column not in pivot_columns]
    basis: list[tuple[int, ...]] = []
    for free in free_columns:
        vector = [0] * columns
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            # Minus equals plus in characteristic two.
            vector[pivot] = rows[row][free]
        basis.append(tuple(vector))
    return tuple(basis)


def _gf4_codeword(message: Sequence[int], generator: Sequence[int]) -> tuple[int, ...]:
    result = [0] * Q_MODULUS
    for shift, coefficient in enumerate(message):
        if coefficient == 0:
            continue
        for index, value in enumerate(generator):
            result[shift + index] = gf4_add(
                result[shift + index], gf4_multiply(coefficient, value)
            )
    return tuple(result)


@dataclass(frozen=True)
class ShortenedCodeProfile:
    support: tuple[int, ...]
    generator_index: int
    dimension: int
    maximum_weight: int
    full_support_codeword_count: int


def shortened_gf4_code_profile(
    support: Iterable[int], generator: Sequence[int], generator_index: int
) -> ShortenedCodeProfile:
    support_t = tuple(sorted(int(q) for q in support))
    support_set = set(support_t)
    outside = [q for q in range(Q_MODULUS) if q not in support_set]
    message_length = Q_MODULUS - (len(generator) - 1)
    equations = []
    for coordinate in outside:
        row = []
        for shift in range(message_length):
            index = coordinate - shift
            row.append(generator[index] if 0 <= index < len(generator) else 0)
        equations.append(row)
    basis = _gf4_nullspace(equations, message_length)

    maximum_weight = 0
    full_count = 0
    for coefficients in itertools.product(range(4), repeat=len(basis)):
        if not any(coefficients):
            continue
        message = [0] * message_length
        for coefficient, vector in zip(coefficients, basis, strict=True):
            if coefficient:
                message = [
                    gf4_add(value, gf4_multiply(coefficient, basis_value))
                    for value, basis_value in zip(message, vector, strict=True)
                ]
        codeword = _gf4_codeword(message, generator)
        if any(codeword[q] for q in outside):
            raise AssertionError("shortened-code nullspace leaked outside the support")
        weight = sum(value != 0 for value in codeword)
        maximum_weight = max(maximum_weight, weight)
        if weight == len(support_t) and all(codeword[q] for q in support_t):
            full_count += 1
    return ShortenedCodeProfile(
        support=support_t,
        generator_index=int(generator_index),
        dimension=len(basis),
        maximum_weight=maximum_weight,
        full_support_codeword_count=full_count,
    )


@dataclass(frozen=True)
class CongruenceProofResult:
    mod3_candidates_per_code: tuple[int, int]
    mod3_support_count: int
    mod3_support_orbit_count: int
    shortened_profile_count: int
    shortened_dimension_counts: tuple[tuple[int, int], ...]
    shortened_maximum_weight_counts: tuple[tuple[int, int], ...]
    compatible_mod2_support_count: int
    profiles: tuple[ShortenedCodeProfile, ...]


@lru_cache(maxsize=1)
def congruence_nonexistence_proof() -> CongruenceProofResult:
    if not gf3_factorization_is_valid():
        raise AssertionError("F3 factorization or irreducibility check failed")
    if not gf4_factorization_is_valid():
        raise AssertionError("F4 factorization or irreducibility check failed")

    candidates_by_code = enumerate_mod3_candidates_numpy()
    for candidates in candidates_by_code:
        for candidate in candidates:
            if mod3_periodic_autocorrelation(candidate.sequence) != (0,) * Q_MODULUS:
                raise AssertionError("enumerated F3 candidate fails zero autocorrelation")

    supports = tuple(
        sorted({candidate.support for code in candidates_by_code for candidate in code})
    )
    canonical_supports = {canonical_support(support) for support in supports}
    generators = gf4_code_generators()
    profiles = tuple(
        shortened_gf4_code_profile(support, generator, generator_index)
        for support in supports
        for generator_index, generator in enumerate(generators)
    )
    compatible_supports = {
        profile.support for profile in profiles if profile.full_support_codeword_count
    }
    return CongruenceProofResult(
        mod3_candidates_per_code=tuple(len(code) for code in candidates_by_code),  # type: ignore[arg-type]
        mod3_support_count=len(supports),
        mod3_support_orbit_count=len(canonical_supports),
        shortened_profile_count=len(profiles),
        shortened_dimension_counts=tuple(sorted(Counter(p.dimension for p in profiles).items())),
        shortened_maximum_weight_counts=tuple(
            sorted(Counter(p.maximum_weight for p in profiles).items())
        ),
        compatible_mod2_support_count=len(compatible_supports),
        profiles=profiles,
    )


def congruence_result_summary(result: CongruenceProofResult) -> dict[str, object]:
    """Return the deterministic, compact machine-readable proof summary."""
    candidates_by_code = enumerate_mod3_candidates_numpy()
    supports_by_code = tuple(
        sorted({candidate.support_mask for candidate in code})
        for code in candidates_by_code
    )
    all_support_masks = sorted(set().union(*(set(values) for values in supports_by_code)))
    return {
        "statement": "No normalized sixth-root phase lift exists for CW(105,36).",
        "method": "finite-field congruence obstruction over F3 and F4",
        "f3_factorization_valid": gf3_factorization_is_valid(),
        "f4_factorization_valid": gf4_factorization_is_valid(),
        "mod3_candidates_per_code": list(result.mod3_candidates_per_code),
        "mod3_supports_per_code": [len(values) for values in supports_by_code],
        "mod3_support_count": result.mod3_support_count,
        "mod3_support_orbit_count": result.mod3_support_orbit_count,
        "shortened_profile_count": result.shortened_profile_count,
        "shortened_dimension_counts": [list(item) for item in result.shortened_dimension_counts],
        "shortened_maximum_weight_counts": [
            list(item) for item in result.shortened_maximum_weight_counts
        ],
        "compatible_mod2_support_count": result.compatible_mod2_support_count,
        "support_masks": all_support_masks,
        "conclusion": "nonexistence" if result.compatible_mod2_support_count == 0 else "inconclusive",
    }


def main() -> int:
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Run the exact F3/F4 congruence obstruction for CW(105,36)."
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSON path; the same summary is always printed to stdout",
    )
    args = parser.parse_args()
    result = congruence_nonexistence_proof()
    summary = congruence_result_summary(result)
    rendered = json.dumps(summary, sort_keys=True, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0 if result.compatible_mod2_support_count == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
