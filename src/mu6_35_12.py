"""Nonexistence of a 12-sparse perfect sixth-root sequence of length 35.

A sequence ``u`` with entries in ``{0} union mu_6`` is *perfect* here when
its periodic Hermitian autocorrelation equals 12 at shift zero and zero at
every nonzero shift.  Such a sequence would generate a circulant complex
generalized weighing matrix ``CGW(35,12;6)``.

The proof has two parts.

1. Reduction modulo ``1-omega`` places the support in one of two ternary
   cyclic [35,12] codes.  Exact enumeration shows that their weight-12
   words have 420 distinct supports, and those supports form one affine
   orbit under ``q -> a*q+b`` on Z_35.
2. For one canonical support, reduction modulo 2 gives an elementary
   contradiction from only shifts 2, 4, 8 and 10 over F_4.

No condition inherited from the CW(105,36) contraction is used.  In
particular, the four positions {0,15,25,30} are *not* fixed to zero.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, Mapping, Sequence

import numpy as np

from .eisenstein_congruences import (
    Mod3Candidate,
    _all_ternary_messages,
    gf3_code_generators,
    gf3_factorization_is_valid,
    gf4_add,
    gf4_inverse,
    gf4_multiply,
    mod3_periodic_autocorrelation,
)

N = 35
WEIGHT = 12
AFFINE_UNITS = tuple(a for a in range(N) if math.gcd(a, N) == 1)
CANONICAL_SUPPORT = (0, 1, 2, 3, 7, 10, 12, 16, 21, 22, 26, 28)
CORE_SHIFTS = (2, 4, 8, 10)


def support_mask(support: Iterable[int]) -> int:
    """Encode a subset of Z_35 as a 35-bit integer."""
    mask = 0
    for q in support:
        value = int(q) % N
        mask |= 1 << value
    return mask


def support_from_mask(mask: int) -> tuple[int, ...]:
    """Decode a 35-bit support mask."""
    value = int(mask)
    if value < 0 or value >> N:
        raise ValueError("support mask is not contained in Z_35")
    return tuple(q for q in range(N) if (value >> q) & 1)


def affine_transform_support(
    support: Iterable[int], multiplier: int, translation: int
) -> tuple[int, ...]:
    """Apply q -> multiplier*q+translation on Z_35 to a support."""
    a = int(multiplier) % N
    b = int(translation) % N
    if math.gcd(a, N) != 1:
        raise ValueError("affine multiplier must be a unit modulo 35")
    return tuple(sorted({(a * int(q) + b) % N for q in support}))


def affine_support_orbit(support: Iterable[int]) -> tuple[tuple[int, ...], ...]:
    """Return the complete AGL(1,Z_35) orbit of a support."""
    support_t = tuple(sorted(int(q) % N for q in support))
    return tuple(
        sorted(
            {
                affine_transform_support(support_t, a, b)
                for a in AFFINE_UNITS
                for b in range(N)
            }
        )
    )


def affine_stabilizer(support: Iterable[int]) -> tuple[tuple[int, int], ...]:
    """Return all affine maps fixing the support setwise."""
    support_t = tuple(sorted(int(q) % N for q in support))
    return tuple(
        (a, b)
        for a in AFFINE_UNITS
        for b in range(N)
        if affine_transform_support(support_t, a, b) == support_t
    )


@lru_cache(maxsize=1)
def enumerate_all_mod3_weight12_candidates_numpy() -> tuple[tuple[Mod3Candidate, ...], ...]:
    """Enumerate all weight-12 words in the two forced ternary cyclic codes.

    This is the generalized enumeration.  Unlike
    ``enumerate_mod3_candidates_numpy`` in ``src.eisenstein_congruences``, it
    imposes no prescribed zero coordinates.
    """
    messages = _all_ternary_messages(12)
    result: list[tuple[Mod3Candidate, ...]] = []
    for code_index, generator in enumerate(gf3_code_generators()):
        matrix = np.zeros((12, N), dtype=np.int16)
        generator_array = np.asarray(generator, dtype=np.int16)
        for shift in range(12):
            matrix[shift, shift : shift + len(generator)] = generator_array
        codewords = (messages @ matrix) % 3
        selected = np.flatnonzero(np.count_nonzero(codewords, axis=1) == WEIGHT)
        result.append(
            tuple(
                Mod3Candidate(
                    code_index=code_index,
                    message=tuple(int(value) for value in messages[index]),
                    sequence=tuple(int(value) for value in codewords[index]),
                )
                for index in selected
            )
        )
    return tuple(result)


def all_mod3_weight12_supports() -> tuple[tuple[int, ...], ...]:
    """Return all supports allowed by the characteristic-three reduction."""
    candidates = enumerate_all_mod3_weight12_candidates_numpy()
    return tuple(sorted({candidate.support for code in candidates for candidate in code}))


def decimal_mask_lines_sha256(masks: Iterable[int]) -> str:
    """Hash the canonical decimal-lines serialization of support masks."""
    payload = "".join(f"{int(mask)}\n" for mask in sorted(int(m) for m in masks))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def gf4_ratio(numerator: int, denominator: int) -> int:
    """Return numerator/denominator in F_4, requiring denominator nonzero."""
    if int(denominator) == 0:
        raise ZeroDivisionError("F4 ratio by zero")
    return gf4_multiply(int(numerator), gf4_inverse(int(denominator)))


def gf4_correlation_on_support(
    values: Mapping[int, int], shift: int, support: Sequence[int] = CANONICAL_SUPPORT
) -> int:
    """Compute one exact periodic Hermitian correlation in F_4.

    ``values`` must assign a nonzero F_4 element (1, 2 or 3) to every support
    coordinate.  Conjugation on F_4 is inversion on F_4^*.
    """
    support_set = set(int(q) for q in support)
    result = 0
    for q in support:
        target = (int(q) + int(shift)) % N
        if target not in support_set:
            continue
        left = int(values[int(q)])
        right = int(values[target])
        if left not in (1, 2, 3) or right not in (1, 2, 3):
            raise ValueError("active F4 values must be nonzero")
        result = gf4_add(result, gf4_ratio(left, right))
    return result


@dataclass(frozen=True)
class FourShiftCertificate:
    """The exact four-shift identity for one F_4 phase assignment."""

    c2: int
    c4: int
    c8: int
    c10: int
    p: int | None
    scaled_c2: int | None
    four_terms_of_c10: int | None
    residual_two_terms: int | None
    forced_c4_value: int | None
    contradiction_triggered: bool


def four_shift_certificate(values: Mapping[int, int]) -> FourShiftCertificate:
    """Evaluate the hand-checkable F_4 contradiction on the canonical support.

    If C_8=0, put

        p = x_2/x_10 = x_28/x_1.

    Then p*C_2 is exactly the sum of four of the six terms in C_10.  Hence
    C_2=C_8=C_10=0 forces

        x_12/x_22 = x_16/x_26,

    equivalently x_12/x_16 = x_22/x_26.  Those two equal terms cancel in
    C_4, leaving x_3/x_7, which is nonzero.  Thus C_4 cannot also vanish.
    """
    required = set(CANONICAL_SUPPORT)
    if set(values) != required:
        missing = sorted(required - set(values))
        extra = sorted(set(values) - required)
        raise ValueError(f"values must be indexed exactly by the canonical support; missing={missing}, extra={extra}")
    if any(int(values[q]) not in (1, 2, 3) for q in CANONICAL_SUPPORT):
        raise ValueError("all canonical-support values must lie in F4^*")

    c2 = gf4_correlation_on_support(values, 2)
    c4 = gf4_correlation_on_support(values, 4)
    c8 = gf4_correlation_on_support(values, 8)
    c10 = gf4_correlation_on_support(values, 10)
    if c8 != 0:
        return FourShiftCertificate(c2, c4, c8, c10, None, None, None, None, None, False)

    p = gf4_ratio(values[2], values[10])
    if p != gf4_ratio(values[28], values[1]):
        raise AssertionError("C8=0 did not force the two nonzero terms to agree")

    scaled_c2 = gf4_multiply(p, c2)
    four_terms = 0
    for left, right in ((0, 10), (28, 3), (2, 12), (26, 1)):
        four_terms = gf4_add(four_terms, gf4_ratio(values[left], values[right]))
    if scaled_c2 != four_terms:
        raise AssertionError("p*C2 identity failed")

    residual = gf4_add(
        gf4_ratio(values[12], values[22]),
        gf4_ratio(values[16], values[26]),
    )
    if c10 != gf4_add(four_terms, residual):
        raise AssertionError("C10 decomposition failed")

    forced_c4_value: int | None = None
    contradiction = False
    if c2 == 0 and c10 == 0:
        if residual != 0:
            raise AssertionError("C2=C8=C10=0 did not annihilate the residual terms")
        left_ratio = gf4_ratio(values[12], values[16])
        right_ratio = gf4_ratio(values[22], values[26])
        if left_ratio != right_ratio:
            raise AssertionError("cross-ratio consequence failed")
        forced_c4_value = gf4_ratio(values[3], values[7])
        if c4 != forced_c4_value or forced_c4_value == 0:
            raise AssertionError("forced nonzero C4 identity failed")
        contradiction = True

    return FourShiftCertificate(
        c2=c2,
        c4=c4,
        c8=c8,
        c10=c10,
        p=p,
        scaled_c2=scaled_c2,
        four_terms_of_c10=four_terms,
        residual_two_terms=residual,
        forced_c4_value=forced_c4_value,
        contradiction_triggered=contradiction,
    )


@dataclass(frozen=True)
class CoreEnumerationResult:
    normalized_assignments: int
    c8_zero: int
    c2_c8_zero: int
    c2_c8_c10_zero: int
    all_four_zero: int
    algebraic_identity_failures: int


def exhaustive_four_shift_core_check() -> CoreEnumerationResult:
    """Exhaustively regression-check the hand proof on 3^11 assignments."""
    values_order = (1, 2, 3)
    support = CANONICAL_SUPPORT
    counts = Counter()
    failures = 0
    for free_values in itertools.product(values_order, repeat=WEIGHT - 1):
        assignment = {support[0]: 1}
        assignment.update(zip(support[1:], free_values, strict=True))
        try:
            certificate = four_shift_certificate(assignment)
        except AssertionError:
            failures += 1
            continue
        counts["total"] += 1
        if certificate.c8 == 0:
            counts["c8"] += 1
        if certificate.c2 == 0 and certificate.c8 == 0:
            counts["c2c8"] += 1
        if certificate.c2 == 0 and certificate.c8 == 0 and certificate.c10 == 0:
            counts["c2c8c10"] += 1
            if not certificate.contradiction_triggered:
                failures += 1
        if all(
            value == 0
            for value in (certificate.c2, certificate.c4, certificate.c8, certificate.c10)
        ):
            counts["all"] += 1
    return CoreEnumerationResult(
        normalized_assignments=counts["total"],
        c8_zero=counts["c8"],
        c2_c8_zero=counts["c2c8"],
        c2_c8_c10_zero=counts["c2c8c10"],
        all_four_zero=counts["all"],
        algebraic_identity_failures=failures,
    )


@dataclass(frozen=True)
class Mu6NonexistenceResult:
    mod3_candidates_per_code: tuple[int, int]
    mod3_supports_per_code: tuple[int, int]
    mod3_support_intersection_size: int
    mod3_support_count: int
    affine_orbit_count: int
    affine_orbit_size: int
    affine_stabilizer: tuple[tuple[int, int], ...]
    canonical_support: tuple[int, ...]
    support_masks_sha256: str
    core: CoreEnumerationResult


@lru_cache(maxsize=1)
def mu6_35_12_nonexistence_proof() -> Mu6NonexistenceResult:
    """Run the complete generalized finite proof."""
    if not gf3_factorization_is_valid():
        raise AssertionError("F3 factorization or irreducibility check failed")

    candidates_by_code = enumerate_all_mod3_weight12_candidates_numpy()
    for candidates in candidates_by_code:
        for candidate in candidates:
            if mod3_periodic_autocorrelation(candidate.sequence) != (0,) * N:
                raise AssertionError("weight-12 F3 codeword failed direct autocorrelation")

    supports_by_code = tuple(
        {candidate.support for candidate in code} for code in candidates_by_code
    )
    intersection = supports_by_code[0] & supports_by_code[1]
    supports = supports_by_code[0] | supports_by_code[1]
    orbit = set(affine_support_orbit(CANONICAL_SUPPORT))
    if supports != orbit:
        missing = sorted(supports - orbit)
        extra = sorted(orbit - supports)
        raise AssertionError(
            f"F3 support set is not the claimed affine orbit: missing={len(missing)}, extra={len(extra)}"
        )
    if min(orbit) != CANONICAL_SUPPORT:
        raise AssertionError("stored support is not affine-canonical")

    core = exhaustive_four_shift_core_check()
    if core.normalized_assignments != 3**11:
        raise AssertionError("unexpected normalized F4 assignment count")
    if core.all_four_zero != 0 or core.algebraic_identity_failures != 0:
        raise AssertionError("four-shift F4 contradiction failed")

    masks = tuple(sorted(support_mask(support) for support in supports))
    return Mu6NonexistenceResult(
        mod3_candidates_per_code=tuple(len(code) for code in candidates_by_code),
        mod3_supports_per_code=tuple(len(values) for values in supports_by_code),
        mod3_support_intersection_size=len(intersection),
        mod3_support_count=len(supports),
        affine_orbit_count=1,
        affine_orbit_size=len(orbit),
        affine_stabilizer=affine_stabilizer(CANONICAL_SUPPORT),
        canonical_support=CANONICAL_SUPPORT,
        support_masks_sha256=decimal_mask_lines_sha256(masks),
        core=core,
    )


def mu6_result_summary(result: Mu6NonexistenceResult) -> dict[str, object]:
    """Return a deterministic machine-readable theorem summary."""
    supports = all_mod3_weight12_supports()
    masks = tuple(sorted(support_mask(support) for support in supports))
    return {
        "statement": "No 12-sparse perfect sixth-root sequence of length 35 exists.",
        "matrix_form": "No circulant CGW(35,12;6) exists.",
        "uses_prescribed_zero_positions": False,
        "f3_factorization_valid": gf3_factorization_is_valid(),
        "mod3_candidates_per_code": list(result.mod3_candidates_per_code),
        "mod3_supports_per_code": list(result.mod3_supports_per_code),
        "mod3_support_intersection_size": result.mod3_support_intersection_size,
        "mod3_support_count": result.mod3_support_count,
        "affine_orbit_count": result.affine_orbit_count,
        "affine_orbit_size": result.affine_orbit_size,
        "affine_stabilizer": [list(item) for item in result.affine_stabilizer],
        "canonical_support": list(result.canonical_support),
        "core_shifts": list(CORE_SHIFTS),
        "normalized_f4_assignments_regression_checked": result.core.normalized_assignments,
        "core_counts": {
            "c8_zero": result.core.c8_zero,
            "c2_c8_zero": result.core.c2_c8_zero,
            "c2_c8_c10_zero": result.core.c2_c8_c10_zero,
            "all_four_zero": result.core.all_four_zero,
            "algebraic_identity_failures": result.core.algebraic_identity_failures,
        },
        "support_masks_sha256_decimal_lines": result.support_masks_sha256,
        "support_masks": list(masks),
        "conclusion": "nonexistence",
    }


def main() -> int:
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Prove nonexistence of a 12-sparse perfect mu_6 sequence of length 35."
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--masks-output", type=Path)
    args = parser.parse_args()

    result = mu6_35_12_nonexistence_proof()
    summary = mu6_result_summary(result)
    rendered = json.dumps(summary, sort_keys=True, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    if args.masks_output is not None:
        masks = summary["support_masks"]
        assert isinstance(masks, list)
        args.masks_output.parent.mkdir(parents=True, exist_ok=True)
        args.masks_output.write_text(
            "# all F3-admissible weight-12 supports in Z_35; decimal masks\n"
            + "".join(f"{int(mask)}\n" for mask in masks),
            encoding="utf-8",
        )
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
