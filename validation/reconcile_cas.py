#!/usr/bin/env python3
"""Reconcile independently generated SageMath and GAP artifacts."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
import sys
from pathlib import Path
from typing import Any


HEX_MASK = re.compile(r"[0-9a-f]{9}")
DECIMAL_MASK = re.compile(r"(?:0|[1-9][0-9]*)")
LENGTH = 35
WEIGHT = 12


class ReconciliationError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReconciliationError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"{path}: top-level JSON must be an object")
    return value


def read_hex_supports(path: Path) -> tuple[bytes, tuple[int, ...]]:
    raw = path.read_bytes()
    require(raw.endswith(b"\n"), f"{path}: support file must end with LF")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as error:
        raise ReconciliationError(f"{path}: support file is not ASCII") from error
    require(lines, f"{path}: support file is empty")
    require(
        all(HEX_MASK.fullmatch(line) for line in lines),
        f"{path}: support masks must be nine-digit lowercase hexadecimal",
    )
    masks = tuple(int(line, 16) for line in lines)
    validate_masks(path, masks)
    return raw, masks


def read_decimal_certificate(path: Path) -> tuple[int, ...]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    require(lines, f"{path}: certificate is empty")
    require(
        all(DECIMAL_MASK.fullmatch(line) for line in lines),
        f"{path}: certificate masks must be decimal integers",
    )
    masks = tuple(int(line) for line in lines)
    validate_masks(path, masks)
    return masks


def validate_masks(path: Path, masks: tuple[int, ...]) -> None:
    require(masks == tuple(sorted(masks)), f"{path}: support masks are not sorted")
    require(len(masks) == len(set(masks)), f"{path}: duplicate support mask")
    require(
        all(0 <= mask < (1 << LENGTH) for mask in masks),
        f"{path}: support mask exceeds {LENGTH} bits",
    )
    require(
        all(mask.bit_count() == WEIGHT for mask in masks),
        f"{path}: support mask does not have weight {WEIGHT}",
    )


def factor_degrees(factors: list[dict[str, Any]]) -> list[int]:
    return [int(factor["degree"]) for factor in factors]


def branch_dimensions(branches: list[dict[str, Any]]) -> list[int]:
    return [int(branch["dimension"]) for branch in branches]


def validate_factor_orbits(
    label: str,
    factors: list[dict[str, Any]],
    orbits: list[list[int]],
    sharp_images: list[int],
) -> None:
    factor_indices = [int(factor["index"]) for factor in factors]
    require(
        len(factor_indices) == len(set(factor_indices)),
        f"{label}: duplicate factor index",
    )
    require(
        all(len(orbit) in (1, 2) for orbit in orbits),
        f"{label}: factor-sharp orbit is not singleton or paired",
    )
    flattened = [int(index) for orbit in orbits for index in orbit]
    require(
        len(flattened) == len(set(flattened)),
        f"{label}: factor-sharp orbits overlap",
    )
    require(
        set(flattened) == set(factor_indices),
        f"{label}: factor-sharp orbits do not partition the factorization",
    )
    require(
        len(sharp_images) == len(factor_indices),
        f"{label}: factor-sharp image list has the wrong length",
    )
    image_by_index = {
        factor_index: int(sharp_images[position])
        for position, factor_index in enumerate(factor_indices)
    }
    require(
        all(image in image_by_index for image in image_by_index.values()),
        f"{label}: factor-sharp image is not a factor index",
    )
    require(
        all(image_by_index[image_by_index[index]] == index for index in factor_indices),
        f"{label}: factor-sharp images are not an involution",
    )
    derived_orbits = {
        frozenset((index, image_by_index[index])) for index in factor_indices
    }
    reported_orbits = {
        frozenset(int(index) for index in orbit) for orbit in orbits
    }
    require(
        reported_orbits == derived_orbits,
        f"{label}: reported factor orbits disagree with sharp images",
    )


def validate_factor_branch_product(
    label: str,
    orbits: list[list[int]],
    branches: list[dict[str, Any]],
) -> None:
    expected = {
        tuple(sorted(choice))
        for choice in itertools.product(*(tuple(orbit) for orbit in orbits))
    }
    actual = {
        tuple(sorted(int(index) for index in branch["factor_indices"]))
        for branch in branches
    }
    require(
        len(actual) == len(branches),
        f"{label}: duplicate factor branch",
    )
    require(
        actual == expected,
        f"{label}: factor branches do not form the complete Cartesian product",
    )


def reconcile(args: argparse.Namespace) -> dict[str, Any]:
    sage_result = read_json(args.sage_result)
    gap_result = read_json(args.gap_result)
    sage_raw, sage_masks = read_hex_supports(args.sage_supports)
    gap_raw, gap_masks = read_hex_supports(args.gap_supports)
    certificate_masks = read_decimal_certificate(args.certificate)

    require(
        sage_raw == gap_raw,
        "SageMath and GAP support serializations differ",
    )
    require(
        sage_masks == certificate_masks,
        "generated support masks differ from the canonical certificate",
    )

    sage_char3 = sage_result["char3"]
    gap_char3 = gap_result["characteristic_3"]
    sage_char2 = sage_result["char2"]
    gap_char2 = gap_result["characteristic_2"]

    require(
        sage_result["status"] == gap_result["status"] == "nonexistent",
        "generated final statuses do not both establish nonexistence",
    )
    require(
        sage_result["external_expected_values_compared"] is False
        and gap_result["external_expected_values_compared"] is False,
        "a validator reports comparing expected values during its search",
    )
    require(
        all(sage_result["completeness"].values()),
        "SageMath completeness certificate is not entirely true",
    )
    require(
        all(gap_result["completeness_certificate"].values()),
        "GAP completeness certificate is not entirely true",
    )

    support_count = len(sage_masks)
    require(
        sage_char3["support_count"]
        == gap_char3["support_count"]
        == sage_char2["support_count_covered"]
        == gap_char2["eliminated_support_count"]
        == support_count,
        "generated support counts disagree",
    )
    require(
        sage_char3["support_masks_sha256"]
        == hashlib.sha256(sage_raw).hexdigest(),
        "SageMath result records the wrong support-list SHA-256",
    )

    sage_degrees3 = factor_degrees(sage_char3["factorization"])
    gap_degrees3 = factor_degrees(gap_char3["factorization"])
    require(
        sage_degrees3 == gap_degrees3,
        "characteristic-three factor degrees disagree",
    )
    sage_dimensions3 = branch_dimensions(sage_char3["branches"])
    gap_dimensions3 = branch_dimensions(gap_char3["branches"])
    require(
        sage_dimensions3 == gap_dimensions3,
        "characteristic-three branch dimensions disagree",
    )
    require(
        [branch["enumerated_words"] for branch in sage_char3["branches"]]
        == [branch["enumerated_words"] for branch in gap_char3["branches"]],
        "characteristic-three branch enumeration sizes disagree",
    )
    sage_orbits3 = [
        *[[int(index)] for index in sage_char3["fixed_factor_indices"]],
        *[
            [int(index) for index in orbit]
            for orbit in sage_char3["paired_factor_orbits"]
        ],
    ]
    validate_factor_orbits(
        "SageMath characteristic three",
        sage_char3["factorization"],
        sage_orbits3,
        [
            int(factor["reciprocal_index"])
            for factor in sage_char3["factorization"]
        ],
    )
    sage_fixed3 = [orbit[0] for orbit in sage_orbits3 if len(orbit) == 1]
    sage_full_branches3 = [
        {
            "factor_indices": [
                *sage_fixed3,
                *[int(index) for index in branch["selected_factor_indices"]],
            ]
        }
        for branch in sage_char3["branches"]
    ]
    validate_factor_branch_product(
        "SageMath characteristic three",
        sage_orbits3,
        sage_full_branches3,
    )
    validate_factor_orbits(
        "GAP characteristic three",
        gap_char3["factorization"],
        gap_char3["factor_sharp_orbits"],
        gap_char3["factor_sharp_images"],
    )
    validate_factor_branch_product(
        "GAP characteristic three",
        gap_char3["factor_sharp_orbits"],
        gap_char3["branches"],
    )

    require(
        sage_char3["affine_group"]["partition_complete"] is True
        and gap_char3["complete"] is True,
        "support or affine partition is incomplete",
    )
    require(
        len(sage_char3["affine_group"]["orbits"])
        == gap_char3["affine_orbit_count"],
        "affine orbit counts disagree",
    )
    require(
        [orbit["orbit_size"] for orbit in sage_char3["affine_group"]["orbits"]]
        == [orbit["size"] for orbit in gap_char3["affine_orbits"]],
        "affine orbit sizes disagree",
    )

    require(
        sage_char2["normalization"]["normalized_assignments_per_support"]
        == 3 ** (WEIGHT - 1),
        "SageMath normalized phase count is not 3^(weight-1)",
    )
    require(
        sage_char2["normalization"]["unnormalized_assignments_per_support"]
        == 3**WEIGHT,
        "SageMath unnormalized phase count is not 3^weight",
    )
    require(
        sage_char2["total_normalized_surviving_support_assignment_pairs"] == 0
        and gap_char2["raw_normalized_branch_occurrences"] == 0
        and gap_char2["unique_normalized_assignment_count"] == 0
        and gap_char2["surviving_support_count"] == 0,
        "a characteristic-two normalized assignment survived",
    )

    gap_branches4 = gap_char2["branches"]
    validate_factor_orbits(
        "GAP characteristic two",
        gap_char2["factorization"],
        gap_char2["factor_sharp_orbits"],
        gap_char2["factor_sharp_images"],
    )
    validate_factor_branch_product(
        "GAP characteristic two",
        gap_char2["factor_sharp_orbits"],
        gap_branches4,
    )
    expected_cases = support_count * len(gap_branches4)
    require(
        gap_char2["support_branch_cases_expected"]
        == gap_char2["support_branch_cases_processed"]
        == expected_cases,
        "GAP did not process the complete support-by-branch Cartesian product",
    )
    histogram = gap_char2["kernel_dimension_histogram_d0_through_d17"]
    require(
        sum(histogram) == expected_cases,
        "GAP kernel-dimension histogram does not cover every case",
    )
    require(
        sum(count * (4**dimension - 1) for dimension, count in enumerate(histogram))
        == gap_char2["nonzero_kernel_vectors_examined"],
        "GAP nonzero shortened-kernel vector count is inconsistent",
    )
    gap_factor_coefficients4 = [
        int(coefficient)
        for factor in gap_char2["factorization"]
        for coefficient in factor["coefficients"]
    ]
    require(
        all(coefficient in (0, 1, 2, 3) for coefficient in gap_factor_coefficients4)
        and any(coefficient in (2, 3) for coefficient in gap_factor_coefficients4),
        "GAP GF(4) factor record contains no coefficient outside the GF(2) subfield",
    )
    require(
        sage_char3["freeze_complete"] is True
        and sage_char2["freeze_complete"] is True
        and gap_char3["complete"] is True
        and gap_char2["complete"] is True
        and gap_char2["support_branch_case_completeness_checked"] is True,
        "a generated stage-level completeness flag is false",
    )

    checks = {
        "all_assertions_passed": True,
        "canonical_certificate_equal_as_exact_integer_masks": True,
        "complete_factor_branch_products_checked": True,
        "complete_support_branch_cartesian_product_checked": True,
        "explicit_gap_gf4_coefficients_checked": True,
        "factor_orbit_partitions_and_involutions_checked": True,
        "generated_results_structurally_reconciled": True,
        "no_expected_values_compared_during_search": True,
        "normalized_phase_count_checked": True,
        "stage_and_validator_completeness_attestations_true": True,
        "support_lists_byte_identical": True,
        "support_serialization_syntax_weight_and_integer_equality_checked": True,
        "zero_characteristic_two_survivors": True,
    }
    return {
        "schema_version": 1,
        "status": sage_result["status"],
        "declared_bit_convention": "bit q represents coordinate q in Z/35Z",
        "supports": {
            "count": support_count,
            "mask_encoding": "9-digit-lowercase-hex",
            "sha256": hashlib.sha256(sage_raw).hexdigest(),
        },
        "characteristic_3": {
            "factor_degrees": sage_degrees3,
            "branch_dimensions": sage_dimensions3,
            "branch_count": len(sage_dimensions3),
        },
        "characteristic_2": {
            "gap_factor_degrees": factor_degrees(gap_char2["factorization"]),
            "gap_branch_count": len(gap_branches4),
            "support_branch_cases": expected_cases,
            "surviving_normalized_assignments": 0,
        },
        "input_sha256": {
            "sage_result": file_sha256(args.sage_result),
            "sage_supports": file_sha256(args.sage_supports),
            "gap_result": file_sha256(args.gap_result),
            "gap_supports": file_sha256(args.gap_supports),
            "canonical_certificate": file_sha256(args.certificate),
        },
        "checks": checks,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sage-result", type=Path, required=True)
    parser.add_argument("--sage-supports", type=Path, required=True)
    parser.add_argument("--gap-result", type=Path, required=True)
    parser.add_argument("--gap-supports", type=Path, required=True)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = reconcile(args)
    except (KeyError, OSError, ReconciliationError, TypeError, json.JSONDecodeError) as error:
        print(f"reconciliation failed: {error}", file=sys.stderr)
        return 1
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized, encoding="utf-8", newline="\n")
    print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
