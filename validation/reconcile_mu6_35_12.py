#!/usr/bin/env python3
"""Reconcile the generalized length-35 theorem outputs exactly."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return data


def read_masks(path: Path) -> list[int]:
    values = [
        int(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    if values != sorted(values) or len(values) != len(set(values)):
        raise ValueError("mask certificate must be sorted and duplicate-free")
    return values


def decimal_lines_hash(masks: list[int]) -> str:
    payload = "".join(f"{mask}\n" for mask in masks).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--direct", type=Path, required=True)
    parser.add_argument("--masks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    py = read_json(args.python)
    cpp = read_json(args.cpp)
    direct = read_json(args.direct)
    masks = read_masks(args.masks)

    expected_common = {
        "uses_prescribed_zero_positions": False,
        "f3_factorization_valid": True,
        "mod3_candidates_per_code": [420, 420],
        "mod3_supports_per_code": [210, 210],
        "mod3_support_intersection_size": 0,
        "mod3_support_count": 420,
        "affine_orbit_size": 420,
        "affine_stabilizer": [[1, 0], [29, 28]],
        "canonical_support": [0, 1, 2, 3, 7, 10, 12, 16, 21, 22, 26, 28],
        "core_shifts": [2, 4, 8, 10],
        "core_counts": {
            "c8_zero": 59049,
            "c2_c8_zero": 15309,
            "c2_c8_c10_zero": 5103,
            "all_four_zero": 0,
            "algebraic_identity_failures": 0,
        },
    }
    for key, value in expected_common.items():
        if py.get(key) != value:
            raise AssertionError(f"Python mismatch for {key}: {py.get(key)!r}")
        if cpp.get(key) != value:
            raise AssertionError(f"C++ mismatch for {key}: {cpp.get(key)!r}")

    if py.get("support_masks") != masks:
        raise AssertionError("Python support list differs from certificate")
    if cpp.get("support_masks") != masks:
        raise AssertionError("C++ support list differs from certificate")

    support_hash = decimal_lines_hash(masks)
    if py.get("support_masks_sha256_decimal_lines") != support_hash:
        raise AssertionError("Python support hash differs from certificate hash")
    if support_hash != "c067600256b37e077ee8a83889ec5e09347c397feea723466fb863ca2685f3d2":
        raise AssertionError("unexpected generalized support certificate hash")

    expected_direct = {
        "support_count": 420,
        "uses_prescribed_zero_positions": False,
        "global_phase_fixed": True,
        "assignments_per_support": 177147,
        "assignments_checked": 74401740,
        "compatible_support_count": 0,
        "zero_correlation_assignments_fixed_phase": 0,
    }
    for key, value in expected_direct.items():
        if direct.get(key) != value:
            raise AssertionError(f"direct checker mismatch for {key}: {direct.get(key)!r}")
    histogram = direct.get("first_failure_shift_counts")
    if not isinstance(histogram, list) or len(histogram) != 17:
        raise AssertionError("direct checker failure histogram must have 17 entries")
    if sum(int(value) for value in histogram) != expected_direct["assignments_checked"]:
        raise AssertionError("direct checker failure histogram does not cover all assignments")

    result = {
        "statement": "No 12-sparse perfect sixth-root sequence of length 35 exists.",
        "reconciled": True,
        "support_count": 420,
        "support_masks_sha256_decimal_lines": support_hash,
        "python_cpp_support_lists_equal": True,
        "four_shift_all_zero_assignments": 0,
        "direct_assignments_checked": 74401740,
        "direct_compatible_assignments": 0,
        "uses_prescribed_zero_positions": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
